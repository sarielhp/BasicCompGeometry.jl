# Performance Remediation Plan — Review 010

- **Source Report**: [`reviews/010_performance.md`](reviews/010_performance.md)
- **Scope**: `src/IpeDraw.jl`, `src/IpeDraw/interface.jl`, `src/IpeDraw/primitives.jl`, `src/IpeDraw/diagrams.jl`, `test/test_ipedraw.jl`, `CHANGELOG.md`

## Triage Summary

| # | Location | Severity | Verdict | Summary |
|---|---|---|---|---|
| 1 | `src/IpeDraw.jl:154-159` (`_apply_tf`) | Moderate | **Genuine** | Value-dependent return type (`Union{Point{2,T}, Point{2,Float64}}`) when `viewport === nothing` and input is non-`Float64`. |
| 2 | `src/IpeDraw.jl:50-52, 107`; `src/IpeDraw/interface.jl:14-16` (`Style`, `IpeCanvas.active_style`, `_style_options`) | Moderate | **Genuine** | Abstract `Style.values::NamedTuple`, abstract `active_style` read without a function barrier, and generator-based `NamedTuple` filtering cause dynamic dispatch, boxing, and per-element allocations. |
| 3 | `src/IpeDraw/primitives.jl:125-149, 364-399, 420-444` (+ adjacent `452-484`) | Moderate | **Genuine** | Path serialization allocates intermediate transformed-point vectors, slices, per-segment strings, and `Vector{String}` buffers before `join`. |
| 4 | `src/IpeDraw/diagrams.jl:5-10, 32-45` (`_escape_ipe_xml`, `draw_label!`) | Moderate | **Genuine** | Per-call 9-entry `Dict` allocation for `anchor`, three sequential `replace` passes in `_escape_ipe_xml`, and intermediate string concatenations in `draw_label!`. |

All 4 findings are genuine performance defects. `src/IpeDraw/output.jl` is in scope but has no findings in `reviews/010_performance.md` and requires no changes.

---

## Implementation Plan

### 1. Type-Stable `_apply_tf` (`src/IpeDraw.jl`)
- Update `_apply_tf(canvas::IpeCanvas, p::Point{2})` to convert `x = Float64(p.x)` and `y = Float64(p.y)` unconditionally before checking `vp === nothing`, returning `Point(x, y)` in the un-transformed branch and `Point(vp.tx + vp.scale * x, vp.ty + y_scale * y)` otherwise.
- Return type is concretely `Point{2,Float64}` for all `Point{2}` coordinate types (`Int`, `Float32`, `Float64`, `Rational`, etc.).

### 2. Concrete `Style{NT}`, Compile-Time Option Filtering, and Function Barrier (`src/IpeDraw.jl`, `src/IpeDraw/interface.jl`, `src/IpeDraw/diagrams.jl`)
- Parameterize `Style`:
  ```julia
  struct Style{NT<:NamedTuple}
      values::NT
  end
  Style(values::NT) where {NT<:NamedTuple} = Style{NT}(values)
  Style(; kwargs...) = Style((; kwargs...))
  ```
- Keep `IpeCanvas.active_style::Style` so `IpeCanvas` remains unparameterized and mutable across `with_style(f, canvas, ...)` scopes without breaking user callbacks typed on `canvas::IpeCanvas`.
- Replace generator-based `NamedTuple` construction in `_style_options` and `_draw_fill_opacity!` with a shared `Base.@constprop :aggressive @inline` tuple key filter `_filter_keys(ks, allowed, ::Val{keep})` and `NamedTuple{ks}(...)` construction.
- Update `_style_options` to accept concrete `active::Style{A}`, `style::S where {S<:Union{Nothing,Style}}`, `overrides::O where {O<:NamedTuple}`, and `::Val{allowed}` (avoiding unbound type parameters for `Aqua.test_unbound_args` compatibility).
- Introduce a function-barrier dispatch layer (`_dispatch_active`, `_draw_styled!`, `_draw_one!`, `_draw_each!`, `_label_styled!`) so `draw!`, `mark!`, and `label!` pass `canvas.active_style` across the barrier once per call (or once per vector in `draw!(canvas, objects::AbstractVector; ...)`), with a type-narrowing fast path when `active isa Style{@NamedTuple{}}`.
- Specialize `_draw_fill_opacity!(renderer::F, ...)` on `F` and extract `_without_key` / `_filter_keys(..., Val(false))` to eliminate generator allocations.

### 3. Single-Pass `IOBuffer` Streaming for Path Primitives (`src/IpeDraw/primitives.jl`)
- Add shared streaming helpers:
  - `_write_attr!(buf::IO, name::String, val::Union{Symbol,String,Nothing})` to write formatted XML attributes directly to `buf` without intermediate strings.
  - `_write_pt!(buf::IO, p::Point{2})` using `@printf(buf, "%.3f %.3f", p.x, p.y)`.
- Refactor `draw_polygon!` to iterate `poly` directly into a single `IOBuffer` (handling `< 2`-point iterables cleanly without allocating when fewer than 2 points are present).
- Refactor `draw_spline!` to stream `<path ...>`, `m`, `c` segments, optional `h`, and `</path>` directly into a single `IOBuffer`.
- Refactor `draw_bspline!` to validate `n >= 3 || throw(ArgumentError(...))` and iterate `pts` directly into a single `IOBuffer` without `tf_pts`, `tf_pts[2:end]`, `pt_strs`, or `lines` arrays (supporting any iterable without assuming `nextind`).
- Update adjacent `draw_polygon_with_holes!` / `_append_subpath!` in the same file to stream directly into `IOBuffer` as well, keeping path serialization uniform.

### 4. Constant Alignment Lookup and Streaming XML Escaping (`src/IpeDraw/diagrams.jl`)
- Move the 9-entry anchor table to `const _LABEL_ALIGNMENTS = Dict{Symbol,Tuple{Symbol,Symbol}}(...)` and extract `_anchor_alignment(anchor::Symbol)`.
- Implement `_write_escaped_ipe_xml!(buf::IO, s::AbstractString)` to escape `&`, `<`, `>` in a single pass directly into `buf` (with a fast `write(buf, s)` path when none of `&`, `<`, `>` appear), and delegate `_escape_ipe_xml(s)` to it so `ipe_group` continues to work identically.
- Update `draw_label!` to stream the `<text ...>...</text>` XML directly into a single `IOBuffer` (writing `$` math wrappers directly to `buf` without intermediate string concatenation).

### 5. Verification & Documentation
- Add regression tests in `test/test_ipedraw.jl` covering:
  - Type stability of `_apply_tf` for `Int` / `Float32` / `Rational` points with and without a viewport.
  - Parameterized `Style{NT}`, inferred `_style_options`, `with_style`, `mark!`, `label!`, `fill_opacity`, `draw!(canvas, vector)`, and low-allocation repeated `draw!`.
  - Exact output, edge cases (`< 2` polygon points, `< 3` B-spline `ArgumentError`, open/closed splines and B-splines, polygons with holes), and bounded allocations on large path inputs.
  - All 9 label anchors, invalid anchor `ArgumentError`, XML entity escaping (`&`, `<`, `>`), and low allocations for repeated anchored labels.
- Update `CHANGELOG.md` under `## [Unreleased]`.
- Run `tools/gate` (and `tools/build`) to confirm all tests and `Aqua` checks pass.
