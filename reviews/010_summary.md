# Post-Work Summary — Review 010 (Performance)

- **Source Report**: [`reviews/010_performance.md`](reviews/010_performance.md)
- **Remediation Plan**: [`reviews/010_performance_plan.md`](reviews/010_performance_plan.md)
- **Status**: Complete — all 4 genuine findings remediated, tested, and verified via `tools/build` and `tools/gate`.

---

## Issue-by-Issue Summary

### 1. Type Instability in `_apply_tf` (`src/IpeDraw.jl`)
- **Status**: Fixed (Genuine defect).
- **Mitigation**: Updated `_apply_tf(canvas::IpeCanvas, p::Point{2})` to convert `x = Float64(p.x)` and `y = Float64(p.y)` unconditionally before checking `canvas.viewport === nothing`, returning `Point(x, y)` (`Point{2,Float64}`) in both branches.
- **Verification**: Added `@inferred` checks in [`test/test_ipedraw.jl`](../test/test_ipedraw.jl) for `Int`, `Float32`, and `Rational{Int}` points both before and after `fit!`. Verified passing in `tools/gate`.

### 2. Abstract `Style.values`, Abstract `active_style` Read, and Dynamic Option Filtering (`src/IpeDraw.jl`, `src/IpeDraw/interface.jl`, `src/IpeDraw/diagrams.jl`)
- **Status**: Fixed (Genuine defect).
- **Mitigation**:
  - Parameterized `struct Style{NT<:NamedTuple}; values::NT; end` with `Style(; kwargs...) = Style((; kwargs...))`.
  - Kept `IpeCanvas.active_style::Style` so `IpeCanvas` remains unparameterized across `with_style` callbacks while routing `draw!`, `mark!`, and `label!` through concrete function barriers (`_dispatch_active`, `_draw_spec`, `_draw_one!`, `_draw_each!`, `_draw_styled!`, `_mark_many!`, `_label_one!`) with a type-narrowing fast path for the default `Style{@NamedTuple{}}`.
  - Replaced generator-based `NamedTuple` filtering in `_style_options` and `_draw_fill_opacity!` with a shared `Base.@constprop :aggressive @inline` tuple key filter `_filter_keys` and `NamedTuple{ks}(...)` construction, and typed `_draw_fill_opacity!(renderer::F, ...)` to specialize on the renderer function.
- **Verification**: Added regression tests in [`test/test_ipedraw.jl`](../test/test_ipedraw.jl) verifying `Style{NT}` concrete typing, `@inferred` `_style_options`, `ArgumentError` on unsupported options, `with_style` / `fill_opacity` behavior, and `draw!(canvas, vector)`. Aqua `Method ambiguity` and `Unbound type parameters` checks pass.

### 3. Path Serialization Intermediate Allocations (`src/IpeDraw/primitives.jl`, `src/IpeDraw.jl`)
- **Status**: Fixed (Genuine defect).
- **Mitigation**:
  - Added direct-to-buffer helpers `_write_uint!`, `_write_fixed3!`, `_write_pt!`, `_write_attr!`, and `_write_path_attrs!`, avoiding intermediate formatted strings and per-coordinate `@printf` buffer allocations.
  - Refactored `draw_polygon!`, `draw_spline!`, `draw_bspline!`, and adjacent `draw_polygon_with_holes!` / `_append_subpath!` to stream directly into a single pre-sized `IOBuffer`.
  - Updated `draw_bspline!` to iterate arbitrary iterables directly without allocating `tf_pts`, `tf_pts[2:end]`, `pt_strs`, or `lines`, and to throw `ArgumentError` for fewer than 3 control points.
- **Verification**: Measured 1,000-vertex `draw_polygon!` allocations drop from 789,936 B to 24,416 B, 999-segment `draw_spline!` from 2,205,024 B to 64,320 B, and 1,000-point `draw_bspline!` from 789,840 B to 24,384 B. Added functional and `< 3 * sizeof(xml)` allocation-bound regression tests in [`test/test_ipedraw.jl`](../test/test_ipedraw.jl).

### 4. Per-Label Alignment `Dict` and Multi-Pass XML Escaping (`src/IpeDraw/diagrams.jl`)
- **Status**: Fixed (Genuine defect).
- **Mitigation**:
  - Extracted the 9-entry anchor mapping to module-level `const _LABEL_ALIGNMENTS = Dict{Symbol,Tuple{Symbol,Symbol}}(...)` with helper `_anchor_alignment(anchor::Symbol)`.
  - Replaced three sequential `replace` passes with single-pass `_write_escaped_ipe_xml!(buf, s)` (using a fast `print(buf, s)` path when none of `&`, `<`, `>` appear), delegated `_escape_ipe_xml` to it, and streamed `draw_label!` directly into a pre-sized `IOBuffer` without intermediate `"$" * ... * "$"` string allocations.
- **Verification**: Added regression tests in [`test/test_ipedraw.jl`](../test/test_ipedraw.jl) covering all 9 anchors, invalid anchor `ArgumentError`, XML entity escaping in `:math` and `:normal` modes, and `_escape_ipe_xml`.

---

## Verification Evidence

- `./tools/build`: Exited 0 (clean precompilation of `BasicCompGeometry` and its extensions).
- `./tools/gate`: Exited 0; all Aqua checks passed (`Method ambiguity`, `Unbound type parameters`, `Undefined exports`, `Compare Project.toml and test/Project.toml`, `Stale dependencies`, `Compat bounds`, `Piracy`, `Persistent tasks`) and `BasicCompGeometry.jl` test suite passed (`156572 / 156572` assertions in 1m30.7s).
- [`CHANGELOG.md`](../CHANGELOG.md): Updated under `## [Unreleased]`.
