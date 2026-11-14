# Performance Remediation Plan #011

- **Source Report**: [reviews/011_performance.md](reviews/011_performance.md)
- **Target Scope**: `src/ReadWrite.jl`, `ext/BBTCairoExt.jl`, regression tests in `test/`
- **Date**: 2026-11-14

---

## 1. Triage Summary

All four findings in [reviews/011_performance.md](reviews/011_performance.md) were reproduced and verified against the current codebase. None are false positives or intentional design trade-offs.

| # | Location | Severity | Triage | Verification Evidence |
|---|---|---|---|---|
| 1 | `src/ReadWrite.jl` (`_parse_numbers`, `read_plt`) | Moderate | **Genuine** | `@code_warntype BasicCompGeometry.ReadWrite.read_plt("x.plt")` infers `PntSeq{_A, Float64} where _A` (non-concrete `UnionAll`) due to runtime `dim` type parameterization; `readlines` materializes all lines; regex `split` + comprehension allocates two temporary vectors per line and parses unused columns. |
| 2 | `src/ReadWrite.jl` (`read_geolife_plt`) | Moderate | **Genuine** | `Iterators.drop(readlines(filename), 6)` eagerly loads the entire file into a `Vector{String}`, and `split(line, ",")` allocates a 7-element `Vector{SubString{String}}` per data line when only fields 1–2 are used. |
| 3 | `src/ReadWrite.jl` (`write_fbin`, `write_ibin`) | Moderate | **Genuine** | `@which write(devnull, adjoint(zeros(Float32, 2, 2)))` resolves to `write(s::IO, A::AbstractArray) @ Base io.jl:844` (scalar element-by-element writes with stride `nvecs`), whereas `StridedArray` resolves to `write(s::IO, A::StridedArray) @ Base io.jl:857` (bulk `unsafe_write`). |
| 4 | `ext/BBTCairoExt.jl` (`BBT.Tree_draw`, `node_draw`) | Moderate | **Genuine** | `node_draw` does not prune recursion when `level >= last(range)`, traversing all $N$ nodes on each of the $d + 1$ pages ($\Theta(d \cdot N)$ / $\Theta(N^2)$ unbalanced); `Clr.parse(Clr.Colorant, ...)` performs abstract-typed string dictionary lookups per node and per page. |

---

## 2. Architectural & Implementation Plan

### Issue 1: Type-stable, streaming `read_plt` (`src/ReadWrite.jl`)
- **Note on Report Suggestion**: The report's proposed `_parse_first_n` captures a mutable `idx` inside an `ntuple(D) do _ ... end` closure, which `@code_warntype` shows introduces a `Core.Box` heap wrapper. We avoid that capture.
- **Design**:
  1. Validate `dim >= 1 || throw(ArgumentError("dim must be >= 1, got $dim"))` in `read_plt`.
  2. Mark `read_plt(filename::AbstractString; dim::Int = 2)` with `Base.@constprop :aggressive` and dispatch to a `Val{D}` function barrier `_read_plt(filename::AbstractString, ::Val{D}) where {D}`.
  3. Stream lines lazily with `eachline(filename)` instead of `readlines(filename)`.
  4. Replace `_parse_numbers` with two small helpers:
     - `_parse_next_float(text, idx, last_idx, dim, raw_line)`: skips `,` and whitespace using UTF-8-safe `firstindex`/`nextind`/`prevind`, throws on premature end-of-line, and parses `SubString(text, start_idx, prevind(text, idx))`.
     - `_parse_point(text, ::Val{D}, raw_line)`: populates a stack-allocated `MVector{D,Float64}(undef)` in a `1:D` loop and returns `Point{D,Float64}(buf)` (zero heap allocations, no `Core.Box`, handles arbitrary `D` without compile-time recursion blowup).
- **Public API**: `read_plt(filename::AbstractString; dim::Int = 2)` signature and behavior unchanged.

### Issue 2: Streaming, allocation-free `read_geolife_plt` (`src/ReadWrite.jl`)
- **Design**:
  1. Replace `readlines(filename)` with `eachline(filename)` inside `Iterators.drop(eachline(filename), 6)`.
  2. Locate the first two commas with `findfirst(',', line)` and `findnext(',', line, y_start)`. Throw on missing first comma; if the second comma is absent (2-column record), slice through `lastindex(line)` for backward compatibility with `split(line, ",")`.
  3. Parse both `SubString` views directly into `point(x, y)` without allocating intermediate arrays.

### Issue 3: Bulk-buffered `write_fbin` and `write_ibin` (`src/ReadWrite.jl`)
- **Design**:
  1. Extract a shared parameterized helper `_write_bin(filename::String, vecs::Matrix{T}) where {T<:Union{Float32,Int32}}` wrapped in `open(filename, "w") do f ... end` for exception-safe resource cleanup.
  2. Rather than `permutedims(vecs)` (which duplicates the entire matrix in memory, e.g. 512 MiB for $10^6 \times 128$ `Float32`), allocate a bounded reusable transpose buffer `Matrix{T}(undef, dim, chunk_rows)` with `chunk_rows = clamp(div(262144, dim), 1, nvecs)` (≤ 1 MiB).
  3. Transpose row-blocks into `buf` and call `write(f, @view(buf[:, 1:count]))`, which dispatches to `write(s::IO, A::StridedArray)` (bulk `unsafe_write`).
  4. Keep `write_fbin` and `write_ibin` public signatures unchanged as thin callers of `_write_bin`.

### Issue 4: $O(N)$ level-indexed `BBT.Tree_draw` with concrete colors (`ext/BBTCairoExt.jl`)
- **Design**:
  1. Precompute concrete color values outside loops:
     - `fill_color = Clr.RGBA{Float64}(1.0, 1.0, 0.0, 0.1)` (`"yellow"` at 10% alpha)
     - `stroke_color = Clr.RGB{Float64}(0.0, 0.0, 1.0)` (`"blue"`)
     - `bg_color = Clr.RGB{Float64}(173 / 255, 216 / 255, 230 / 255)` (`"lightblue"`)
  2. Extract top-level helpers `_draw_bbt_bg`, `_draw_bbt_node`, and `_collect_bbt_levels!` so `BBT.Tree_draw` stays flat (`nesting <= 2`) and free of closure-capture overhead.
  3. Collect non-leaf nodes per depth level in one $O(N)$ DFS pass into `levels::Vector{Vector{NodeT}}`, then render Page 1 (`0:d`) and each level-pair page (`i:(i+1)` for `i = (d-1):-1:0`) directly from `levels` in $O(N)$ total work.
  4. Guard `isnothing(tree.root)` with `ArgumentError`.

---

## 3. Verification & Regression Test Plan

1. **Unit & Regression Tests**:
   - In `test/test_curve_ops.jl`:
     - Verify `@inferred read_plt(f)` returns `PntSeq{2,Float64}` and `@inferred read_plt(h; dim=3)` returns `PntSeq{3,Float64}` (type stability).
     - Verify UTF-8 / trailing-comma / invalid `dim < 1` / malformed GeoLife record handling and multi-line streaming correctness.
     - Verify bounded allocations for `read_geolife_plt` and `read_plt` on repeated calls (warmed `@allocated`).
   - In `test/test_read_write.jl`:
     - Verify `write_fbin` / `write_ibin` round-trip accuracy for non-trivial matrix sizes, empty matrices (`0 x 3`, `3 x 0`), and chunk-boundary row counts (e.g. `> chunk_rows` so multi-chunk writing is exercised).
   - In `test/test_canvas.jl` (runs under `tools/gate` where `Cairo` and `Colors` are test dependencies):
     - Add a regression test exercising `BBT.Tree_draw` on a 2D `BBT.Tree` (both shallow/empty-sublevel and multi-level expanded trees) producing a valid multi-page PDF.
2. **Quality Gates**:
   - Run `./tools/build` and `./tools/gate` (`Pkg.test()` including `Aqua.test_all`).
   - Update `CHANGELOG.md` under `## [Unreleased]`.
   - Author post-work summary in `reviews/011_summary.md`.
