# Post-Work Summary — Review #011 (Performance)

- **Review Report**: [reviews/011_performance.md](reviews/011_performance.md)
- **Plan**: [reviews/011_performance_plan.md](reviews/011_performance_plan.md)
- **Date**: 2026-11-14
- **Status**: Complete

---

## Issue-by-Issue Disposition

### 1. `read_plt` / `_parse_numbers` (`src/ReadWrite.jl`)
- **Status**: Genuine defect — **Fixed**.
- **Mitigation**:
  - Replaced `_parse_numbers` with UTF-8-safe `_parse_next_float` and `_parse_point`, using a stack-allocated `MVector{D,Float64}` inside a `1:D` loop (avoiding the `Core.Box` capture present in an `ntuple` closure) and parsing only the first `D` coordinates via `SubString` views.
  - Added `dim >= 1` validation and `Base.@constprop :aggressive` on `read_plt`, dispatching to a `Val{D}` function barrier `_read_plt` that streams lines with `eachline(filename)` instead of `readlines(filename)`.
- **Verification**:
  - `@code_warntype BasicCompGeometry.ReadWrite.read_plt("x.plt")` and `_read_plt("x.plt", Val(2))` confirm concrete `PntSeq{2,Float64}` inference with zero `Core.Box` captures.
  - Added regression tests in [test/test_curve_ops.jl](test/test_curve_ops.jl) checking `@inferred read_plt(...)` for 2D and 3D, `dim < 1` `ArgumentError`, UTF-8 whitespace / trailing commas / extra columns, and `(@allocated _parse_point(...)) == 0`.

### 2. `read_geolife_plt` (`src/ReadWrite.jl`)
- **Status**: Genuine defect — **Fixed**.
- **Mitigation**:
  - Switched from `readlines(filename)` to `eachline(filename)` inside `Iterators.drop(..., 6)`.
  - Replaced `split(line, ",")` with `findfirst`/`findnext` comma indexing and `SubString` slicing, throwing on missing commas while preserving support for 2-column records.
- **Verification**:
  - `@code_warntype BasicCompGeometry.ReadWrite.read_geolife_plt("x.plt")` confirms concrete `PntSeq{2,Float64}` inference with `SubString{String}` views and no array allocations per line.
  - Added 2-column and malformed-record regression tests in [test/test_curve_ops.jl](test/test_curve_ops.jl).

### 3. `write_fbin` and `write_ibin` (`src/ReadWrite.jl`)
- **Status**: Genuine defect — **Fixed**.
- **Mitigation**:
  - Extracted shared `_write_bin(filename::String, vecs::Matrix{T}) where {T<:Union{Float32,Int32}}` using `open(filename, "w") do f ... end`.
  - Replaced `write(f, vecs')` (which dispatched to scalar `write(::IO, ::AbstractArray)`) with chunked transposition into a bounded (≤ 1 MiB) contiguous buffer written via `write(f, @view(buf[:, 1:count]))` (`write(::IO, ::StridedArray)` bulk `unsafe_write`).
- **Verification**:
  - Benchmarked on a $100{,}000 \times 128$ `Float32` matrix: reduced write time from `0.219 s` to `0.057 s` while bounding buffer allocation to `1.004 MiB` (vs `48.833 MiB` for full `permutedims`).
  - Added multi-chunk (`5000 × 64` > `262144` elements) and empty-matrix (`0 × 3`) round-trip tests in [test/test_read_write.jl](test/test_read_write.jl).

### 4. `BBT.Tree_draw` / `node_draw` (`ext/BBTCairoExt.jl`)
- **Status**: Genuine defect — **Fixed**.
- **Mitigation**:
  - Precomputed concrete color constants `_BBT_FILL_COLOR` (`Colors.RGBA{Float64}`), `_BBT_STROKE_COLOR` (`Colors.RGB{Float64}`), and `_BBT_BG_COLOR` (`Colors.RGB{Float64}`), removing repeated `Colors.parse(Colors.Colorant, ...)` lookups inside rendering loops.
  - Extracted top-level helpers `_collect_bbt_levels!`, `_draw_bbt_levels!`, `_draw_bbt_bg`, and `_draw_bbt_node` to pre-collect internal nodes by depth in a single $O(N)$ pass and render all pages in $O(N)$ total work without $\Theta(d \cdot N)$ unpruned traversals. Added `ArgumentError` guard when `tree.root === nothing`.
- **Verification**:
  - Added `BBT.Tree_draw (BBTCairoExt)` regression tests in [test/test_canvas.jl](test/test_canvas.jl) exercising multi-level expanded trees and single-point (depth 1) trees.

---

## Quality Gate Results

- `./tools/build`: passed (`using BasicCompGeometry`).
- `./tools/gate`: passed all **156,588 / 156,588** tests (including `Aqua.test_all`, `CairoExt`, `BBTCairoExt`, `LaTeXCairoExt`, and all newly added regression tests).
- [CHANGELOG.md](CHANGELOG.md): updated under `## [Unreleased]`.
