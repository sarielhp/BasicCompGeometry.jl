# Code Review Remediation Summary #004 (Lens: performance)

- **Review Document**: `reviews/004_performance.md`
- **Architectural Plan**: `reviews/004_performance_plan.md`
- **Date**: 2026-11-13
- **Focus Area**: Performance Audit (`VisibilityPolygon.jl` and `MinkowskiSum.jl`)
- **Status**: Complete & Verified

---

## 1. Executive Summary

All 6 issues reported in `reviews/004_performance.md` were triaged, reproduced, and confirmed as genuine performance or type-stability defects. Minimal, idiomatic fixes were implemented in `src/VisibilityPolygon.jl` and `src/MinkowskiSum.jl`, with regression unit tests added to `test/test_visibility_polygon.jl` and `test/test_minkowski_sum.jl`.

### Measured Improvements (Before vs. After)
- `visibility_polygon` on 40-vertex `Float64` star polygon: **656,456 B -> 4,896 B** (~134x reduction); on 200-vertex star polygon: **22,824 B** (previously ~15.0 MB reported).
- `polygon_kernel` on 40-vertex `Float64` star polygon: **385,384 B -> 5,024 B** (~77x reduction); on 200-vertex star polygon: **20,576 B** (previously ~7.0 MB reported).
- `_minkowski_convex(sq, sq, FloatingKernel())`: **3,744 B -> 864 B**.
- `minkowski_sum(sq, sq)`: **10,784 B -> 1,120 B**.
- `minkowski_sum(s1, s2)` on `Segment{2, Float64}`: **8,832 B -> 960 B**.
- `is_convex(sq)`: **3,520 B -> 128 B**.
- `minkowski_sum_decomposition` on two 40-vertex star polygons: **5,114,472 B -> 274,152 B** (~18.7x reduction).
- `@code_warntype` confirms zero `Core.Box` slots and zero `::Any` variables across `visibility_polygon`, `_visibility_polygon_impl`, `polygon_kernel`, `_polygon_kernel_impl`, `_ray_edge_intersect`, `is_convex`, `_is_convex_impl`, `_minkowski_convex`, `minkowski_sum_decomposition`, and `_minkowski_sum_decomposition_impl`.

---

## 2. Detailed Findings, Mitigations & Verification

### Finding 1: `Core.Box` Captures & Untyped Dispatch in `visibility_polygon`
- **Location**: `src/VisibilityPolygon.jl` (`visibility_polygon`)
- **Status**: **Mitigated**
- **Remediation**:
  - Resolved coordinate type `R` via `_vis_coord_type(kernel, _scalar_type(poly), eltype(q))` and dispatched to function barrier `_visibility_polygon_impl(poly, q, kernel, R)`.
  - Eliminated `u_start` / `u_end` closure captures using `_find_ray_index` and `_sector_ray_range`.
  - Extracted `_sorted_unique_rays`, `_closest_sector_hit`, `_is_closer_sector_hit`, `_collect_deduped_points!`, and `_prune_collinear`; removed intermediate `edges` tuple vector allocation.
- **Verification**:
  - `@inferred` passes for `Float64`, `Float32`, `Int64`, and `Rational{Int64}` inputs; 40-vertex star allocation guard (`< 32_768` bytes) passes in `test/test_visibility_polygon.jl`.

### Finding 2: `Core.Box` Capture & Per-Iteration Allocations in `polygon_kernel`
- **Location**: `src/VisibilityPolygon.jl` (`polygon_kernel`)
- **Status**: **Mitigated**
- **Remediation**:
  - Resolved `R` via `_vis_coord_type` and dispatched to function barrier `_polygon_kernel_impl(poly, kernel, R)`.
  - Extracted `_bounding_rect`, `_clip_polygon_halfplane!`, `_halfplane_edge_intersect`, `_halfplane_thresh`, and `_halfplane_inside` (evaluating `turn_sign` once per vertex).
  - Reused `next_poly` and `deduped` scratch buffers across the `1:n` edge loop via buffer swapping.
- **Verification**:
  - `@inferred` passes for `Float64`, `Float32`, `Int64`, and `Rational{Int64}`; 40-vertex star allocation guard (`< 32_768` bytes) passes in `test/test_visibility_polygon.jl`.

### Finding 3: Return-Type Instability in `_ray_edge_intersect` & `InexactError` on Integer Inputs
- **Location**: `src/VisibilityPolygon.jl` (`visibility_polygon`, `polygon_kernel`, `_ray_edge_intersect`)
- **Status**: **Mitigated**
- **Remediation**:
  - `_vis_coord_type` promotes non-`ExactKernel` coordinate types via `float(promote_type(T1, T2))` so integer inputs compute intersections in `Float64` without throwing `InexactError`.
  - Initialized `best_t1` / `best_t2` in `_closest_sector_hit` using `kernel isa ExactKernel ? typemax(R) : R(Inf)`.
  - Updated `_ray_edge_intersect` to use `clamp(s, zero(t), one(t))` and convert hit results to `(true, R(t), pt)` with `Point{2, R}` coordinates.
- **Verification**:
  - Regression tests in `test/test_visibility_polygon.jl` verify `@inferred` `Tuple{Bool, Float32, Point{2, Float32}}` on both hit and miss paths of `_ray_edge_intersect`, `Float64` promotion on `Int64` polygons, and `Rational{Int64}` preservation under `ExactKernel()`.

### Finding 4: `Core.Box` Captures & Excessive Allocations in `_minkowski_convex` and `minkowski_sum_decomposition`
- **Location**: `src/MinkowskiSum.jl` (`_minkowski_convex`, `minkowski_sum_decomposition`)
- **Status**: **Mitigated**
- **Remediation**:
  - Added `_minkowski_coord_type` and function barriers `_minkowski_convex_impl` and `_minkowski_sum_decomposition_impl`.
  - Eliminated `v1`, `v2`, `edges1`, and `edges2` allocations by cyclically indexing `p1`/`p2` directly in `_minkowski_merge_ccw!`.
  - Represented triangulated pieces via zero-allocation `_TrianglePieces{R} <: AbstractVector{SVector{3, Point{2, R}}}`, dispatched concretely through `_dispatch_p1!` / `_dispatch_p2!` / `_minkowski_pair_loop!`, and reused `scratch_raw` / `scratch_dedup` buffers across triangle pairs.
- **Verification**:
  - Regression tests in `test/test_minkowski_sum.jl` verify `@inferred` types for `Float64`, `Int64`, and `Rational{Int64}`, plus allocation bounds (`< 4096` bytes for `minkowski_sum(sq, sq)` and `< 1_000_000` bytes for 40-vertex star `minkowski_sum_decomposition`).

### Finding 5: Three `Core.Box` Captures & Redundant Allocations in `is_convex`
- **Location**: `src/MinkowskiSum.jl` (`is_convex`)
- **Status**: **Mitigated**
- **Remediation**:
  - Dispatched `is_convex` into `_is_convex_impl(poly, kernel, R)`, deduplicating in one pass via `_collect_deduped_points!` and reversing in-place (`reverse!(deduped)`).
  - Extracted `_has_all_ccw_turns`, `_lowest_left_index`, and `_edges_angle_monotone` to verify angular monotonicity without allocating `edges` or `nz_edges`.
- **Verification**:
  - Regression tests in `test/test_minkowski_sum.jl` verify `@inferred is_convex(sq)` (`Float64` and `Rational{Int64}`) and allocation bound `< 1024` bytes.

### Finding 6: Abstract `Point{2}[...]` Containers in `minkowski_sum` Overloads
- **Location**: `src/MinkowskiSum.jl` (`minkowski_sum` overloads for `Point{2}` and `Segment{2}`)
- **Status**: **Mitigated**
- **Remediation**:
  - Added `@inline _pts_vec(s::Segment{2, T}) where {T} = Point{2, T}[s.p, s.q]` and `@inline _pts_vec(p::Point{2, T}) where {T} = Point{2, T}[p]` and replaced all `Point{2}[...]` calls in the five `Point{2}`/`Segment{2}` overloads.
- **Verification**:
  - Regression tests in `test/test_minkowski_sum.jl` verify `@inferred` return types on all `Segment{2}` / `Point{2}` overloads and allocation bound `< 4096` bytes for `minkowski_sum(seg1, seg2)`.

---

## 3. Quality Gate & Metrics Verification

- **Build (`tools/build`)**: Passed cleanly.
- **Full Test Suite (`tools/gate`, `Pkg.test()`)**: Passed (`156,383 / 156,383` assertions in `BasicCompGeometry.jl`, plus all `Aqua` package checks).
- **Sizing & Nesting**: All new and modified functions in `src/VisibilityPolygon.jl` and `src/MinkowskiSum.jl` are `<= 43` lines with maximum nesting depth `<= 3` (`<= 4` threshold).
- **Scope & Anti-Oscillation**: Unrelated modules (`src/MonotoneTriangulation.jl`, `src/SweepTriangulation.jl`, etc.) were left untouched. Public APIs and exports are unchanged. User-visible changes are recorded in `CHANGELOG.md` under `## [Unreleased]`.
