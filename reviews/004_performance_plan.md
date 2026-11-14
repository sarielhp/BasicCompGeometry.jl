# Architectural Plan: Performance Review #004 Remediation

- **Review Document**: `reviews/004_performance.md`
- **Focus Area**: Performance Audit (Type Stability & Allocations in `VisibilityPolygon.jl` and `MinkowskiSum.jl`)
- **Target Files**:
  - `src/VisibilityPolygon.jl`
  - `src/MinkowskiSum.jl`
  - `test/test_visibility_polygon.jl`
  - `test/test_minkowski_sum.jl`
  - `CHANGELOG.md`
- **Output Artifacts**:
  - `reviews/004_performance_plan.md`
  - `reviews/004_summary.md`

---

## 1. Triage of Findings

All 6 reported findings in `reviews/004_performance.md` have been triaged and reproduced against the current codebase via `@code_warntype`, `@allocated`, and integer/float32 test invocations. All 6 are confirmed as **genuine defects**:

### Finding 1: `Core.Box` Captures & Untyped Dispatch in `visibility_polygon`
- **Location**: `src/VisibilityPolygon.jl:57-173`
- **Severity**: Major
- **Triage Assessment**: **GENUINE DEFECT**
- **Verified Reproduction**:
  - `visibility_polygon(pts40, point(0.0, 0.0))` on a 40-vertex `Float64` star polygon allocates **656,456 bytes**.
  - Conditional reassignment `R = Rational{BigInt}` captured by the comprehension `[point(R(p.x), R(p.y)) for p in poly]` boxes `R` in `Core.Box`, and `u_start`/`u_end` assigned in `elseif` branches and captured by `findfirst` closures create two more `Core.Box` slots.
- **Remediation**:
  - Resolve coordinate type `R` cleanly via `_vis_coord_type(kernel, _scalar_type(poly), eltype(q))` and pass `R` as a `::Type{R}` static parameter into a function barrier `_visibility_polygon_impl`.
  - Extract small, shallow-nested helpers: `_find_ray_index` (closure-free loop), `_sector_ray_range`, `_closest_sector_hit`, `_is_closer_sector_hit`, `_dedup_consecutive!`, and `_prune_collinear`.
  - Iterate directly over consecutive vertex pairs of `poly_conv` without allocating an intermediate `edges` tuple vector.

### Finding 2: `Core.Box` Capture & Per-Iteration Allocations in `polygon_kernel`
- **Location**: `src/VisibilityPolygon.jl:301-395`
- **Severity**: Major
- **Triage Assessment**: **GENUINE DEFECT**
- **Verified Reproduction**:
  - `polygon_kernel(pts40)` on a 40-vertex `Float64` star polygon allocates **385,384 bytes**.
  - Conditional reassignment of `R` before the comprehension boxes `R`, causing dynamic dispatch throughout the $O(n \cdot m)$ Sutherland–Hodgman clipping loop. Additionally, each outer iteration allocates two new vectors (`next_poly`, `deduped`), evaluates `turn_sign` twice per vertex, and duplicates the segment-halfplane intersection calculation.
- **Remediation**:
  - Resolve `R` via `_vis_coord_type` and dispatch into `_polygon_kernel_impl(poly, kernel, R)`.
  - Extract `_bounding_rect`, `_clip_polygon_halfplane!`, `_halfplane_edge_intersect`, and `_halfplane_inside`.
  - Preallocate and reuse `next_poly` and `deduped` scratch buffers across the `1:n` loop via buffer swapping, and reuse `_dedup_consecutive!` and `_prune_collinear`.

### Finding 3: Return-Type Instability in `_ray_edge_intersect` & `InexactError` on Integer Inputs
- **Location**: `src/VisibilityPolygon.jl:124-127, 454-492`
- **Severity**: Major
- **Triage Assessment**: **GENUINE DEFECT**
- **Verified Reproduction**:
  - `@code_warntype` on `_ray_edge_intersect(..., Float32)` returns `Union{Tuple{Bool, Float32, SVector{2, Float32}}, Tuple{Bool, Float32, SVector{2, Float64}}}` because `clamp(s, 0.0, 1.0)` promotes `Float32` to `Float64` on the hit branch.
  - Calling `visibility_polygon` or `polygon_kernel` on `Vector{Point{2, Int64}}` with `FloatingKernel()` throws `InexactError(:Int64, ...)`.
  - `best_t1` / `best_t2` were initialized as `typemax(Rational{BigInt})` or `Inf` (`Float64`) regardless of `R`, introducing `Union{Rational{BigInt}, Rational{Int64}}` for `Rational{Int64}` under `ExactKernel()`.
- **Remediation**:
  - Have `_vis_coord_type(::AbstractKernel, T1, T2)` promote via `float(promote_type(T1, T2))` for non-exact kernels while preserving `Rational` subtypes (or `Rational{BigInt}`) for `ExactKernel`.
  - Initialize `best_t1 = best_t2 = kernel isa ExactKernel ? typemax(R) : R(Inf)` so `best_t1`/`best_t2` always have type `R`.
  - In `_ray_edge_intersect`, use `clamp(s, zero(t), one(t))` and convert `t` and `pt` coordinates to `R`, guaranteeing homogeneous `Tuple{Bool, R, Point{2, R}}` return types across all branches.

### Finding 4: `Core.Box` Captures & Excessive Allocations in `_minkowski_convex` and `minkowski_sum_decomposition`
- **Location**: `src/MinkowskiSum.jl:181-263, 346-383`
- **Severity**: Major
- **Triage Assessment**: **GENUINE DEFECT**
- **Verified Reproduction**:
  - `_minkowski_convex(sq, sq, FloatingKernel())` allocates **3,744 bytes**; `minkowski_sum(sq, sq)` allocates **10,784 bytes**; `minkowski_sum_decomposition(pts40, pts40)` allocates **5,114,472 bytes**.
  - Reassignments of `p1`, `p2` in `_minkowski_convex` and `R` in `minkowski_sum_decomposition` create `Core.Box` slots and untyped `Any`/`Vector` dispatch across the $O(n \cdot m)$ loop.
  - `_minkowski_convex` allocates up to 11 temporary vectors per call (`p1`, `p2`, dedups, `v1`, `v2`, `edges1`, `edges2`, `result`, `deduped`, `final_poly`), and `minkowski_sum_decomposition` allocates per-triangle `Vector`s.
- **Remediation**:
  - Introduce `_minkowski_coord_type` (preserving `Int` under `FloatingKernel()` since Minkowski sums of convex polygons only add coordinates, and widening to `Rational{BigInt}` for non-rational inputs under `ExactKernel()`).
  - Use function barriers `_minkowski_convex_impl` and `_minkowski_sum_decomposition_impl` parameterized on `::Type{R}`.
  - Replace `v1`, `v2`, `edges1`, `edges2` with direct cyclic indexing over `p1`/`p2` in `_minkowski_merge_ccw!`.
  - In `minkowski_sum_decomposition`, wrap triangulated indices with a zero-allocation `_TrianglePieces{R} <: AbstractVector{SVector{3, Point{2, R}}}` view, dispatch cleanly into `_minkowski_pair_loop!` without `Union` types, and reuse `scratch_raw` / `scratch_dedup` buffers across all `(c1, c2)` pairs.

### Finding 5: Three `Core.Box` Captures & Redundant Allocations in `is_convex`
- **Location**: `src/MinkowskiSum.jl:57-115`
- **Severity**: Moderate
- **Triage Assessment**: **GENUINE DEFECT**
- **Verified Reproduction**:
  - `is_convex(sq)` allocates **3,520 bytes** due to `Core.Box` slots for `R`, `deduped`, and `idx_min`, plus two vector allocations (`edges`, `nz_edges`) for edge-angle monotonicity checking.
- **Remediation**:
  - Dispatch `is_convex` into `_is_convex_impl(poly, kernel, R)`, deduplicate directly in one pass via `_collect_deduped_points!`, reverse in-place (`reverse!(deduped)`), and extract `_has_all_ccw_turns`, `_lowest_left_index`, and `_edges_angle_monotone` to verify monotonicity with zero intermediate vector allocations.

### Finding 6: Abstract `Point{2}[...]` Containers in `minkowski_sum` Overloads
- **Location**: `src/MinkowskiSum.jl:304-317`
- **Severity**: Moderate
- **Triage Assessment**: **GENUINE DEFECT**
- **Verified Reproduction**:
  - `minkowski_sum(s1, s2)` on `Segment{2, Float64}` inputs allocates **8,832 bytes** because `Point{2}[...]` creates `Vector{SVector{2}}` with abstract `eltype`.
- **Remediation**:
  - Add `@inline _pts_vec(s::Segment{2, T}) where {T} = Point{2, T}[s.p, s.q]` and `@inline _pts_vec(p::Point{2, T}) where {T} = Point{2, T}[p]`, and use them in all five `Point{2}`/`Segment{2}` overloads.

---

## 2. Quality, Architecture & Compatibility Constraints

- **Public API**: Preserve all exported names, signatures, keyword arguments, and return semantics (except fixing `InexactError` on integer inputs to `visibility_polygon`/`polygon_kernel` by returning `Float64` points under `FloatingKernel()`).
- **Code Structure**: Keep every new and modified function small, single-purpose, with nesting depth $\le 4$, and share `_scalar_type`, `_collect_deduped_points!`, `_dedup_consecutive!`, and `_prune_collinear` across `VisibilityPolygon.jl` and `MinkowskiSum.jl` (included consecutively in `src/BasicCompGeometry.jl`).
- **Anti-Oscillation**: Leave `src/MonotoneTriangulation.jl`, `src/SweepTriangulation.jl`, and all other unrelated files untouched.

---

## 3. Verification & Regression Test Plan

1. **Regression Unit Tests**:
   - In `test/test_visibility_polygon.jl`:
     - Verify `Float32` type stability (`@inferred`, element type `Point{2, Float32}`) for `visibility_polygon`, `polygon_kernel`, and `_ray_edge_intersect`.
     - Verify `Int64` inputs under `FloatingKernel()` succeed without `InexactError` for `visibility_polygon` and `polygon_kernel` and return `Point{2, Float64}`.
     - Verify `Rational{Int64}` inputs under `ExactKernel()` preserve `Point{2, Rational{Int64}}` and type stability (`@inferred`).
     - Verify `@inferred` return types on `visibility_polygon` and `polygon_kernel` for `Float64` and `ExactKernel()`, plus allocation sanity bounds on a 40-vertex star polygon.
   - In `test/test_minkowski_sum.jl`:
     - Verify `@inferred` return types for `is_convex`, `minkowski_sum` (polygons and `Segment{2}` / `Point{2}` overloads), and `minkowski_sum_decomposition`.
     - Verify `Rational{Int64}` and `Int64` preservation under `ExactKernel()` / `FloatingKernel()` respectively.
     - Verify allocation bounds on `minkowski_sum(s1, s2)`, `minkowski_sum(sq, sq)`, `is_convex(sq)`, and `minkowski_sum_decomposition(star40, star40)`.
2. **Quality Gates**:
   - Run `tools/build` and `tools/gate` (`Pkg.test()`) after edits and confirm all tests pass.
   - Re-measure `@code_warntype` and `@allocated` on all cases from `reviews/004_performance.md`.
