# Architectural Plan: Performance Review #007 Remediation

- **Review Document**: `reviews/007_performance.md`
- **Focus Area**: Performance Audit (Type Stability, Compile-Time Dimension Dispatch, and Heap Allocation Elimination)
- **Target Files**:
  - `src/Points.jl`
  - `src/Segments.jl`
  - `src/Transforms2D.jl`
  - `src/BBoxes.jl`
  - `test/runtests.jl`
  - `test/test_mat_pnt_seq.jl`
  - `CHANGELOG.md`
- **Output Artifacts**:
  - `reviews/007_performance_plan.md`
  - `reviews/007_summary.md`

---

## 1. Triage of Findings

All 7 findings in `reviews/007_performance.md` have been evaluated against the current codebase and verified via Julia type-inference and allocation checks. All 7 are **genuine defects**:

### Finding 1: Union Return Type in `dist_subspace(p, q, dir)`
- **Location**: `src/Points.jl:56-62`
- **Severity**: Moderate
- **Triage Assessment**: **GENUINE DEFECT**
- **Verified Reproduction**:
  - `Base.return_types(BasicCompGeometry.dist_subspace, (Point{2,Int}, Point{2,Int}, Point{2,Int}))` returns `Any[Union{Float64, Int64}]` because the positive branch returns `sqrt(d_sq)` (`Float64`) while the nonpositive branch returns `zero(T)` (`Int64`).
- **Remediation**:
  - Replace `return d_sq > 0 ? sqrt(d_sq) : zero(T)` with `return sqrt(max(zero(d_sq), d_sq))` so both paths share the exact return type of `sqrt(d_sq)`.

### Finding 2: Uninferred `D` and Heap Allocation in `rand_gaussian` / `rand_point`
- **Location**: `src/Points.jl:199-208`
- **Severity**: Moderate
- **Triage Assessment**: **GENUINE DEFECT**
- **Verified Reproduction**:
  - `Base.return_types(d -> rand_point(d), (Int,))` infers `SVector{_A, Float64} where _A` and `Base.return_types(d -> rand_gaussian(d), (Int,))` infers `Any` when `D` is non-constant.
  - Even for constant `D = 3`, `rand_gaussian(3)` allocates 80 bytes per call due to `randn(T, D)` creating a temporary heap `Vector{T}`.
- **Remediation**:
  - Add `Val{D}` overloads: `rand_gaussian(::Val{D}, ::Type{T} = Float64) where {D,T<:Real}` using `Point{D,T}(ntuple(_ -> randn(T), Val(D)))` (zero heap allocations) and `rand_point(::Val{D}, ::Type{T} = Float64) where {D,T<:Real} = rand(Point{D,T})`.
  - Retain the `D::Int` convenience overloads for public API compatibility, delegating `rand_gaussian(D::Int, ::Type{T} = Float64)` to `rand_gaussian(Val(D), T)` and specializing on `::Type{T}` in both `rand_gaussian` and `rand_point`.

### Finding 3: Value-Dependent Union Return Type in `nn_point(s_p, s_q, qr)`
- **Location**: `src/Segments.jl:122-130`
- **Severity**: Moderate
- **Triage Assessment**: **GENUINE DEFECT**
- **Verified Reproduction**:
  - `Base.return_types(nn_point, (Point{2,Int}, Point{2,Int}, Point{2,Int}))` returns `Union{SVector{2,Float64}, SVector{2,Int64}}`.
  - `Base.return_types(nn_point, (Point{2,Float32}, Point{2,Float32}, Point{2,Float32}))` returns `Union{SVector{2,Float32}, SVector{2,Float64}}` because `clamp(t, 0.0, 1.0)` widens `Float32` to `Float64` on the non-degenerate branch.
- **Remediation**:
  - Compute `R = typeof(one(T) / one(T))`, convert `s_p`, `s_q`, `qr` to `Point{D,R}`, check `iszero(d_sq) && return p`, and clamp with `clamp(t, zero(t), one(t))`. This yields a uniform concrete `Point{D,R}` return type for `Int`, `Float32`, `Float64`, and `Rational`.

### Finding 4: Union Return Type in `bisection_point(seg, p, q)`
- **Location**: `src/Segments.jl:210-225`
- **Severity**: Moderate
- **Triage Assessment**: **GENUINE DEFECT**
- **Verified Reproduction**:
  - `Base.return_types(BasicCompGeometry.bisection_point, (Segment{2,Int}, Point{2,Int}, Point{2,Int}))` returns `Union{Tuple{Bool,Float64,SVector{2,Float64}}, Tuple{Bool,Float64,SVector{2,Int64}}}`.
  - For `Float32`, it returns `Union{Tuple{Bool,Float32,SVector{2,Float32}}, Tuple{Bool,Float64,SVector{2,Float32}}}`.
- **Remediation**:
  - Promote coordinates to `R = typeof(one(T) / one(T))`, use `R(2)`, `zero(R)`, `one(R)`, and return `(false, zero(R), sp)` on the parallel branch.
  - Guard the parallel tolerance check with `R <: AbstractFloat ? abs(denom) < R(1e-12) : iszero(denom)` so exact `Rational` types do not convert `1e-12` into a large-denominator rational.

### Finding 5: Heap Matrix/Vector Allocations in `dist_segment_segment`
- **Location**: `src/Segments.jl:146-188`
- **Severity**: Moderate
- **Triage Assessment**: **GENUINE DEFECT**
- **Verified Reproduction**:
  - Warmed `dist_segment_segment` allocates 592 bytes per call from constructing `c`, `m`, and `m \ c`, and widens `Float32` via `0.0`/`1.0` clamp bounds.
- **Remediation**:
  - Convert endpoints to `Point{D,R}` with `R = typeof(one(T) / one(T))`, solve the 2x2 normal equations via scalar Cramer's rule with zero heap allocations, and extract `_min_endpoint_segment_dist` to avoid duplicating the 4-way endpoint-segment fallback.
  - Use a scale-invariant degeneracy check `Δ <= tol * (A * C)` (with `tol = R <: AbstractFloat ? max(R(1e-12), 4 * eps(R)) : zero(R)`) rather than `max(..., one(float(Δ)))` so small nonparallel segments (`length <= 1e-3`) are not misclassified as parallel.

### Finding 6: Heap Allocations in `apply_transform(m, p)`
- **Location**: `src/Transforms2D.jl:44-48`
- **Severity**: Moderate
- **Triage Assessment**: **GENUINE DEFECT**
- **Verified Reproduction**:
  - Warmed `apply_transform(m, p)` allocates 160 bytes per point from `[p[1], p[2], one(T)]` and `m * ...`.
- **Remediation**:
  - Mark `@inline`, widen `m::Matrix{T}` to `m::AbstractMatrix{T}`, verify `size(m) == (3, 3)` (`DimensionMismatch` otherwise), and compute `x` and `y` directly via `muladd` using `p.x` and `p.y` (per `AGENTS.md` 2D coordinate access conventions).

### Finding 7: Vector Splatting and `O(n)` Allocations in `union_bbox(objects::AbstractVector)`
- **Location**: `src/BBoxes.jl:73`
- **Severity**: Moderate
- **Triage Assessment**: **GENUINE DEFECT**
- **Verified Reproduction**:
  - `union_bbox(boxes)` allocates 111,744 bytes for 100 boxes and 13,234,680 bytes (~13.2 MB) for 10,000 boxes due to splatting `objects...` and allocating `map(bbox, objects)`.
- **Remediation**:
  - Implement `union_bbox(objects::AbstractVector)` directly as a single-pass streaming reduction via a concrete-typed helper `_union_bbox_iter(objects, result::BBox{D,T}, state)`, calling `bbox(object)` once per element and promoting `result` only when an element widens `T`.

---

## 2. Implementation & Quality Verification Plan

1. **Apply Targeted Fixes**:
   - Edit `src/Points.jl`, `src/Segments.jl`, `src/Transforms2D.jl`, and `src/BBoxes.jl` with minimal, well-factored changes keeping all functions compact and nesting depth `<= 4`.
   - Preserve all public APIs and exported names.
2. **Add Regression Unit Tests**:
   - `test/runtests.jl`: add `@inferred` and type-stability tests for `dist_subspace` (`Int`, `Float32`, `Float64`, `Rational`), `rand_point(Val(D))` and `rand_gaussian(Val(D))` (zero allocations + inference), `nn_point` (`Int`, `Float32`, `Rational` degenerate and nondegenerate), `bisection_point` (`Int`, `Float32`, `Rational` parallel and nonparallel), `dist_segment_segment` (zero allocations + small nonparallel segments + `Float32`/`Int`), and `apply_transform` (zero allocations + `SMatrix`/`Matrix` + dimension check).
   - `test/test_mat_pnt_seq.jl`: add `union_bbox` tests for large vectors (bounded allocations, concrete inference), heterogeneous element types, empty input (`ArgumentError`), and mismatched dimensions (`DimensionMismatch`).
3. **Update `CHANGELOG.md`**:
   - Document user-visible fixes and performance improvements under `## [Unreleased]`.
4. **Verify Quality Gate & Build**:
   - Run `tools/build` and `tools/gate` (`Pkg.test()`), ensuring 100% pass with no skipped assertions.
   - Author `reviews/007_summary.md` and commit all changes.
