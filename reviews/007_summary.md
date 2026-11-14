# Post-Work Summary: Performance Review #007

- **Review Document**: `reviews/007_performance.md`
- **Plan Document**: `reviews/007_performance_plan.md`
- **Focus Area**: `performance`
- **Status**: Completed & Verified

---

## 1. Executive Summary

All 7 findings in `reviews/007_performance.md` were triaged, confirmed as genuine defects, fixed with minimal, type-stable implementations, and verified with targeted regression tests and the full `tools/build` and `tools/gate` checks:

- `dist_subspace`, `nn_point`, and `bisection_point`: eliminated value-dependent union return types across `Int`, `Float32`, `Float64`, and `Rational{Int}`.
- `rand_point` and `rand_gaussian`: added `Val{D}` overloads for compile-time dimension inference and replaced heap-allocated `randn(T, D)` with stack-allocated `ntuple` construction (0 heap allocations vs. 80 bytes previously for `rand_gaussian(3)`).
- `dist_segment_segment`: replaced heap `Vector`/`Matrix` allocation and `m \ c` solve (592 bytes per query) with zero-allocation scalar Cramer's rule and scale-invariant degeneracy checking.
- `apply_transform`: eliminated augmented-vector heap allocations (160 bytes per point -> 0 bytes) using inline `muladd` with `p.x`/`p.y` and widened `Matrix{T}` to `AbstractMatrix{T}`.
- `union_bbox(objects::AbstractVector)`: eliminated runtime vector splatting and intermediate `map(bbox, objects)` allocation via single-pass streaming reduction (`_union_bbox_iter`), reducing 10,000-box allocations from 13,234,680 bytes to 96 bytes.

---

## 2. Detailed Findings, Mitigations & Verification

### Finding 1: Union Return Type in `dist_subspace(p, q, dir)`
- **Location**: `src/Points.jl:56-62` (`dist_subspace`)
- **Status**: **Mitigated**
- **Mitigation**:
  - Replaced `d_sq > 0 ? sqrt(d_sq) : zero(T)` with `sqrt(max(zero(d_sq), d_sq))` so both nonpositive and positive cases return `typeof(sqrt(d_sq))`.
- **Verification**:
  - `Base.return_types` confirmed `Float64` for `Point{2,Int}` (previously `Union{Float64, Int64}`), and `@inferred` tests pass for `Int`, `Float32`, and `Rational{Int}` in `test/runtests.jl`.

### Finding 2: Uninferred `D` and Heap Allocation in `rand_gaussian` / `rand_point`
- **Location**: `src/Points.jl:194-216` (`rand_gaussian`, `rand_point`)
- **Status**: **Mitigated**
- **Mitigation**:
  - Added `rand_gaussian(::Val{D}, ::Type{T} = Float64)` using `Point{D,T}(ntuple(_ -> randn(T), Val(D)))` and `rand_point(::Val{D}, ::Type{T} = Float64) = rand(Point{D,T})`.
  - Retained `D::Int` overloads for API compatibility, delegating `rand_gaussian(D::Int, ::Type{T})` to `rand_gaussian(Val(D), T)` and specializing on `::Type{T}`.
- **Verification**:
  - Verified concrete `SVector{3,Float64}` inference and `@allocated == 0` for `rand_point(Val(3))`, `rand_gaussian(Val(3))`, and `rand_gaussian(3)` in `test/runtests.jl`.

### Finding 3: Value-Dependent Union Return Type in `nn_point(s_p, s_q, qr)`
- **Location**: `src/Segments.jl:122-135` (`nn_point`)
- **Status**: **Mitigated**
- **Mitigation**:
  - Converted endpoints and query to `Point{D,R}` with `R = typeof(one(T) / one(T))`, used `iszero(d_sq) && return p`, and clamped with `clamp(t, zero(t), one(t))`.
- **Verification**:
  - `@inferred` tests in `test/runtests.jl` pass for degenerate and nondegenerate segments across `Int`, `Float32`, and `Rational{Int}`.

### Finding 4: Union Return Type in `bisection_point(seg, p, q)`
- **Location**: `src/Segments.jl:219-236` (`bisection_point`)
- **Status**: **Mitigated**
- **Mitigation**:
  - Promoted segment and bisector points to `Point{D,R}` (`R = typeof(one(T) / one(T))`), returned `(false, zero(R), sp)` on parallel segments, and used `R <: AbstractFloat ? abs(denom) < R(1e-12) : iszero(denom)` to avoid converting `1e-12` into a large-denominator `Rational`.
- **Verification**:
  - `@inferred` tests in `test/runtests.jl` pass for both parallel and nonparallel segments across `Int`, `Float32`, and `Rational{Int}`.

### Finding 5: Heap Matrix/Vector Allocations in `dist_segment_segment`
- **Location**: `src/Segments.jl:144-198` (`_min_endpoint_segment_dist`, `dist_segment_segment`)
- **Status**: **Mitigated**
- **Mitigation**:
  - Replaced `c`, `m`, and `m \ c` with scalar Cramer's rule over promoted `Point{D,R}` coordinates, extracted `_min_endpoint_segment_dist`, and checked `Δ <= tol * (A * C)` (`tol = R <: AbstractFloat ? max(R(1e-12), 4 * eps(R)) : zero(R)`).
- **Verification**:
  - Verified `@allocated == 0` inside a function, exact `0.0` distance for `1e-4`-length crossing segments, and `@inferred` `Float32` and `Int` results in `test/runtests.jl`.

### Finding 6: Heap Allocations in `apply_transform(m, p)`
- **Location**: `src/Transforms2D.jl:44-55` (`apply_transform`)
- **Status**: **Mitigated**
- **Mitigation**:
  - Widened `m::Matrix{T}` to `m::AbstractMatrix{T}`, checked `size(m) == (3, 3)`, and evaluated `x` and `y` inline via `muladd` using `p.x` and `p.y`.
- **Verification**:
  - Verified `@allocated == 0`, `SMatrix{3,3,Float64}` inference, and `DimensionMismatch` on a `2x2` matrix in `test/runtests.jl`.

### Finding 7: Vector Splatting and `O(n)` Allocations in `union_bbox(objects::AbstractVector)`
- **Location**: `src/BBoxes.jl:73-105, 418` (`_union_bbox_iter`, `union_bbox`, `bound!`)
- **Status**: **Mitigated**
- **Mitigation**:
  - Implemented `union_bbox(objects::AbstractVector)` as a single-pass streaming reduction via `_union_bbox_iter` (calling `bbox` once per element and promoting `BBox{D,T}` only if a wider coordinate type is encountered) and generalized `bound!(bb::BBox{D,T}, other::BBox{D})`.
- **Verification**:
  - Measured 96 bytes allocated on 10,000 `BBox2F` objects (down from 13,234,680 bytes); added inference, allocation bound (`<= 1024` bytes), heterogeneous promotion, `ArgumentError`, and `DimensionMismatch` tests in `test/test_mat_pnt_seq.jl`.

---

## 3. Quality Gate & Metrics Verification

- **Build (`tools/build`)**: Passed (`exit 0`).
- **Full Test Suite (`tools/gate`, `Pkg.test()`)**: Passed (`156,449 / 156,449` assertions in `BasicCompGeometry.jl`, plus all 8 `Aqua` testsets).
- **Sizing & Nesting**: All new and modified functions are concise with nesting depth `<= 2` (within the `<= 4` limit).
- **Scope & Anti-Oscillation**: All unrelated modules and prior review remediations (rounds 004–006) remain untouched; public APIs and exported names are unchanged; changes are documented in `CHANGELOG.md` under `## [Unreleased]`.
