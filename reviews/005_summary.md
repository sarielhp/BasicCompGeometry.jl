# Code Review Remediation Summary #005 (Lens: performance)

- **Review Document**: `reviews/005_performance.md`
- **Architectural Plan**: `reviews/005_performance_plan.md`
- **Date**: 2026-11-13
- **Focus Area**: Performance Audit (`RMQ.jl`, `BBT.jl`, and `WSPD.jl`)
- **Status**: Complete & Verified

---

## 1. Executive Summary

All 4 findings reported in `reviews/005_performance.md` were triaged, reproduced, and confirmed as genuine performance, type-stability, or scalar-type conversion defects. Targeted, idiomatic fixes were implemented in `src/RMQ.jl`, `src/BBT.jl`, and `src/WSPD.jl`, with regression unit tests added to `test/test_chain_incidence.jl` and `test/runtests.jl`.

### Measured Improvements (Before vs. After)
- `BasicCompGeometry._sparse_levels` on 100,000 `Float64` keys: **82,810,040 B -> 12,554,664 B** (~6.6x reduction).
- `build_rmq(Vector{Float64})`: `Union{BlockRMQ, SparseRMQ}` -> concrete `BlockRMQ{Vector{Float64}, typeof(isless)}` (and `SparseRMQ{...}` when `backend = :sparse`).
- `build_chain_index` + 1,000 `leftmost(w)` queries inside the same function: **43,632 B -> 9,152 B** (with 10,000 queries also allocating **9,152 B**, i.e., **0 B per query**).
- `BBT.get_min_orig_index` on an uncached, fully expanded 2,000-point `BBT.Tree`: **34,160 B -> 0 B**.
- `WSPD.separation` on `WSPair{2, Float32, Float32, ...}`: `Union{Float32, Float64}` -> concrete `Float32`; `approx_diameter(ps32, 0.1)` no longer throws `MethodError` and returns `Float32`.

---

## 2. Detailed Findings, Mitigations & Verification

### Finding 1: `Core.Box` Capture & Missing `less::F` Specialization in `_sparse_levels`
- **Location**: `src/RMQ.jl` (`_sparse_levels`, `_better_index`, `_block_masks`)
- **Status**: **Mitigated**
- **Remediation**:
  - Parameterized `_sparse_levels(keys::V, less::F, minima::Vector{Int}) where {V,F}` (and `_better_index`, `_block_masks` on `keys::V, less::F where {V,F}`) so Julia specializes on the comparator function type `F`.
  - Replaced the `half`-capturing comprehension with a preallocated `row = Vector{Int}(undef, m)` filled by an `@inbounds for i in 1:m` loop.
- **Verification**:
  - Allocation bound check on `_sparse_levels` over 10,000 keys passes in `test/test_chain_incidence.jl`.
  - Measured 100,000-key allocation dropped from 82,810,040 B to 12,554,664 B.

### Finding 2: Value-Dependent Union Return Type & Missing `less::F` Specialization in `build_rmq`
- **Location**: `src/RMQ.jl` (`build_rmq`, `_build_rmq`)
- **Status**: **Mitigated**
- **Remediation**:
  - Extracted `@inline` backend constructors `_build_rmq(keys::V, less::F, ::Val{:blocked})` and `_build_rmq(keys::V, less::F, ::Val{:sparse})`.
  - Updated `build_rmq(keys::V; less::F = isless, backend::Union{Symbol,Val} = Val(:blocked)) where {V<:AbstractVector,F}` with `Base.@constprop :aggressive`, validating `Val{:blocked}` / `Val{:sparse}` and throwing `ArgumentError` for unknown backends.
- **Verification**:
  - `@inferred` tests for default, `:blocked`, `:sparse`, and custom `less` RMQ construction, invalid backend `ArgumentError` checks, concrete `ChainIndex` inference, and identical `1,000` vs `10,000` query loop allocations pass in `test/test_chain_incidence.jl`.

### Finding 3: Dynamic Reflection (`getfield`/`setfield!`) & Generator Allocation in `BBT.get_min_max_orig_index`
- **Location**: `src/BBT.jl` (`get_min_orig_index`, `get_max_orig_index`, `_range_ext_orig_index`)
- **Status**: **Mitigated**
- **Remediation**:
  - Replaced `get_min_max_orig_index` with a shared `@inline` helper `_range_ext_orig_index(tree::Tree, r::UnitRange{Int}, better::F)::Int where {F}` (guarding empty ranges with `ArgumentError` and scanning with `@inbounds`).
  - Rewrote `get_min_orig_index` and `get_max_orig_index` to access `node.min_index` and `node.max_index` directly with local `left, right = node.left, node.right` narrowing.
- **Verification**:
  - `@inferred` and `@allocated == 0` checks on uncached expanded `BBT.Tree` instances pass in `test/runtests.jl` (`BBT (Bounding Box Tree)` testset).

### Finding 4: Type-Parameter & Return-Type Instability in `WSPD`
- **Location**: `src/WSPD.jl` (`separation`, `PD`, `init`, `push_pair!`, `top_refine!`, `top_diam_ub`) and `src/BBT.jl` (`Tree_init_inner`)
- **Status**: **Mitigated**
- **Remediation**:
  - Updated `separation(pair::WSPair{D,T,S,V})::S` to `@inline`, returning `S(MAX_SEP)` when `pair.dist <= zero(S)`.
  - Parameterized `mutable struct PD{D,T,S,V,P<:AbsPntSeq{D,T}}` on concrete `pnts::P` and updated `push_pair!`, `top_refine!`, and `top_diam_ub` signatures (plus `zero(S)` comparisons in `top_refine!`).
  - Updated `WSPD.init(pnts::P, sep::Real) where {D,T,P<:AbsPntSeq{D,T}}` to derive `S = typeof(tree.root.diam)` and `V = typeof(Points(pnts))`, converting `sep` via `S(sep)`.
  - Updated `BBT.Tree_init_inner` to compute `S = typeof(dist(zero(Point{D,T}), zero(Point{D,T})))` without a runtime `isempty(pnts)` branch so both `BBT.Tree_init` and `WSPD.init` infer a concrete `S` for `Float32`.
- **Verification**:
  - `@inferred WSPD.init`, `isconcretetype(typeof(W32))`, `@inferred WSPD.separation`, `@inferred WSPD.expand!`, and `approx_diameter` on `Float32` point sequences pass in `test/runtests.jl`.

---

## 3. Quality Gate & Metrics Verification

- **Build (`tools/build`)**: Passed cleanly (`exit 0`).
- **Full Test Suite (`tools/gate`, `Pkg.test()`)**: Passed (`156,403 / 156,403` assertions in `BasicCompGeometry.jl`, plus all `Aqua` package checks).
- **Sizing & Nesting**: All new and modified functions in `src/RMQ.jl`, `src/BBT.jl`, and `src/WSPD.jl` are compact (`<= 19` lines each) with shallow nesting (`<= 3`, below the `<= 4` limit).
- **Scope & Anti-Oscillation**: Unrelated modules (`src/Diameter.jl`, `src/MetricSpace.jl`, `src/VisibilityPolygon.jl`, `src/MinkowskiSum.jl`, triangulation modules, etc.) were left untouched. Public APIs and exports remain unchanged; user-visible changes are recorded in `CHANGELOG.md` under `## [Unreleased]`.
