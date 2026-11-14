# Architectural Plan: Performance Review #005 Remediation

- **Review Document**: `reviews/005_performance.md`
- **Focus Area**: Performance Audit (Type Stability, Specialization, and Allocations in `RMQ.jl`, `BBT.jl`, and `WSPD.jl`)
- **Target Files**:
  - `src/RMQ.jl`
  - `src/BBT.jl`
  - `src/WSPD.jl`
  - `test/test_chain_incidence.jl`
  - `test/test_bbt_ann.jl`
  - `test/runtests.jl`
  - `CHANGELOG.md`
- **Output Artifacts**:
  - `reviews/005_performance_plan.md`
  - `reviews/005_summary.md`

---

## 1. Triage of Findings

All 4 findings in `reviews/005_performance.md` have been triaged and reproduced against the current codebase via `Base.return_types`, `@allocated`, and `Float32` test invocations. All 4 are confirmed as **genuine defects**:

### Finding 1: `Core.Box` Capture & Missing `less::F` Specialization in `_sparse_levels`
- **Location**: `src/RMQ.jl:36-50`
- **Severity**: Major
- **Triage Assessment**: **GENUINE DEFECT**
- **Verified Reproduction**:
  - Calling `BasicCompGeometry._sparse_levels(keys, isless, collect(1:100_000))` on a 100,000-element `Vector{Float64}` allocates **82,810,040 bytes** after warmup.
  - The comprehension `[ _better_index(keys, less, previous[i], previous[i+half]) for i = 1:length(minima)-2half+1 ]` captures `half` across `half *= 2` reassignments, wrapping `half` in a `Core.Box` (`Any`) and producing an untyped `Base.Generator`. Additionally, `_sparse_levels(keys, less, minima::Vector{Int})` does not parameterize `less::F`, so Julia's specialization heuristic skips specializing on `less` and dynamically dispatches `_better_index` on every element.
- **Remediation**:
  - Parameterize `_sparse_levels(keys::V, less::F, minima::Vector{Int}) where {V,F}` (as well as `_better_index` and `_block_masks` on `keys::V, less::F where {V,F}`).
  - Replace the comprehension with a preallocated `row = Vector{Int}(undef, m)` and a tight `@inbounds for i in 1:m` loop.

### Finding 2: Value-Dependent Union Return Type & Missing `less::F` Specialization in `build_rmq`
- **Location**: `src/RMQ.jl:90-100` (and downstream `src/ChainIndex.jl:77-88`)
- **Severity**: Major
- **Triage Assessment**: **GENUINE DEFECT**
- **Verified Reproduction**:
  - `Base.return_types(BasicCompGeometry.build_rmq, (Vector{Float64},))` infers `Union{BlockRMQ{Vector{Float64}, typeof(isless)}, SparseRMQ{Vector{Float64}, typeof(isless)}}`.
  - Downstream `build_chain_index` builds four RMQs, inferring a 16-way `UnionAll` over `ChainIndex{V,A,B,C,D}`, which exceeds Julia's union-splitting threshold. Calling `build_chain_index` and querying `leftmost(w)` 1,000 times inside the same function allocates **43,632 bytes**.
- **Remediation**:
  - Extract `@inline` backend helpers `_build_rmq(keys::V, less::F, ::Val{:blocked})` and `_build_rmq(keys::V, less::F, ::Val{:sparse})`.
  - Annotate `build_rmq(keys::V; less::F = isless, backend::Union{Symbol,Val} = Val(:blocked)) where {V<:AbstractVector,F}` with `Base.@constprop :aggressive`, validating `b = backend isa Val ? backend : Val(backend)` is `Val{:blocked}` or `Val{:sparse}` before dispatching to `_build_rmq(keys, less, b)`.
  - Preserve full backward compatibility for callers passing `backend = :blocked`, `backend = :sparse`, or invalid symbols (throwing `ArgumentError`).

### Finding 3: Dynamic Reflection (`getfield`/`setfield!`) & Generator Allocation in `BBT.get_min_max_orig_index`
- **Location**: `src/BBT.jl:94-129`
- **Severity**: Moderate
- **Triage Assessment**: **GENUINE DEFECT**
- **Verified Reproduction**:
  - Building and fully expanding a `BBT.Tree` over 2,000 `Point2F` points and calling `BBT.get_min_orig_index(tree, tree.root)` once after warmup allocates **34,160 bytes**.
  - `get_min_max_orig_index` uses dynamic `getfield(node, field_name)` / `setfield!(node, field_name, val)` with a runtime `Symbol`, inferring a 7-way `Union` over every field type of `Node`, boxing `Int` on `setfield!`, and allocating a `Base.Generator` at leaves/unexpanded nodes.
- **Remediation**:
  - Replace `get_min_max_orig_index` with a shared `@inline` zero-allocation helper `_range_ext_orig_index(tree::Tree, r::UnitRange{Int}, better::F)::Int where {F}` (guarding `isempty(r)` with `ArgumentError`, then scanning `first(r)` and `(first(r)+1):last(r)` with `@inbounds`).
  - Implement `get_min_orig_index` and `get_max_orig_index` using direct field accesses (`node.min_index`, `node.max_index`) and local bindings `left, right = node.left, node.right` for flow-sensitive narrowing.

### Finding 4: Type-Parameter & Return-Type Instability in `WSPD`
- **Location**: `src/WSPD.jl:31-38, 51-58, 130-143`
- **Severity**: Moderate
- **Triage Assessment**: **GENUINE DEFECT** (plus correctness crash for `Float32` / non-`Float64` `sep`)
- **Verified Reproduction**:
  - `Base.return_types(WSPD.separation, (WSPD.WSPair{2,Float32,Float32,Vector{Point{2,Float32}}},))` infers `Union{Float32, Float64}` because `MAX_SEP = 1e18` is a `Float64` literal.
  - Calling `approx_diameter(ps32, 0.1)` on `Float32` points throws `MethodError: Cannot convert an object of type BBT.Tree{2, Float32, Float32, ...} to an object of type BBT.Tree{2, Float32, Float64, ...}` because `WSPD.init(pnts, sep::S)` binds `S` from `sep` instead of from `tree.root.diam`.
  - `PD{D,T,S,V}` stores `pnts::AbsPntSeq{D,T}` with an abstract field type.
- **Remediation**:
  - Update `separation(pair::WSPair{D,T,S,V})::S` to `@inline`, returning `S(MAX_SEP)` when `pair.dist <= zero(S)` and `pair.max_sides_diam / pair.dist` otherwise.
  - Parameterize `mutable struct PD{D,T,S,V,P<:AbsPntSeq{D,T}}` with `pnts::P` and update internal method signatures (`push_pair!`, `top_refine!`, `top_diam_ub`) accordingly.
  - Update `init(pnts::P, sep::Real) where {D,T,P<:AbsPntSeq{D,T}}` to construct `tree = BBT.Tree_init(pnts)`, derive `S = typeof(tree.root.diam)` and `V = typeof(Points(pnts))`, and convert `sep` via `S(sep)`. Also use `zero(S)` instead of `0.0` in `top_refine!` so `Float32` hot paths stay in `S`, and make `BBT.Tree_init_inner` compute `S = typeof(dist(zero(Point{D,T}), zero(Point{D,T})))` without the runtime `isempty(pnts)` branch so both `BBT.Tree_init` and `WSPD.init` are concretely inferred for `Float32`.

---

## 2. Quality, Architecture & Compatibility Constraints

- **Public API**: Preserve all exported names (`build_rmq`, `range_argmin`, `RMQIndex`, `ChainIndex`, `get_min_orig_index`, `get_max_orig_index`, `WSPair`, `PD`, `WSPD.init`, `approx_diameter`, etc.) and keyword compatibility (`backend = :blocked` / `:sparse`). `WSPD.init(pnts, sep)` now accepts any `sep::Real` matching or convertible to `typeof(tree.root.diam)`, fixing the `MethodError` on `Float32` points and integer `sep`.
- **Code Structure**: Keep every new and modified function small, single-purpose, with nesting depth $\le 4$, and share `_range_ext_orig_index` across `get_min_orig_index` and `get_max_orig_index` to avoid duplicated range-scan loops.
- **Scope & Anti-Oscillation**: Keep edits strictly scoped to `src/RMQ.jl`, `src/BBT.jl`, `src/WSPD.jl`, their corresponding tests, and `CHANGELOG.md`. Leave `src/Diameter.jl`, `src/MetricSpace.jl`, `src/VisibilityPolygon.jl`, `src/MinkowskiSum.jl`, and triangulation modules untouched.

---

## 3. Verification & Regression Test Plan

1. **Regression Unit Tests**:
   - `test/test_chain_incidence.jl`:
     - Verify `@inferred` return types for `build_rmq(keys)`, `build_rmq(keys; backend = :sparse)`, `build_rmq(keys; less = custom_less)`, and `build_chain_index(points)`.
     - Verify `build_rmq(keys; backend = :invalid)` and `build_rmq(keys; backend = Val(:invalid))` throw `ArgumentError`.
     - Verify allocation bounds for `_sparse_levels` and for repeated `leftmost(w)` / `bounds(w)` queries on a `ChainIndex` built inside a function.
   - `test/runtests.jl` (inside `BBT` testset):
     - Verify `@inferred get_min_orig_index(tree, tree.root)` and `@inferred get_max_orig_index(tree, tree.root)`, and zero-byte allocation after warmup on an expanded `BBT.Tree` (freshly built/un-cached tree allocation guard as well).
   - `test/test_bbt_ann.jl` / `test/runtests.jl` (`WSPD` / `Diameter` testsets):
     - Verify `WSPD.init` and `WSPD.expand!` / `WSPD.expand_bichromatic!` work on `Float32` point sequences with `Float64` or `Int` `sep` (`approx_diameter(ps32, 0.1)` succeeds and is `@inferred` / matches `exact_diameter` tolerance).
     - Verify `@inferred WSPD.separation(pair)` returns `Float32` for `Float32` trees (both `dist == 0` and `dist > 0` pairs) and `Float64` for `Float64` trees, and `isconcretetype(typeof(W))` with `typeof(W.pnts)` concrete.
2. **Quality Gates**:
   - Run `tools/build` and `tools/gate` (`Pkg.test()`) after edits and confirm all tests pass.
   - Re-measure `@code_warntype` / `Base.return_types` and `@allocated` on all 4 findings from `reviews/005_performance.md`.
