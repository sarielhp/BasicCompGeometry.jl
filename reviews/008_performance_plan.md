# Architectural Plan: Performance Remediation Cycle #008

## 1. Triage of Findings (`reviews/008_performance.md`)

All 9 findings from `reviews/008_performance.md` were inspected against the source code and confirmed as **genuine defects**. None are false positives or intentional design choices.

| # | Severity | Location | Summary of Defect | Disposition |
|---|---|---|---|---|
| 1 | Major | `src/Arrangements.jl` (`vertical_decomposition`) | `new_active = _ActiveTrap[]` uses abstract `UnionAll` element type, forcing dynamic dispatch and field boxing across the slab sweep | Fix |
| 2 | Major | `src/Arrangements.jl` (`_tc_link_neighbors!`) | $O(n^2) = O(k^4)$ all-pairs nested loop + repeated `any(...)` scan over `vert_line_x` inside inner loop | Fix |
| 3 | Major | `src/Arrangements.jl` (`_vd_walls`) | Linear `findfirst` scan on sorted `unique_x` ($O(k^4)$ total) and four $O(k)$ temporary vector allocations per vertex | Fix |
| 4 | Major | `src/Bundles.jl` (`_unique_points`) | $O(m^2)$ all-pairs scan `any(q -> dist_sq(p, q) <= tol_sq, res)` | Fix |
| 5 | Moderate | `src/Arrangements.jl` (`_vd_slab_lines`, `_vd_find_active`, `vertical_decomposition`) | Per-slab allocations (`slab_lines`, `unique_slab`, `bounds`) and $O(k^2)$-per-slab linear `_vd_find_active` + `deleteat!` | Fix |
| 6 | Moderate | `src/PntSeqs.jl` (`at`, `sample_uniformly`) | `sample_uniformly` calls `at` (and thus `prefix_lengths`) $n$ times ($O(nm)$ time/allocations); `at` returns `Union{Point{D,T}, Point{D,Float64}}` when `T != Float64` | Fix |
| 7 | Moderate | `src/PntSeqHierarchy.jl` (`SimplificationHierarchy`, `_fill!`, `_p_extract!`) | `plt::Palette{T}` lacks `R<:SplitRule` parameter (abstract field); `_fill!` and `_p_extract!` return `0.0` instead of `zero(T)` | Fix |
| 8 | Moderate | `src/PntSeqs.jl` (`dist(::AbsPntSeq, i, j)`) | Returns `0.0::Float64` when `i == j` instead of `zero(float(T))`, causing `Union{Float32, Float64}` | Fix |
| 9 | Moderate | `src/LocalIncidence.jl` (`vertical_ray_segment`, `first_positive_vertical_contact`) | Unconditional `_exact` (`BigInt`/`Rational{BigInt}`) widening before `x`-span rejection | Fix |

---

## 2. Targeted Remediation Plan

### 2.1 `src/Bundles.jl` — `_unique_points` (Finding 4)
- Sort point indices by `(p.x, p.y)` via `sortperm`, pre-size `res` with `sizehint!`, and scan backwards in `res` stopping early as soon as `(p.x - q.x) > abs_tol`.
- Generalize signature parameter from `Point{2,Float64}` to `Point{2,T}` to keep callers in `Arrangements.jl` type-generic.

### 2.2 `src/Arrangements.jl` — Vertical Decomposition & Neighbor Linking (Findings 1, 2, 3, 5)
1. **Concrete double-buffered active trapezoids (Finding 1)**: Preallocate `active_traps = _ActiveTrap{R_type}[]` and `new_active = _ActiveTrap{R_type}[]` outside the slab loop, reuse `new_active` via `empty!`, construct `_ActiveTrap{R_type}(...)`, and swap buffers at the end of each slab.
2. **Binary search & zero-allocation wall bounds (Finding 3)**: Add `_nearest_x_idx(unique_x, x_target, tol)` using `searchsortedfirst` on sorted `unique_x`. In `_vd_walls`, replace `findfirst` with `_nearest_x_idx`, and replace the four temporary vectors (`non_vert`, `y_vals`, `above_y`, `below_y`) with a direct running `wall_top`/`wall_bot` loop over `lines`.
3. **Scratch buffer reuse & two-pointer slab matching (Finding 5)**:
   - Reuse preallocated `slab_lines`, `unique_slab`, and `bounds` across slabs using `_vd_slab_lines!`.
   - Extract `_vd_advance_slab!` to match bottom-to-top `active_traps` against bottom-to-top `bounds` pairs with an $O(k)$ two-pointer scan, eliminating `_vd_find_active` and `deleteat!`, and retiring non-extending trapezoids directly into `completed`.
4. **Indexed neighbor linking (Finding 2)**: In `_tc_link_neighbors!`, sort `vert_line_x` once and binary-search `xr - tol` before the inner loop, and sort trapezoid indices by `x_left` (`order`, `x_lefts`) so each trapezoid only inspects candidates with `x_left in [xr - tol, xr + tol]`.

### 2.3 `src/PntSeqs.jl` — Type Stability & $O(m+n)$ Resampling (Findings 6, 8)
1. **`dist(::AbsPntSeq{D,T}, i::Int, j::Int)` (Finding 8)**: Return `zero(float(T))` when `i == j` so return type matches `dist(P.pnts[i], P.pnts[j])`.
2. **`at` and `sample_uniformly` (Finding 6)**:
   - In `at(pnt_seq::AbsPntSeq{D,T}, t::Real)`, standardize all return branches to `Point{D,float(T)}` (preserving `Float32` for `Float32` sequences, `Float64` for `Float64`/`Int` sequences, and `BigFloat` for `BigFloat` sequences).
   - In `sample_uniformly(P::AbsPntSeq{D,T}, n::Int)`, compute `prefix_lengths(P)` once and advance a monotone segment pointer `seg` in a single $O(m + n)$ pass into a preallocated `Vector{Point{D,float(T)}}`, preserving `PntSeq{D,Float32}` output on `Float32` curves while fixing `InexactError` on integer sequences.

### 2.4 `src/PntSeqHierarchy.jl` — Concrete Palette Field & Base-Case Types (Finding 7)
- Parameterize `SimplificationHierarchy{D,T,R<:SplitRule,M<:HierarchyMeasure}` on `R` so `plt::Palette{T,R}` is concretely typed.
- Replace `0.0` base-case returns in `_fill!` and `_p_extract!` with `zero(eltype(errs))` and `zero(eltype(pal))`.

### 2.5 `src/LocalIncidence.jl` — Early Span Rejection Before `BigInt` Promotion (Finding 9)
- Add `@inline _check_exact_scalar(c)` (in `src/ChainIndex.jl`, where `_ExactScalar` is defined) so non-exact inputs (e.g., `Float64`) still throw `ArgumentError` unconditionally while concrete `Int`/`Rational` inputs incur zero runtime cost.
- In `vertical_ray_segment`, check exactness via `_check_exact_scalar` and reject out-of-span `x` (and out-of-ray-side `y`) before calling `_exact(...)` to allocate `BigInt` / `Rational{BigInt}`.

---

## 3. Verification & Regression Tests
- Add regression tests in:
  - `test/test_bundles.jl`: `_unique_points` / `vertices(...; unique=true)` correctness and subquadratic scaling.
  - `test/test_arrangements.jl`: `vertical_decomposition` & `trapezoid_chains` correctness on crossing/vertical/coincident lines and allocation/scaling bounds.
  - `test/runtests.jl`: `@inferred` type stability for `dist`, `at`, and `sample_uniformly` across `Float32`/`Float64`/`Int`, and allocation bound for `sample_uniformly`.
  - `test/test_pnt_seq_hierarchy.jl`: concrete `isconcretetype(fieldtype(typeof(h), :plt))` and `@inferred` checks for `_fill!`, `_p_extract!`, and `palette_extract` on `Float32`.
  - `test/test_chain_incidence.jl`: error check on out-of-span inexact coordinates and bounded allocations for `first_positive_vertical_contact` on non-spanning integer segments.
- Run `tools/build` and `tools/gate` to verify all unit tests and precompilation pass cleanly.
- Update `CHANGELOG.md` under `[Unreleased]` with user-visible improvements.
