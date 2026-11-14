# Post-Work Summary: Performance Review #006

- **Review Document**: `reviews/006_performance.md`
- **Plan Document**: `reviews/006_performance_plan.md`
- **Focus Area**: `performance`
- **Status**: Completed & Verified

---

## 1. Executive Summary

All 4 findings in `reviews/006_performance.md` were triaged, confirmed as genuine defects, remediated with minimal, type-stable fixes, and verified with targeted regression tests and the full `tools/gate` + `tools/build` checks:

- `LongestConvexSubset.compute_largest_convex_subset` (60 points): allocations dropped from **192,478,832 B** (~0.405 s) to **173,584 B** (~0.00182 s), and both `LCSContext.PS` and `RotateOrder.PS` now have concrete field types.
- `MVBB.approx_diam` and `MVBB.approx_diam_subspace`: priority-queue keys now infer concretely as `Tuple{Node{D,T,S,V},Node{D,T,S,V}}` (`isconcretetype(keytype(pq)) == true`); on 2,000 3D points with `ε = 0.005`, `approx_diam` runtime dropped from **0.0422 s** (174,024 B) to **0.000193 s** (110,696 B) and `approx_diam_subspace` dropped from **0.0176 s** (243,336 B) to **0.00558 s** (160,136 B).
- `ConvexHull3D.convex_hull_3d` (400 sphere points): allocations dropped from **4,639,224 B** to **2,630,768 B** (~43% reduction) by eliminating per-horizon-edge conflict concatenation and `unique!` churn.

---

## 2. Detailed Findings, Mitigations & Verification

### Finding 1: Abstract `PS::AbsPntSeq{2,T}` Fields in `LCSContext` and `RotateOrder`
- **Location**: `src/LongestConvexSubset.jl:45-58` (`LCSContext`, `RotateOrder`)
- **Status**: **Mitigated**
- **Mitigation**:
  - Parameterized `mutable struct LCSContext{T,P<:AbsPntSeq{2,T}}` and `struct RotateOrder{T,P<:AbsPntSeq{2,T}} <: Base.Order.Ordering` with `PS::P`.
  - Updated `pt_init` and `lcs_init` to construct `RotateOrder(PS, PS[loc])` and `LCSContext(PS, pt, lens)` so Julia preserves the concrete sequence subtype `P`.
- **Verification**:
  - `isconcretetype(fieldtype(typeof(lct), :PS))` and `isconcretetype(fieldtype(typeof(ro), :PS))` pass in `test/test_longest_convex_subset.jl` (`"Concrete context types and context reuse"`).

### Finding 2: Repeated `O(n^2)` Context Initialization and Sorting in `compute_largest_convex_subset`
- **Location**: `src/LongestConvexSubset.jl:64-93, 334-340` (`pt_init`, `lcs_init`, `_reset_context!`, `compute_largest_convex_subset`)
- **Status**: **Mitigated**
- **Mitigation**:
  - Hoisted `lct_upper = lcs_init(work_ps)` and `lct_lower = lcs_init(ps_reflected)` before the `for target = 1:n` loop.
  - Added `_reset_context!(lct::LCSContext)` to clear `lct.Lens` via `empty!` (retaining allocated capacity) and reset `pt.f_terminal = false` across `lct.turns` before each `compute_dp_table!` pass.
  - Preallocated `tail`, `tvals`, `jump` in `pt_init` and `sizehint!(lens, n)` in `lcs_init`.
- **Verification**:
  - Verified exact output equivalence against per-target `lcs_init` across 25-point and 30-case randomized sweeps, and added an allocation bound test (`@allocated compute_largest_convex_subset(ps) < 2_000_000` on 60 points; measured 175,728 B) in `test/test_longest_convex_subset.jl`.

### Finding 3: Non-Concrete `PriorityQueue{Tuple{Node,Node},T}` in `approx_diam` and `approx_diam_subspace`
- **Location**: `src/MVBB.jl:246, 352` (`approx_diam`, `approx_diam_subspace`, `_init_node_pair_pq`, `_push_child_pairs!`)
- **Status**: **Mitigated**
- **Mitigation**:
  - Narrowed `root = tree.root` via `root === nothing && return pp` and extracted `_init_node_pair_pq(root::N, d0, ::Type{T}) where {N<:Node,T}` to build `PriorityQueue{Tuple{N,N},T}(Base.Order.Reverse)` with concrete `Node{D,T,S,V}` keys.
  - Parameterized `_push_child_pairs!(pq, n1, n2, threshold, max_dist_fn::F) where {F}` to specialize on the bound function/closure.
- **Verification**:
  - Verified `isconcretetype(keytype(pq))`, `isconcretetype(typeof(pq))`, and `@inferred` return types for both `approx_diam` and `approx_diam_subspace` in `test/test_mvbb.jl` (`"MVBB PriorityQueue type stability"`).

### Finding 4: Per-Horizon-Edge Conflict Concatenation & `unique!` Churn in `ConvexHull3D`
- **Location**: `src/ConvexHull3D.jl:183-201, 245-283` (`convex_hull_3d`, `_replace_visible_faces!`, `_assign_face_conflicts!`)
- **Status**: **Mitigated**
- **Mitigation**:
  - Allocated `seen_stamp = zeros(Int, n_in)` once in `convex_hull_3d` with a monotonic integer `stamp` counter passed through `_replace_visible_faces!`.
  - Extracted `_assign_face_conflicts!` to scan `(f.visible_points, f_nonvis.visible_points)` directly and deduplicate candidates in `O(1)` via `seen_stamp[pt] == stamp` without temporary vector concatenation or `unique!` allocation.
- **Verification**:
  - Verified hull validity (`verify_convex_hull_3d == (true, "Hull is valid")`) on 400-point and 2,000-point spherical shells and added an allocation bound test (`@allocated ... < 3_500_000` on 400 sphere points; measured 2.54–2.73 MB) in `test/test_convex_hull_3d.jl`.

---

## 3. Quality Gate & Metrics Verification

- **Build (`./tools/build`)**: Passed (`exit 0`).
- **Full Test Suite (`./tools/gate`, `Pkg.test()`)**: Passed (`156,417 / 156,417` assertions in `BasicCompGeometry.jl`, plus all `Aqua` package quality checks).
- **Sizing & Nesting**: All new and modified functions (`pt_init`, `lcs_init`, `_reset_context!`, `compute_largest_convex_subset`, `_init_node_pair_pq`, `_push_child_pairs!`, `_assign_face_conflicts!`, `_replace_visible_faces!`) are concise and maintain nesting depth `<= 3` (within the `<= 4` limit).
- **Scope & Anti-Oscillation**: Unrelated modules (`src/RMQ.jl`, `src/BBT.jl`, `src/WSPD.jl`, `src/VisibilityPolygon.jl`, `src/MinkowskiSum.jl`, triangulation modules, etc.) were left untouched. Public APIs and exported names are unchanged; changes are recorded in `CHANGELOG.md` under `## [Unreleased]`.
