# Architectural Plan: Performance Review #006 Remediation

- **Review Document**: `reviews/006_performance.md`
- **Focus Area**: Performance Audit (Type Specialization, Redundant Work Elimination, and Conflict-List Deduplication Allocations)
- **Target Files**:
  - `src/LongestConvexSubset.jl`
  - `src/MVBB.jl`
  - `src/ConvexHull3D.jl`
  - `test/test_longest_convex_subset.jl`
  - `test/test_mvbb.jl`
  - `test/test_convex_hull_3d.jl`
  - `CHANGELOG.md`
- **Output Artifacts**:
  - `reviews/006_performance_plan.md`
  - `reviews/006_summary.md`

---

## 1. Triage of Findings

All 4 findings in `reviews/006_performance.md` have been evaluated against the current codebase and verified via Julia type/allocation inspection. All 4 are **genuine defects**:

### Finding 1: Abstract `PS::AbsPntSeq{2,T}` Fields in `LCSContext` and `RotateOrder`
- **Location**: `src/LongestConvexSubset.jl:45-58`
- **Severity**: Moderate
- **Triage Assessment**: **GENUINE DEFECT**
- **Verified Reproduction**:
  - `isconcretetype(fieldtype(typeof(lct), :PS))` and `isconcretetype(fieldtype(typeof(ro), :PS))` both return `false` for `PntSeq{2,Float64}`.
  - Every `Base.Order.lt(o::RotateOrder, a, b)` comparison during angular sorting and every `lct.PS[...]` lookup in the DP inner loop reads an abstract `AbsPntSeq{2,T}` field, causing dynamic dispatch.
- **Remediation**:
  - Parameterize `mutable struct LCSContext{T,P<:AbsPntSeq{2,T}}` and `struct RotateOrder{T,P<:AbsPntSeq{2,T}} <: Base.Order.Ordering` with `PS::P`.
  - Construct `RotateOrder(PS, PS[loc])` in `pt_init` and `LCSContext(PS, pt, lens)` in `lcs_init` so type inference preserves the concrete sequence type `P`.

### Finding 2: Repeated `O(n^2)` Context Initialization and Angular Sorting Inside `compute_largest_convex_subset`
- **Location**: `src/LongestConvexSubset.jl:64-93, 334-340`
- **Severity**: Moderate
- **Triage Assessment**: **GENUINE DEFECT**
- **Verified Reproduction**:
  - `compute_largest_convex_subset` on 60 random 2D points allocates **192,478,832 bytes (~192.5 MB)** and takes **~0.405 s** because `lcs_init(work_ps)` and `lcs_init(ps_reflected)` are rebuilt and re-sorted from scratch for every `target in 1:n`.
  - `PointTurn.tail` and the sorted angular topology depend only on the point set, not on `target`.
- **Remediation**:
  - Hoist `lct_upper = lcs_init(work_ps)` and `lct_lower = lcs_init(ps_reflected)` before the `for target = 1:n` loop.
  - Add a compact `_reset_context!(lct::LCSContext)` helper that clears `lct.Lens` via `empty!` (retaining allocated capacity) and resets `pt.f_terminal = false` on each `PointTurn` before each `compute_dp_table!` pass.
  - Preallocate `tail`, `tvals`, `jump` in `pt_init` (instead of incremental `push!`) and `sizehint!(lens, n)` in `lcs_init`.
  - Verify exact equivalence of results across random and structured inputs, plus drastic allocation/time reduction.

### Finding 3: Non-Concrete `PriorityQueue{Tuple{Node,Node},T}` in `approx_diam` and `approx_diam_subspace`
- **Location**: `src/MVBB.jl:246` and `src/MVBB.jl:352`
- **Severity**: Moderate
- **Triage Assessment**: **GENUINE DEFECT**
- **Verified Reproduction**:
  - Both `isconcretetype(BBT.Node)` and `isconcretetype(Tuple{BBT.Node, BBT.Node})` are `false`.
  - `tree.root` has declared type `Union{Nothing, Node{D,T,S,V}}`, so without narrowing `root !== nothing`, even `N = typeof(root)` fails to infer concretely.
- **Remediation**:
  - Narrow `root = tree.root` with `root === nothing && return pp`, bind `N = typeof(root)` (inferred concretely as `Node{D,T,S,V}`), and construct `PriorityQueue{Tuple{N,N},T}(Base.Order.Reverse)`.
  - Parameterize `_push_child_pairs!(pq, n1, n2, threshold, max_dist_fn::F) where {F}` so the comparator/bound closure in `approx_diam_subspace` specializes cleanly.

### Finding 4: Per-Horizon-Edge Conflict List Concatenation and `unique!` Churn in `ConvexHull3D`
- **Location**: `src/ConvexHull3D.jl:245-283` (lines 272-273)
- **Severity**: Moderate
- **Triage Assessment**: **GENUINE DEFECT**
- **Verified Reproduction**:
  - On 400 points sampled from a sphere, `convex_hull_3d` allocates **4,639,224 bytes (~4.64 MB)**, with every horizon edge allocating `[f.visible_points; f_nonvis.visible_points]` and running `unique!(candidates)`.
- **Remediation**:
  - Allocate a single `seen_stamp = zeros(Int, n_in)` scratch vector once in `convex_hull_3d` and maintain a monotonically increasing integer `stamp` (O(1) per horizon face, avoiding even per-face `empty!` or word-array clearing).
  - Extract `_assign_face_conflicts!` (keeping nesting depth `<= 3` and small single-purpose functions) to iterate `(f.visible_points, f_nonvis.visible_points)` directly, deduplicating via `seen_stamp[pt] == stamp` with zero per-face allocations.

---

## 2. Implementation & Quality Verification Plan

1. **Edit Core Modules**:
   - Update `src/LongestConvexSubset.jl`, `src/MVBB.jl`, and `src/ConvexHull3D.jl` with minimal, targeted changes keeping public APIs and exported names unchanged.
   - Ensure all new and modified functions remain small, single-purpose, and within nesting depth `<= 4`.
2. **Add Regression Unit Tests**:
   - `test/test_longest_convex_subset.jl`: verify concrete `PS` field types in `LCSContext` and `RotateOrder`, correctness equivalence on larger point sets, and bounded allocations for `compute_largest_convex_subset` (confirming hoisting eliminates the `O(n^3)` allocation blowup).
   - `test/test_mvbb.jl`: verify concrete `PriorityQueue` key inference in `approx_diam` and `approx_diam_subspace` and correctness on sphere/cube point sets.
   - `test/test_convex_hull_3d.jl`: verify conflict deduplication correctness (shared visible points across horizon edges) and bounded allocations on spherical shells.
3. **Update `CHANGELOG.md`**:
   - Record performance improvements under `## [Unreleased]` -> `### Changed`.
4. **Run Gates & Checks**:
   - Run `tools/build` and `tools/gate` (`Pkg.test()`), verifying 100% pass with zero skipped assertions.
   - Author `reviews/006_summary.md` and commit all changes cleanly.
