# Architectural Plan: Performance Review #001 Remediation

- **Review Document**: `reviews/001_performance.md`
- **Focus Area**: Performance Audit (Monotone & Sweep Triangulation)
- **Target Files**:
  - `src/MonotoneTriangulation.jl`
  - `src/SweepTriangulation.jl`
  - `test/test_monotone_triangulation.jl`
  - `CHANGELOG.md`
- **Output Artifacts**:
  - `reviews/001_performance_plan.md`
  - `reviews/001_summary.md`

---

## 1. Triage of Findings

All 4 reported findings from `reviews/001_performance.md` have been triaged against the codebase and confirmed as genuine performance defects:

### Finding 1: Stack Buffer Reallocation in `triangulate_monotone`
- **Location**: `src/MonotoneTriangulation.jl:173-175`
- **Severity**: Moderate
- **Triage Assessment**: **GENUINE DEFECT**
- **Analysis**:
  In `triangulate_monotone`, during the sweep across an x-monotone polygon, when current vertex `v` lies on the opposite chain from `last(stack)`, the inner loop drains the stack down to 1 element, pops the final element, and then executes:
  ```julia
  stack = Int[order[j-1], v]
  ```
  This discards the previously allocated vector buffer and allocates a fresh 2-element `Vector{Int}` on the heap. On alternating polygons (combs, staircases, zigzags), this triggers $O(k)$ heap allocations and buffer resizings in the inner triangulation loop.
- **Resolution**:
  Pre-size `stack` with `sizehint!(stack, k)` upon creation. On chain alternation, pop the last element and push `order[j-1]` and `v` into the existing stack without discarding the buffer:
  ```julia
  pop!(stack)
  push!(stack, order[j-1], v)
  ```

### Finding 2: Redundant Geometric Predicate Re-evaluation via `make_ccw`
- **Location**: `src/MonotoneTriangulation.jl:155-161, 178-185`
- **Severity**: Moderate
- **Triage Assessment**: **GENUINE DEFECT**
- **Analysis**:
  In `triangulate_monotone`, triangle orientations are deterministically known from chain classification (`is_upper`):
  1. For cross-chain vertices: When `is_upper[v]` is true, vertex `v` is on the upper chain and the popped stack elements `s_next`, `s_top` are on the lower chain in increasing $x$-order ($x(s_{next}) < x(s_{top}) \le x(v)$), forming CCW triangle `(v, s_next, s_top)`. When `is_upper[v]` is false, `(v, s_top, s_next)` is CCW.
  2. For same-chain vertices: `turn = orient(kernel, vertices[a], vertices[q], vertices[v])` is already evaluated to test reflexivity. If `is_upper[v]` is true, `turn < 0` (CW), meaning CCW order is `(a, v, q)`. If `is_upper[v]` is false, `turn > 0` (CCW), meaning CCW order is `(a, q, v)`.
  3. For the terminal vertex `v_last`: When `is_upper[s_top]` is true, CCW order is `(v_last, s_top, s_next)`; otherwise `(v_last, s_next, s_top)`.
  Currently, `push!(triangles, make_ccw(...))` is called in all three cases, re-evaluating `orient_sign(kernel, ...)` on every emitted triangle. This doubles the predicate evaluations in hot loops, causing severe slowdowns under `ExactKernel()`.
  Note: This exact direct orientation logic is already proven and used in `src/SeidelTriangulation.jl:852-892`.
- **Resolution**:
  Remove `make_ccw` and emit CCW triangles directly:
  - Cross-chain: `push!(triangles, is_upper[v] ? (v, s_next, s_top) : (v, s_top, s_next))`
  - Same-chain: `push!(triangles, is_upper[v] ? (a, v, q) : (a, q, v))`
  - Terminal: `push!(triangles, is_upper[s_top] ? (v_last, s_top, s_next) : (v_last, s_next, s_top))`

### Finding 3: Intermediate Array Allocation in `_check_distinct_vertices`
- **Location**: `src/MonotoneTriangulation.jl:421-432`
- **Severity**: Moderate
- **Triage Assessment**: **GENUINE DEFECT**
- **Analysis**:
  `_check_distinct_vertices` allocates an intermediate `Vector{Tuple{T, T}}` of length $n$ via a list comprehension `[(vertices[i][1], vertices[i][2]) for i in 1:n]` purely to pass to `sortperm`. For large polygons validated in `_validate_ccw_polygon` (called by ear-clipping and Seidel algorithms), this creates unnecessary heap allocations and GC pressure.
- **Resolution**:
  Pass `by = p -> (p[1], p[2])` to `sortperm` directly on `vertices`:
  ```julia
  order = sortperm(vertices; alg = MergeSort, by = p -> (p[1], p[2]))
  ```
  This completely removes the intermediate vector allocation.

### Finding 4: Heap-Allocated `prev_v` Vector and Unhinted `diags` in `_sweep_partition_diagonals`
- **Location**: `src/SweepTriangulation.jl:268, 289`
- **Severity**: Moderate
- **Triage Assessment**: **GENUINE DEFECT**
- **Analysis**:
  In `_sweep_partition_diagonals`, an $n$-element vector `prev_v = [mod1(i - 1, n) for i in 1:n]` is allocated on the heap. Its only use in the entire file is at line 292: `e_in = prev_v[v]`.
  Since `mod1(v - 1, n)` for $v \in 1:n$ is simply `v == 1 ? n : v - 1`, this is an $O(1)$ arithmetic operation.
  Furthermore, `diags = Tuple{Int,Int}[]` is initialized without `sizehint!(diags, n - 3)`, causing repeated dynamic array growth as partition diagonals are discovered.
- **Resolution**:
  Replace `prev_v` with an `@inline _prev_v(i::Int, n::Int) = i == 1 ? n : i - 1` helper (or direct ternary expression `v == 1 ? n : v - 1`).
  Add `sizehint!(diags, n - 3)` to pre-allocate diagonal storage.

---

## 2. Invariant & Standards Compliance

1. **Cognitive Complexity & Function Sizing**:
   - `triangulate_monotone`: Nesting depth remains $\le 3$. Lines of code reduced by removing `make_ccw`.
   - `_check_distinct_vertices`: Nesting depth 2, single loop, no complexity increase.
   - `_sweep_partition_diagonals`: Nesting depth unaffected, removes 1 heap allocation, sizing unchanged.
2. **Performance of Hot Loops**:
   - Eliminates redundant predicate calls in `triangulate_monotone` ($k-2$ predicates saved per polygon).
   - Eliminates $O(k)$ vector reallocations in `triangulate_monotone`.
   - Eliminates $O(n)$ vector allocation in `_check_distinct_vertices`.
   - Eliminates $O(n)$ vector allocation in `_sweep_partition_diagonals`.
3. **API & Export Stability**:
   - All function signatures, exported symbols, and return types remain strictly identical.
   - Triangle CCW orientation guarantee is preserved.

---

## 3. Implementation Steps

1. **Step 1: Modify `src/MonotoneTriangulation.jl`**:
   - Pre-size `stack` with `sizehint!(stack, k)` in `triangulate_monotone`.
   - Update chain alternation to reuse `stack` via `pop!(stack)` and `push!(stack, order[j-1], v)`.
   - Update triangle emission in cross-chain, same-chain, and terminal branches to emit CCW tuples directly based on `is_upper`. Remove unused `make_ccw`.
   - In `_check_distinct_vertices`, replace comprehension with `sortperm(vertices; alg = MergeSort, by = p -> (p[1], p[2]))`.
2. **Step 2: Modify `src/SweepTriangulation.jl`**:
   - In `_sweep_partition_diagonals`, remove `prev_v = [mod1(i - 1, n) for i in 1:n]`.
   - Define `@inline _prev_v(i::Int, n::Int) = i == 1 ? n : i - 1` and replace `prev_v[v]` with `_prev_v(v, n)`.
   - Add `sizehint!(diags, n - 3)` after initializing `diags`.
3. **Step 3: Add Regression Unit Tests in `test/test_monotone_triangulation.jl`**:
   - Add unit tests verifying:
     - Allocation count and time on alternating chain polygons (comb/zigzag) in `triangulate_monotone`.
     - Direct orientation correctness under both `FloatingKernel` and `ExactKernel`.
     - Allocation reduction in `_check_distinct_vertices` on large polygons.
     - Allocation reduction in `_sweep_partition_diagonals` / `triangulate_sweep`.
4. **Step 4: Quality Gate & Validation**:
   - Run `julia --project=. test/test_monotone_triangulation.jl`.
   - Run `tools/gate` (full suite of 156k+ tests) and `tools/build`.
   - Update `CHANGELOG.md` under `[Unreleased]`.
   - Write `reviews/001_summary.md`.
