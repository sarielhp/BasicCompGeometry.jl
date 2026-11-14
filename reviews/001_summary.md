# Code Review Remediation Summary #001 (Lens: performance)

- **Review Document**: `reviews/001_performance.md`
- **Architectural Plan**: `reviews/001_performance_plan.md`
- **Date**: 2026-11-13
- **Focus Area**: Performance Audit (Monotone & Sweep Triangulation)
- **Status**: Complete & Verified

---

## 1. Executive Summary

All 4 issues reported in `reviews/001_performance.md` were evaluated and confirmed to be genuine performance defects. Targeted, minimal, and idiomatic fixes were implemented across `src/MonotoneTriangulation.jl` and `src/SweepTriangulation.jl`. Dedicated regression tests covering correctness, exact CCW orientation, and memory allocations were added to `test/test_monotone_triangulation.jl`.

### Performance Highlights
- **`triangulate_monotone` on alternating polygons**:
  - `FloatingKernel`: Allocations reduced from **4,040 allocations (291 KiB)** down to **38 allocations (151 KiB)** (~99% reduction in allocations, ~3x faster).
  - `ExactKernel`: Allocations reduced from **280,170 allocations (9.1 MiB)** down to **173 allocations (155 KiB)**; execution time dropped from **6.48 ms down to 0.067 ms (~97x speedup)** due to eliminating redundant `orient_sign` predicate calls.
- **`_check_distinct_vertices` on 20,000 vertices**:
  - Memory allocated reduced from **547 KiB down to 234 KiB (~57% reduction)**.
  - Execution time improved from **8.05 ms (87.5% GC overhead) down to 0.65 ms (0% GC overhead)** (~12x faster).
- **`_sweep_partition_diagonals`**:
  - Eliminated $O(n)$ heap allocation for `prev_v` array via inline modular arithmetic (`v == 1 ? n : v - 1`).
  - Added `sizehint!(diags, n - 3)` to avoid dynamic reallocation during plane sweep.

---

## 2. Detailed Findings & Mitigations

### Finding 1: Stack Buffer Reallocation in `triangulate_monotone`
- **Location**: `src/MonotoneTriangulation.jl:163-176`
- **Status**: **Mitigated**
- **Root Cause**: When a sweep event transitioned between upper and lower boundary chains, `stack = Int[order[j-1], v]` discarded the previously allocated stack and allocated a new heap vector, triggering $O(k)$ allocations and buffer resizings on alternating polygons.
- **Remediation**:
  - Pre-allocated stack capacity: `stack = Int[]; sizehint!(stack, k); push!(stack, order[1], order[2])`.
  - Reused existing buffer across chain alternations: replaced allocation with `pop!(stack); push!(stack, order[j-1], v)`.
- **Verification**:
  - Regression unit test in `test/test_monotone_triangulation.jl` verifying that `triangulate_monotone` on a 2002-vertex alternating polygon completes with under 200,000 bytes allocated (154,704 bytes measured vs. ~291,570 bytes previously).

### Finding 2: Redundant Geometric Predicate Evaluation via `make_ccw`
- **Location**: `src/MonotoneTriangulation.jl:155-161, 178-185`
- **Status**: **Mitigated**
- **Root Cause**: `push!(triangles, make_ccw(...))` re-evaluated `orient_sign` on every triangle produced during same-chain reflex steps, cross-chain steps, and terminal steps. Because the orientation is deterministically known from chain classification (`is_upper`), calling `make_ccw` doubled predicate calculations in the inner triangulation loop.
- **Remediation**:
  - Removed local closure `make_ccw`.
  - Emitted counterclockwise triangles directly based on `is_upper`:
    - Cross-chain: `push!(triangles, is_upper[v] ? (v, s_next, s_top) : (v, s_top, s_next))`
    - Same-chain: `push!(triangles, is_upper[v] ? (a, v, q) : (a, q, v))`
    - Terminal: `push!(triangles, is_upper[s_top] ? (v_last, s_top, s_next) : (v_last, s_next, s_top))`
- **Verification**:
  - Regression unit test validating that 100% of emitted triangles are strictly CCW under both `FloatingKernel` and `ExactKernel`.
  - ExactKernel allocation dropped from 280k to 173 allocations, with a 97x runtime improvement.

### Finding 3: Intermediate Array Allocation in `_check_distinct_vertices`
- **Location**: `src/MonotoneTriangulation.jl:421-432`
- **Status**: **Mitigated**
- **Root Cause**: `_check_distinct_vertices` created an intermediate `Vector{Tuple{T, T}}` via a comprehension `[(vertices[i][1], vertices[i][2]) for i in 1:n]` solely to pass to `sortperm`.
- **Remediation**:
  - Extracted points vector if input is `AbsPntSeq` via `Points(vertices)`.
  - Replaced intermediate tuple array allocation with direct `sortperm(pts; alg = MergeSort, by = p -> (p[1], p[2]))` for `AbstractVector` inputs, and `sort!(collect(1:n); alg = MergeSort, by = i -> (vertices[i][1], vertices[i][2]))` for arbitrary indexable sequences.
- **Verification**:
  - Regression unit tests on 20,000 vertices demonstrating memory reduction from 547 KiB to 234 KiB.
  - Regression tests verifying duplicate detection throws `ArgumentError` on both `Vector{Point{2, Float64}}` and `PntSeq{2, Float64}`.

### Finding 4: Heap-Allocated `prev_v` Vector and Unhinted `diags` in `_sweep_partition_diagonals`
- **Location**: `src/SweepTriangulation.jl:268, 289`
- **Status**: **Mitigated**
- **Root Cause**: `prev_v = [mod1(i - 1, n) for i in 1:n]` allocated an $n$-element vector whose elements were accessed only once per event. `diags = Tuple{Int,Int}[]` was unhinted and reallocated repeatedly during diagonal insertion.
- **Remediation**:
  - Replaced `prev_v[v]` with an inline cyclic index calculation `e_in = v == 1 ? n : v - 1`.
  - Pre-allocated diagonal capacity with `sizehint!(diags, n - 3)`.
- **Verification**:
  - Regression tests verifying that `triangulate_sweep` and `_sweep_partition_diagonals` produce valid monotone partitions and triangulation across circular and star polygons.

---

## 3. Invariants & Quality Standards

- **Triage First**: All 4 defects were triaged, assessed as genuine, and fully resolved.
- **Targeted Fixes & Regression Coverage**: Fixes were confined strictly to the reported performance issues. All 4 fixes have dedicated regression unit tests in `test/test_monotone_triangulation.jl`.
- **Quality Gate**:
  - Built and precompiled cleanly with `tools/build`.
  - Quality gate `tools/gate` executed and passed all tests cleanly with zero regressions.
- **Language Sizing & Complexity Tiers**:
  - Maximum nesting depth in modified functions is $\le 3$ (well below the maximum of 4).
  - Code size in `triangulate_monotone` was reduced by removing `make_ccw`.
  - Public APIs, exported names, and return contracts remain 100% stable and unchanged.
- **Changelog**: Changes documented under `## [Unreleased]` in `CHANGELOG.md`.
