# Post-Work Summary: Performance Remediation Cycle #008

- **Source Report**: `reviews/008_performance.md`
- **Plan**: `reviews/008_performance_plan.md`
- **Date**: 2026-11-14

## Issue-by-Issue Status & Mitigation

| # | Severity | Location | Status | Mitigation & Verification |
|---|---|---|---|---|
| 1 | Major | `src/Arrangements.jl` (`vertical_decomposition`) | Fixed | Replaced `new_active = _ActiveTrap[]` with concrete `_ActiveTrap{R_type}[]` double-buffering (`active_traps, new_active = new_active, active_traps`) and `_ActiveTrap{R_type}(...)` construction. Verified with allocation test on 40-line decomposition (`< 5 MB`, measured ~1.28 MiB) and 160-line benchmark (~23.3 MiB / 0.054 s vs ~2.795 GiB / 0.617 s). |
| 2 | Major | `src/Arrangements.jl` (`_tc_link_neighbors!`) | Fixed | Sorted `vert_line_x` once with `searchsortedfirst` lookup outside the candidate loop, and sorted trapezoid indices by `x_left` (`order`, `x_lefts`) to check only `x_left in [xr - tol, xr + tol]`. Verified in `test/test_arrangements.jl` (`< 1 MB` on 40-line arrangement; 160-line chain linking ~9.4 ms vs ~1.06 s). |
| 3 | Major | `src/Arrangements.jl` (`_vd_walls`) | Fixed | Added `_nearest_x_idx` binary search on sorted `unique_x` replacing linear `findfirst`, and replaced the four per-vertex `O(k)` temporary vectors with direct running `wall_top`/`wall_bot` tracking. Verified by `test/test_arrangements.jl`. |
| 4 | Major | `src/Bundles.jl` (`_unique_points`) | Fixed | Replaced `O(m^2)` linear scan with `sortperm` by `(p.x, p.y)`, backward `x`-window early break, and earliest-representative tracking (`rep`) preserving first-seen order; generalized element type parameter `T`. Verified in `test/test_bundles.jl` and `test/test_hyperbolas.jl`. |
| 5 | Moderate | `src/Arrangements.jl` (`_vd_slab_lines!`, `_vd_advance_slab!`, `vertical_decomposition`) | Fixed | Preallocated and reused `slab_lines`, `unique_slab`, and `bounds` across slabs; replaced `_vd_find_active` + `deleteat!` with an `O(k)` two-pointer merge in `_vd_advance_slab!`. Verified across `test/test_arrangements.jl` and 50 adversarial arrangements. |
| 6 | Moderate | `src/PntSeqs.jl` (`at`, `sample_uniformly`) | Fixed | Standardized `at` return branches on `Point{D,float(T)}`; rewrote `sample_uniformly` as a single-pass `O(m + n)` scan reusing `prefix_lengths(P)` once into a preallocated `Vector{Point{D,float(T)}}`. Verified via `@inferred` and `< 100 KB` allocation test for 1,000-point `sample_uniformly` (~24 KB vs ~8.1 MB) in `test/runtests.jl`. |
| 7 | Moderate | `src/PntSeqHierarchy.jl` (`SimplificationHierarchy`, `_fill!`, `_p_extract!`) | Fixed | Added `R<:SplitRule` type parameter to `SimplificationHierarchy{D,T,R,M}` so `plt::Palette{T,R}` is concrete; replaced `0.0` base-case returns in `_fill!` and `_p_extract!` with `zero(eltype(errs))` and `zero(eltype(pal))`. Verified with `isconcretetype` and `@inferred` tests on `Float32` in `test/test_pnt_seq_hierarchy.jl`. |
| 8 | Moderate | `src/PntSeqs.jl` (`dist(::AbsPntSeq, i, j)`) | Fixed | Changed `i == j` branch from `0.0` to `zero(float(T))`. Verified with `@inferred dist(p32, 1, 1) === 0.0f0` in `test/runtests.jl`. |
| 9 | Moderate | `src/LocalIncidence.jl` (`vertical_ray_segment`, `first_positive_vertical_contact`) | Fixed | Added compile-time `_check_exact_scalar` (in `src/ChainIndex.jl`) and early `x`-span / ray-side rejection in `vertical_ray_segment` before `_exact` (`BigInt`/`Rational{BigInt}`) promotion. Verified with `ArgumentError` on out-of-span inexact coordinates and `0` allocations on non-spanning segments in `test/test_chain_incidence.jl`. |

## Quality Gate & Build Verification

- **`tools/build`**: Passed (exit code 0, clean precompilation).
- **`tools/gate`** (`JULIA_PKG_OFFLINE=true`): Passed all **156,479 / 156,479** assertions (including Aqua checks and all newly added regression tests) in 1m26.4s.
- **`CHANGELOG.md`**: Recorded user-visible fixes and performance improvements under `## [Unreleased]`.
