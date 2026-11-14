# Post-Work Summary: Performance Remediation Cycle #009

- **Source Report**: `reviews/009_performance.md`
- **Plan**: `reviews/009_performance_plan.md`
- **Date**: 2026-11-14

## Issue-by-Issue Status & Mitigation

| # | Severity | Location | Status | Mitigation & Verification |
|---|---|---|---|---|
| 1 | Moderate | `src/Bezier.jl` (`bbox(::CubicBezier)`, `bbox(::CubicSpline)`) | Fixed | Replaced `zeros(T, D)`, `_coord_roots_in_unit_interval` (`Vector{T}`), full `D`-D `b(t)[k]` calls, and coordinate splatting with stack-allocated `_eval_1d_bezier`, `_update_coord_extrema`, `_coord_minmax`, and `_bezier_extrema` via `ntuple(..., Val(D))`; updated `bbox(::CubicSpline)` to accumulate segment extrema in-place with `bound!`. Verified `@inferred bbox(b2) isa BBox{2,Float64}` (`96` bytes vs `880` bytes), `@inferred bbox(b3) isa BBox{3,Float64}`, `@inferred bbox(b32) isa BBox{2,Float32}`, and `@inferred bbox(sp) isa BBox{2,Float64}` (`96` bytes vs `73,600` bytes) in `test/test_bezier.jl`. |
| 2 | Moderate | `src/Bezier.jl` (`interpolate_natural_spline`) | Fixed | Extracted `_solve_natural_tridiag!` reusing two `m`-element work buffers (`b`, `r`) across dimensions into `P1_coords = ntuple(_ -> Vector{T}(undef, m), Val(D))`, preallocated `segs = Vector{CubicBezier{D,T}}(undef, m)`, and constructed `p1`/`p2` with `ntuple(..., Val(D))` instead of `[...]...` comprehension splatting. Verified `@inferred interpolate_natural_spline(knots)` (`CubicSpline{2,Float64}`) and allocation drop from `72,528` bytes to `10,184` bytes on 101 knots (`< 16,000` test bound) in `test/test_bezier.jl`. |
| 3 | Moderate | `src/Bezier.jl` (`flatten(::CubicSpline)`) | Fixed | Preallocated `all_pts = Point{D,T}[s.segments[1].p0]` with `sizehint!` and called `_flatten_rec!(all_pts, seg, tol_T, max_depth)` directly per segment (adding `max_depth::Int = 10` keyword), eliminating per-segment `PntSeq` allocations and `pts[2:end]` slices. Verified `@inferred flatten(sp)`, exact point equivalence against per-segment flattening, and allocation drop from `28,600` bytes to `1,888` bytes (`< 4,000` test bound) in `test/test_bezier.jl`. |
| 4 | Moderate | `src/CurvePolygon.jl` (`intersect_line_curve`, `tangent_intersections_cp`, `is_inside`) | Fixed | Added stack-allocated `isbits` `_LineCurveHits` with `_line_curve_hits`, `_circle_arc_hit`, `_combine_hits`, `_dispatch_line_curve_hits`, and `_curve_length`; kept `intersect_line_curve` APIs via `_hits_to_vector`; updated `is_inside` and `tangent_intersections_cp` to accumulate running arc length `acc_len` and consume hits directly without `Vector` or `prefix_lengths` allocations. Verified `@inferred` and `0`-byte allocations (down from `18,400` and `4,240` bytes) on a 100-curve `Segment2F`/`CircleArc2F` `CurvePolygon2D` in `test/test_curve_polygon.jl`. |
| 5 | Moderate | `src/CurvePolygon.jl` (`point_on(::AbsPntSeq; closed=...)`) | Fixed | Standardized all return branches on `Point{D,float(T)}` and extracted `_point_on_closed` using two-pass index wraparound and zero-length edge guarding, eliminating `vcat(pnts, [pnts[1]])` and `zeros(m)`. Verified `@inferred` return types across `Float64`, `Float32`, and `Int64` plus `0`-byte allocations (down from `352` bytes) in `test/test_curve_polygon.jl`. |
| 6 | Moderate | `src/Ellipses.jl` (`bbox(::EllipticArc)`) | Fixed | Extracted `_update_arc_extrema` and tracked scalar `min_x`, `max_x`, `min_y`, `max_y` directly from `start_point`, `end_point`, and critical angles instead of allocating `pts = Point{2,T}[...]` and four generator passes. Verified `@inferred bbox(arc)` (`96` bytes vs `352` bytes) for `Float64` and `Float32` in `test/test_ellipses.jl`. |

## Quality Gate & Build Verification

- **`./tools/build`**: Passed (exit code 0, clean package precompilation).
- **`./tools/gate`**: Passed (exit code 0), running Aqua checks (`Method ambiguity`, `Unbound type parameters`, `Undefined exports`, `Compare Project.toml and test/Project.toml`, `Stale dependencies`, `Compat bounds`, `Piracy`, `Persistent tasks`) and passing **156,513 / 156,513** assertions in `1m26.2s`.
- **`CHANGELOG.md`**: Recorded fixes and performance improvements under `## [Unreleased]`.
- **Status**: Remediated & Verified
- **Footprint**: 9 files changed, 499 insertions(+), 197 deletions(-)
- **Differential Audit**: Clean (0 defects in diff)
