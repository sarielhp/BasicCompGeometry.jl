# Architectural Plan: Performance Remediation Cycle #009

## 1. Triage of Findings (`reviews/009_performance.md`)

All 6 findings in `reviews/009_performance.md` were reproduced against the current codebase with `@allocated` and `@code_warntype` and confirmed as **genuine defects**. None are false positives or intentional design choices.

| # | Severity | Location | Summary of Defect | Disposition |
|---|---|---|---|---|
| 1 | Moderate | `src/Bezier.jl` (`bbox(::CubicBezier)`, `bbox(::CubicSpline)`) | Per-coordinate `Vector{T}` / `zeros(T, D)` allocations, full `D`-D evaluation in `b(t)[k]`, and `...` splatting (`Core._apply_iterate` -> `ANY`) breaking `BBox{D,T}` return-type inference for both `CubicBezier` and `CubicSpline` | Fix |
| 2 | Moderate | `src/Bezier.jl` (`interpolate_natural_spline`) | Per-dimension 5-vector allocations in Thomas solver, `[... ]...` comprehension splatting for `p1`/`p2` (`ANY` + dynamic dispatch on `CubicBezier`/`push!`), and unpreallocated `segs` | Fix |
| 3 | Moderate | `src/Bezier.jl` (`flatten(::CubicSpline)`) | Per-segment `flatten(seg)` `PntSeq` allocation plus non-view `pts[2:end]` slice copy on every segment `i >= 2` | Fix |
| 4 | Moderate | `src/CurvePolygon.jl` (`intersect_line_curve`, `tangent_intersections_cp`, `is_inside`) | Per-curve `Vector` allocations in `intersect_line_curve` (`res`, `candidate_pts`), `prefix_lengths` + `intersections` allocations in `tangent_intersections_cp`, and 5-way `Curve2D` union dynamic dispatch | Fix |
| 5 | Moderate | `src/CurvePolygon.jl` (`point_on(::AbsPntSeq; closed=...)`) | Value-dependent `Union{Point{D,T}, Point{D,Float64}}` return type for `T !== Float64`, plus `vcat(pnts, [pnts[1]])` and `zeros(m)` allocations when `closed=true` | Fix |
| 6 | Moderate | `src/Ellipses.jl` (`bbox(::EllipticArc)`) | Allocates `pts = Point{2,T}[...]` and runs four generator `minimum`/`maximum` passes instead of scalar extrema tracking | Fix |

---

## 2. Targeted Remediation Plan

### 2.1 `src/Bezier.jl` — Zero-Allocation Extrema, Natural Spline, and Spline Flattening (Findings 1, 2, 3)
1. **Type-stable, allocation-free Bézier `bbox` (Finding 1)**:
   - Add `@inline _eval_1d_bezier(v0, v1, v2, v3, t)`, `@inline _update_coord_extrema`, `_coord_minmax`, and `_bezier_extrema(b::CubicBezier{D,T})` using `ntuple(..., Val(D))` to compute `(mini::Point{D,T}, maxi::Point{D,T})` on the stack with zero heap allocations.
   - Implement `bbox(b::CubicBezier{D,T}) = BBox(_bezier_extrema(b)...)` and update `bbox(s::CubicSpline{D,T})` to accumulate segment extrema in-place onto `bb` via `bound!(bb, mini); bound!(bb, maxi)`, eliminating per-segment `BBox` allocations and restoring concrete `BBox{D,T}` inference.
   - Remove the now-unused `_coord_roots_in_unit_interval` helper.
2. **Concrete, zero-splat `interpolate_natural_spline` (Finding 2)**:
   - Preallocate `P1_coords = ntuple(_ -> Vector{T}(undef, m), Val(D))` and reuse two `m`-element scratch buffers `b` and `r` across dimensions via helper `_solve_natural_tridiag!`.
   - Preallocate `segs = Vector{CubicBezier{D,T}}(undef, m)` and construct `p1` and `p2` via `ntuple(..., Val(D))` directly into `CubicBezier{D,T}(p0, p1, p2, p3)`.
3. **In-place `flatten(::CubicSpline)` (Finding 3)**:
   - Preallocate `all_pts = Point{D,T}[s.segments[1].p0]` with `sizehint!(all_pts, length(s.segments) + 1)` and append recursively via `_flatten_rec!(all_pts, seg, tol_T, max_depth)` (supporting `max_depth::Int = 10` keyword to match `flatten(::CubicBezier)`).

### 2.2 `src/CurvePolygon.jl` — Zero-Allocation Curve Hits and Type-Stable `point_on` (Findings 4, 5)
1. **Stack-allocated `_LineCurveHits` for `Segment2F` and `CircleArc2F` (Finding 4)**:
   - Define `isbits` struct `_LineCurveHits` (up to 2 intersections) and zero-allocation `_line_curve_hits(line, seg::Segment2F)` and `_line_curve_hits(line, arc::CircleArc2F)` with `_circle_arc_hit` and `_combine_hits` helpers.
   - Keep `intersect_line_curve` public methods returning `Vector{Tuple{Point2F, Float64}}` by wrapping `_line_curve_hits` via `_hits_to_vector`.
   - Add `_dispatch_line_curve_hits` and `_curve_length` with explicit `Curve2D` subtype branches (`Segment2F`, `CircleArc2F`, `ParabolaArc`, `Parabola`, `Hyperbola`) so Julia avoids 5-way union-splitting fallbacks to dynamic dispatch while preserving full `Curve2D` support (`ParabolaArc`, etc.).
   - Rewrite `is_inside` and `tangent_intersections_cp` to iterate `_dispatch_line_curve_hits` directly without allocating `res`, `candidate_pts`, `intersections`, or `prefix_lengths(CP)` (tracking running `acc_len` and normalizing once at the end).
2. **Type-stable, zero-allocation `point_on(::AbsPntSeq; closed=...)` (Finding 5)**:
   - Standardize all return paths to `Point{D,R}` where `R = float(T)` (matching `at(pnt_seq, t)`), resolving `Union{Point{D,T}, Point{D,Float64}}` instability on integer sequences while preserving `Float32` on `Float32` sequences.
   - Extract `_point_on_closed(pnts, t, ::Type{Point{D,R}})` using index wraparound over `n_edges` and a zero-length edge guard (`seg_len == zero(R)`), eliminating `vcat(pnts, [pnts[1]])` and `zeros(m)`.

### 2.3 `src/Ellipses.jl` — Scalar Extrema in `bbox(::EllipticArc)` (Finding 6)
- Extract `_update_arc_extrema` to avoid duplicating the critical-angle check between `t_x` and `t_y`.
- Track scalar `min_x`, `max_x`, `min_y`, `max_y` directly in `bbox(arc::EllipticArc{T})`, eliminating the `pts` vector and four generator passes.

---

## 3. Verification & Quality Standards

- **Regression tests**:
  - `test/test_bezier.jl`: add `@inferred` and `@allocated` tests for `bbox(::CubicBezier)`, `bbox(::CubicSpline)`, `flatten(::CubicSpline)`, and `interpolate_natural_spline` (including 2D/3D and `Float32`/`Float64`).
  - `test/test_curve_polygon.jl`: add `@inferred` and `@allocated` tests for `is_inside`, `tangent_intersections_cp` (on `Segment2F`/`CircleArc2F` `CurvePolygon2D`), and `point_on(...; closed=true)` across `Float64`, `Float32`, and `Int64` (including zero-length edge and already-closed `PntSeq` cases).
  - `test/test_ellipses.jl`: add `@inferred` and `@allocated` checks for `bbox(::EllipticArc)` on both `Float64` and `Float32`.
- **Quality gate**: Run `./tools/build` and `./tools/gate` after edits; verify 100% pass rate, zero regressions, shallow nesting (`<= 4`), concise helper sizing, and preserved public APIs.
- **Changelog**: Record user-visible fixes and performance improvements in `CHANGELOG.md` under `## [Unreleased]`.
