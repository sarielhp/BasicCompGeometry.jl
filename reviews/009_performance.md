# Systems Code Review Report #009 (Lens: performance)

- **Date**: 2026-11-14
- **Auditor**: Gemini 3.8 Flash (Tier 0 Workhorse) via `review-audit`
- **Focus Lens**: `performance`
- **Model Tier**: `TIER0`
- **Backend**: `gemini`
- **Scope**: `src/Spheres.jl src/Ellipses.jl src/Hyperbolas.jl src/Parabolas.jl src/Bezier.jl src/CurvePolygon.jl`
- **Status**: Action Required

---

Auditing via Gemini Flash (bws run) [Profile: performance]...
[SEVERITY]: Moderate
[LOCATION]: [src/Bezier.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/Bezier.jl#L108-L121), [src/Bezier.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/Bezier.jl#L184-L191)
[ROOT CAUSE]: `bbox(b::CubicBezier{D,T})` allocates heap vectors `min_coords = zeros(T, D)`, `max_coords = zeros(T, D)`, `vals = T[v0, v3]` per coordinate, and a `Vector{T}` in `_coord_roots_in_unit_interval`, while evaluating all `D` coordinates of `b(t)` only to index `[k]`. Splatting `Point{D,T}(min_coords...)` and `Point{D,T}(max_coords...)` routes through `Core._apply_iterate` returning `Any`, breaking return-type inference (`BBox` instead of `BBox{D,T}`) for both `bbox(::CubicBezier)` and `bbox(::CubicSpline)` and forcing dynamic calls to `BBox` and `bound!`.
[FAILURE TRACE]:
- `@code_warntype bbox(CubicBezier(Point2F(0,0), Point2F(0.3,1.0), Point2F(0.7,-0.5), Point2F(1,0)))` infers `Core._apply_iterate(Base.iterate, SVector{2, Float64}, min_coords)::ANY` and `Body::BBox` (uninferred `UnionAll`), allocating `880` bytes per call.
- `@code_warntype bbox(sp)` on a 100-segment `CubicSpline{2,Float64}` infers `bb::BBox` throughout the loop and allocates `78,208` bytes.
[REMEDIATION]:
```julia
@inline function _eval_1d_bezier(v0::T, v1::T, v2::T, v3::T, t::T) where {T}
    u = one(T) - t
    return u * (u * (u * v0 + t * v1) + t * (u * v1 + t * v2)) +
           t * (u * (u * v1 + t * v2) + t * (u * v2 + t * v3))
end

function _coord_minmax(v0::T, v1::T, v2::T, v3::T) where {T}
    lo, hi = minmax(v0, v3)
    A = 3 * (-v0 + 3 * v1 - 3 * v2 + v3)
    B = 6 * (v0 - 2 * v1 + v2)
    C = 3 * (v1 - v0)
    if abs(A) < 1e-12
        if abs(B) > 1e-12
            t = -C / B
            if zero(T) < t < one(T)
                val = _eval_1d_bezier(v0, v1, v2, v3, t)
                lo = min(lo, val)
                hi = max(hi, val)
            end
        end
    else
        disc = B * B - 4 * A * C
        if disc >= zero(T)
            s = sqrt(disc)
            for t in ((-B - s) / (2 * A), (-B + s) / (2 * A))
                if zero(T) < t < one(T)
                    val = _eval_1d_bezier(v0, v1, v2, v3, t)
                    lo = min(lo, val)
                    hi = max(hi, val)
                end
            end
        end
    end
    return lo, hi
end

function bbox(b::CubicBezier{D,T}) where {D,T}
    extrema_tup = ntuple(k -> _coord_minmax(b.p0[k], b.p1[k], b.p2[k], b.p3[k]), Val(D))
    mini = Point{D,T}(ntuple(k -> extrema_tup[k][1], Val(D)))
    maxi = Point{D,T}(ntuple(k -> extrema_tup[k][2], Val(D)))
    return BBox(mini, maxi)
end
```

[SEVERITY]: Moderate
[LOCATION]: [src/Bezier.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/Bezier.jl#L267-L312)
[ROOT CAUSE]: In `interpolate_natural_spline`, lines 305–308 construct `p1` and `p2` inside `for i = 1:m` by allocating two `Vector{T}` comprehensions (`[P1_coords[k][i] for k = 1:D]` and `[2 * pts[i+1][k] - P1_coords[k][i+1] for k = 1:D]`) and splatting them with `...` into `Point{D,T}(...)`. This forces `Core._apply_iterate` returning `ANY`, making `p1::ANY`, `p2::ANY`, and dynamically dispatching `CubicBezier(p0, p1, p2, p3)` and `push!(segs, ...)` on every segment. In addition, the `for k = 1:D` Thomas loop allocates five `m`-element vectors (`K`, `a`, `b`, `c`, `r`) per coordinate even though the tridiagonal coefficients are identical across dimensions, and `segs` is grown via `push!` without preallocation.
[FAILURE TRACE]:
- Calling `interpolate_natural_spline(pts)` on a 101-point `PntSeq{2,Float64}` allocates `72,528` bytes.
- `@code_warntype` on the inner body confirms `p1::ANY`, `p2::ANY`, and unspecialized calls `CubicBezier(p0, p1, p2, p3)::CubicBezier{2}` and `push!(segs, ...)` on every segment iteration.
[REMEDIATION]:
```julia
function interpolate_natural_spline(pts::AbsPntSeq{D,T}; closed::Bool = false) where {D,T}
    n = cardin(pts)
    n >= 3 || return interpolate_catmull_rom(pts; closed = closed)
    m = n - 1
    P1_coords = ntuple(_ -> Vector{T}(undef, m), Val(D))
    b = Vector{T}(undef, m)
    r = Vector{T}(undef, m)

    for k = 1:D
        fill!(b, 4 * one(T))
        b[1] = 2 * one(T)
        b[m] = 7 * one(T)
        r[1] = pts[1][k] + 2 * pts[2][k]
        for i = 2:(m-1)
            r[i] = 4 * pts[i][k] + 2 * pts[i+1][k]
        end
        r[m] = 8 * pts[m][k] + pts[m+1][k]

        for i = 2:m
            a_i = i == m ? 2 * one(T) : one(T)
            w = a_i / b[i-1]
            b[i] -= w
            r[i] -= w * r[i-1]
        end
        p1k = P1_coords[k]
        p1k[m] = r[m] / b[m]
        for i = (m-1):-1:1
            p1k[i] = (r[i] - p1k[i+1]) / b[i]
        end
    end

    segs = Vector{CubicBezier{D,T}}(undef, m)
    for i = 1:m
        p0 = pts[i]
        p3 = pts[i+1]
        p1 = Point{D,T}(ntuple(k -> P1_coords[k][i], Val(D)))
        p2 = if i < m
            Point{D,T}(ntuple(k -> 2 * pts[i+1][k] - P1_coords[k][i+1], Val(D)))
        else
            Point{D,T}(ntuple(k -> (pts[n][k] + P1_coords[k][m]) / 2, Val(D)))
        end
        segs[i] = CubicBezier{D,T}(p0, p1, p2, p3)
    end
    return CubicSpline(segs; closed = false)
end
```

[SEVERITY]: Moderate
[LOCATION]: [src/Bezier.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/Bezier.jl#L193-L205)
[ROOT CAUSE]: `flatten(s::CubicSpline{D,T}; tol::Real = 1e-3)` calls `flatten(seg; tol = tol).pnts` for every segment (allocating a temporary `Vector{Point{D,T}}` and `PntSeq` wrapper per segment), and for every `i > 1` slices `pts[2:end]` without `@view`, allocating a second temporary `Vector{Point{D,T}}` before copying into `all_pts`.
[FAILURE TRACE]:
- Calling `flatten(sp)` on a 100-segment `CubicSpline{2,Float64}` allocates `150,056` bytes.
- `@code_warntype` on the inner body shows `Base.getindex(pts, 2:lastindex(pts))::Vector{SVector{2, Float64}}` allocating a slice copy on each segment `i >= 2`.
[REMEDIATION]:
```julia
function flatten(s::CubicSpline{D,T}; tol::Real = 1e-3, max_depth::Int = 10) where {D,T}
    isempty(s.segments) && return PntSeq{D,T}(Point{D,T}[])
    all_pts = Point{D,T}[s.segments[1].p0]
    sizehint!(all_pts, length(s.segments) + 1)
    tol_T = T(tol)
    for seg in s.segments
        _flatten_rec!(all_pts, seg, tol_T, max_depth)
    end
    return PntSeq(all_pts)
end
```

[SEVERITY]: Moderate
[LOCATION]: [src/CurvePolygon.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/CurvePolygon.jl#L208-L270), [src/CurvePolygon.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/CurvePolygon.jl#L278-L325), [src/CurvePolygon.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/CurvePolygon.jl#L354-L370)
[ROOT CAUSE]: `intersect_line_curve(line, seg)` allocates `res = Tuple{Point2F,Float64}[]` on every call (even on 0 intersections), and `intersect_line_curve(line, arc)` allocates both `res` and `candidate_pts = Point2F[]` on every call despite having at most two intersections. Both `is_inside(q, cp)` and `tangent_intersections_cp(CP, Q, alpha)` invoke `intersect_line_curve` in `for c in cp.curves`, and `tangent_intersections_cp` additionally allocates `intersections = Tuple{Point2F,Float64}[]` via `push!` only to immediately reduce it to `(u_pt, u_t)` and `(v_pt, v_t)`.
[FAILURE TRACE]:
- On a 100-curve `CurvePolygon2DF` (50 `Segment2F`, 50 `CircleArc2F`), a single warmed-up `is_inside(Point2F(0.5, 0.5), cp100)` query allocates `18,400` bytes and `tangent_intersections_cp(cp100, ps, 0.4)` allocates `4,240` bytes.
[REMEDIATION]:
```julia
@inline function _emit_line_curve_intersections(f, line::Line{2,Float64}, seg::Segment2F)
    dp = seg.q - seg.p
    det_m = line.u[1] * (-dp[2]) - line.u[2] * (-dp[1])
    abs(det_m) < 1e-12 && return nothing
    diff = seg.p - line.p
    t = (line.u[1] * diff[2] - line.u[2] * diff[1]) / det_m
    if 0.0 <= t <= 1.0
        f(at(seg, t), t)
    end
    return nothing
end

@inline function _emit_line_curve_intersections(f, line::Line{2,Float64}, arc::CircleArc2F)
    u_norm_sq = dot(line.u, line.u)
    u_norm_sq == 0.0 && return nothing
    t_proj = dot(arc.center - line.p, line.u) / u_norm_sq
    proj = line.p + t_proj * line.u
    d_sq = dist_sq(arc.center, proj)
    r_sq = arc.radius^2
    d_sq > r_sq + 1e-12 && return nothing

    h = sqrt(max(0.0, r_sq - d_sq))
    u_unit = line.u / sqrt(u_norm_sq)
    d_theta = arc.theta2 - arc.theta1
    if d_theta < 0
        d_theta = mod(d_theta, 2pi)
    end

    n_cand = h < 1e-12 ? 1 : 2
    for idx = 1:n_cand
        pt = n_cand == 1 ? proj : (idx == 1 ? proj + h * u_unit : proj - h * u_unit)
        angle = atan(pt[2] - arc.center[2], pt[1] - arc.center[1])
        rel_angle = mod(angle - arc.theta1, 2pi)
        if rel_angle <= d_theta + 1e-9 || (2pi - rel_angle) <= 1e-9
            local_t = clamp(rel_angle / d_theta, 0.0, 1.0)
            f(at(arc, local_t), local_t)
        end
    end
    return nothing
end
```
Use `_emit_line_curve_intersections` directly inside `is_inside` and `tangent_intersections_cp` (updating `crossings` or `min_t_pos`/`max_t_neg` in place without allocating `res`, `candidate_pts`, or `intersections`).

[SEVERITY]: Moderate
[LOCATION]: [src/CurvePolygon.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/CurvePolygon.jl#L127-L158)
[ROOT CAUSE]: When `closed == true`, `point_on(pnt_seq::AbsPntSeq{D,T}, t::Real; closed::Bool = false)` returns `Point{D,T}` on boundary branches (`n == 1`, `target <= 0.0`, `target >= total_len`, `i <= 1`) while the interpolation branch returns `convex_comb(..., local_t::Float64)` of type `Point{D,Float64}`, creating a value-dependent return type whenever `T !== Float64`. Additionally, `vcat(pnts, [pnts[1]])` allocates two temporary vectors (`[pnts[1]]` and the concatenated copy) alongside `lens = zeros(m)` on every `closed == true` call.
[FAILURE TRACE]:
- `@code_warntype BasicCompGeometry.point_on(PntSeq([Point2I(0,0), Point2I(1,0), Point2I(1,1), Point2I(0,1)]), 0.5; closed=true)` infers `Body::Union{SVector{2, Float64}, SVector{2, Int64}}`.
- Calling `BasicCompGeometry.point_on(PntSeq([Point2F(0,0), Point2F(1,0), Point2F(1,1), Point2F(0,1)]), 0.5; closed=true)` allocates `352` bytes for a 4-point sequence.
[REMEDIATION]:
```julia
function point_on(pnt_seq::AbsPntSeq{D,T}, t::Real; closed::Bool = false) where {D,T}
    n = length(pnt_seq)
    n == 0 && throw(ArgumentError("Empty point sequence"))
    R = float(T)
    n == 1 && return Point{D,R}(Points(pnt_seq)[1])

    !closed && return Point{D,R}(at(pnt_seq, t))

    pnts = Points(pnt_seq)
    need_close = pnts[1] != pnts[end]
    n_edges = need_close ? n : n - 1
    total_len = 0.0
    for i = 1:n_edges
        p_next = (need_close && i == n) ? pnts[1] : pnts[i+1]
        total_len += dist(pnts[i], p_next)
    end

    target = clamp(Float64(t), 0.0, 1.0) * total_len
    target <= 0.0 && return Point{D,R}(pnts[1])
    target >= total_len && return Point{D,R}(need_close ? pnts[1] : pnts[end])

    acc = 0.0
    for i = 1:n_edges
        p_next = (need_close && i == n) ? pnts[1] : pnts[i+1]
        seg_len = dist(pnts[i], p_next)
        if acc + seg_len >= target
            local_t = seg_len == 0.0 ? 0.0 : (target - acc) / seg_len
            return Point{D,R}(convex_comb(pnts[i], p_next, local_t))
        end
        acc += seg_len
    end
    return Point{D,R}(need_close ? pnts[1] : pnts[end])
end
```

[SEVERITY]: Moderate
[LOCATION]: [src/Ellipses.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/Ellipses.jl#L226-L254)
[ROOT CAUSE]: `bbox(arc::EllipticArc{T})` allocates `pts = Point{2,T}[start_point(arc), end_point(arc)]`, conditionally `push!`es up to 4 critical points, and then runs four separate `minimum`/`maximum` generator passes over `pts` instead of tracking scalar extrema directly.
[FAILURE TRACE]:
- `@allocated bbox(EllipticArc(Ellipse(Point2F(0.0, 0.0), 2.0, 1.0, 0.3), 0.0, 1.5))` reports `352` bytes per call.
[REMEDIATION]:
```julia
function bbox(arc::EllipticArc{T}) where {T}
    p1 = start_point(arc)
    p2 = end_point(arc)
    min_x, max_x = minmax(p1[1], p2[1])
    min_y, max_y = minmax(p1[2], p2[2])

    e = arc.ellipse
    ca, sa = cos(e.angle), sin(e.angle)

    denom_x = e.r_major * ca
    if abs(denom_x) > 1e-12
        t_x = atan(-e.r_minor * sa, denom_x)
        for t in (t_x, t_x + π)
            if _angle_in_ccw_span(t, arc.alpha1, arc.alpha2)
                p = point_on(e, t)
                min_x = min(min_x, p[1]); max_x = max(max_x, p[1])
                min_y = min(min_y, p[2]); max_y = max(max_y, p[2])
            end
        end
    end

    denom_y = e.r_major * sa
    if abs(denom_y) > 1e-12
        t_y = atan(e.r_minor * ca, denom_y)
        for t in (t_y, t_y + π)
            if _angle_in_ccw_span(t, arc.alpha1, arc.alpha2)
                p = point_on(e, t)
                min_x = min(min_x, p[1]); max_x = max(max_x, p[1])
                min_y = min(min_y, p[2]); max_y = max(max_y, p[2])
            end
        end
    end

    return BBox(Point{2,T}(min_x, min_y), Point{2,T}(max_x, max_y))
end
```
