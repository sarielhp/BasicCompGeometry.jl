# Curve type and CurvePolygon2D

"""
    Curve2D{T}

Type alias for a 2D curve primitive, which can be a `Segment{2, T}`, `CircleArc{T}`, or `Hyperbola{T}`.
"""
const Curve2D{T} = Union{Segment{2,T},CircleArc{T},Hyperbola{T},Parabola{T},ParabolaArc{T}}
const Curve2DF = Curve2D{Float64}

"""
    CurvePolygon2D{T}

A 2D polygon defined by a closed sequence of curves (`Curve2D{T}`) connecting consecutive vertices.
"""
struct CurvePolygon2D{T}
    curves::Vector{Curve2D{T}}
end

function CurvePolygon2D(curves::AbstractVector)
    curves_vec = Vector{Curve2DF}(curves)
    return CurvePolygon2D{Float64}(curves_vec)
end

const CurvePolygon2DF = CurvePolygon2D{Float64}

# eltype for curves
Base.eltype(::Segment{2,T}) where {T} = T
Base.eltype(::CircleArc{T}) where {T} = T

"""
    geom_length(arc::CircleArc)

Calculate the arc length of a `CircleArc`.
The angular span is `mod(theta2 - theta1, 2pi)` (or `theta2 - theta1` if `theta2 >= theta1`).
"""
function geom_length(arc::CircleArc{T}) where {T}
    d_theta = arc.theta2 - arc.theta1
    if d_theta < 0
        d_theta = mod(d_theta, T(2pi))
    end
    return arc.radius * d_theta
end

"""
    geom_length(cp::CurvePolygon2D)

Compute the total length of the curves forming a `CurvePolygon2D`.
"""
function geom_length(cp::CurvePolygon2D{T}) where {T}
    len = 0.0
    for c in cp.curves
        len += geom_length(c)
    end
    return len
end

"""
    prefix_lengths(cp::CurvePolygon2D)

Return a vector containing the cumulative lengths after each curve of `cp`, starting with 0.0.
"""
function prefix_lengths(cp::CurvePolygon2D{T}) where {T}
    n = length(cp.curves)
    v = Vector{Float64}(undef, n + 1)
    v[1] = 0.0
    for i in 1:n
        v[i+1] = v[i] + geom_length(cp.curves[i])
    end
    return v
end

"""
    at(arc::CircleArc, t::Real)

Return the point on `arc` at relative parameter `t` in `[0, 1]`.
`t=0` corresponds to `theta1` and `t=1` corresponds to `theta2`.
"""
function at(arc::CircleArc{T}, t::Real) where {T}
    d_theta = arc.theta2 - arc.theta1
    if d_theta < 0
        d_theta = mod(d_theta, T(2pi))
    end
    angle = arc.theta1 + T(t) * d_theta
    return arc.center + Point{2,T}(arc.radius * cos(angle), arc.radius * sin(angle))
end

"""
    point_on(cp::CurvePolygon2D, t::Real)

Return the point on `cp` at relative length `t` in `[0, 1]`, traversing counterclockwise.
"""
function point_on(cp::CurvePolygon2D{T}, t::Real) where {T}
    n = length(cp.curves)
    n == 0 && throw(ArgumentError("Empty CurvePolygon2D"))

    plens = prefix_lengths(cp)
    total_len = plens[end]
    total_len == 0.0 && return at(cp.curves[1], 0.0)

    target = clamp(Float64(t), 0.0, 1.0) * total_len
    if target <= 0.0
        return at(cp.curves[1], 0.0)
    elseif target >= total_len
        return at(cp.curves[end], 1.0)
    end

    i = searchsortedfirst(plens, target) - 1
    i = clamp(i, 1, n)

    curve_len = plens[i+1] - plens[i]
    if curve_len == 0.0
        local_t = 0.0
    else
        local_t = (target - plens[i]) / curve_len
    end

    return at(cp.curves[i], local_t)
end

function _point_on_closed(
    pnts::AbstractVector{Point{D,T}},
    t::Real,
    ::Type{Point{D,R}},
) where {D,T,R}
    n = length(pnts)
    need_close = pnts[1] != pnts[end]
    n_edges = need_close ? n : n - 1
    total_len = zero(R)
    for i in 1:n_edges
        p_next = (need_close && i == n) ? pnts[1] : pnts[i+1]
        total_len += R(dist(pnts[i], p_next))
    end

    target = clamp(R(t), zero(R), one(R)) * total_len
    target <= zero(R) && return Point{D,R}(pnts[1])
    end_pt = Point{D,R}(need_close ? pnts[1] : pnts[end])
    target >= total_len && return end_pt

    acc = zero(R)
    for i in 1:n_edges
        p_next = (need_close && i == n) ? pnts[1] : pnts[i+1]
        seg_len = R(dist(pnts[i], p_next))
        if acc + seg_len >= target
            local_t = seg_len == zero(R) ? zero(R) : (target - acc) / seg_len
            return Point{D,R}(convex_comb(Point{D,R}(pnts[i]), Point{D,R}(p_next), local_t))
        end
        acc += seg_len
    end
    return end_pt
end

"""
    point_on(pnt_seq::AbsPntSeq, t::Real)

Return the point on `pnt_seq` at relative distance `t` in `[0, 1]`, traversing along the sequence (or closed polygon).
If treated as a closed polygon, pass `point_on(pnt_seq, t; closed=true)`.
"""
function point_on(pnt_seq::AbsPntSeq{D,T}, t::Real; closed::Bool = false) where {D,T}
    n = length(pnt_seq)
    n == 0 && throw(ArgumentError("Empty point sequence"))
    R = float(T)
    n == 1 && return Point{D,R}(pnt_seq[1])
    !closed && return Point{D,R}(at(pnt_seq, t))
    return _point_on_closed(Points(pnt_seq), t, Point{D,R})
end

"""
    direction(alpha::Real)

Return a 2D unit vector `Point{2, Float64}` making angle `alpha` with the positive x-axis.
"""
@inline function direction(alpha::Real)
    return Point{2,Float64}(cos(Float64(alpha)), sin(Float64(alpha)))
end

"""
    tangent_line(poly::AbsPntSeq{2, T}, alpha::Real)

Compute the tangent line to a 2D convex polygon `poly` that has outer normal `v = direction(alpha)`.
Returns `(line, q)`, where `line` is a `Line{2, Float64}` (with point `q` and direction vector orthogonal to `v`) 
and `q` is the extreme point in `poly` in the direction of `v`.
"""
function tangent_line(poly::AbsPntSeq{2,T}, alpha::Real) where {T}
    v = direction(alpha)
    pnts = Points(poly)

    # Find extreme point q maximizing dot(p, v)
    max_dot = -Inf
    best_idx = 1
    for (i, p) in enumerate(pnts)
        d = dot(p, v)
        if d > max_dot
            max_dot = d
            best_idx = i
        end
    end
    q = Point{2,Float64}(pnts[best_idx])

    # Line direction orthogonal to v (counterclockwise rotation by 90 deg: (-v[2], v[1]))
    u = Point{2,Float64}(-v[2], v[1])
    line = Line(q, u)
    return line, q
end

function tangent_line(pnts::AbstractVector{Point{2,T}}, alpha::Real) where {T}
    return tangent_line(PntSeq(pnts), alpha)
end

struct _LineCurveHits
    n::Int
    p1::Point2F
    t1::Float64
    p2::Point2F
    t2::Float64
end

const _NO_HITS = _LineCurveHits(0, Point2F(0.0, 0.0), 0.0, Point2F(0.0, 0.0), 0.0)

@inline _one_hit(p1::Point2F, t1::Float64) =
    _LineCurveHits(1, p1, t1, Point2F(0.0, 0.0), 0.0)

@inline function _hits_to_vector(hits::_LineCurveHits)
    res = Tuple{Point2F,Float64}[]
    hits.n >= 1 && push!(res, (hits.p1, hits.t1))
    hits.n >= 2 && push!(res, (hits.p2, hits.t2))
    return res
end

@inline function _line_curve_hits(line::Line{2,Float64}, seg::Segment2F)
    dp = seg.q - seg.p
    det_m = line.u.x * (-dp.y) - line.u.y * (-dp.x)
    abs(det_m) < 1e-12 && return _NO_HITS
    diff = seg.p - line.p
    t = (line.u.x * diff.y - line.u.y * diff.x) / det_m
    if 0.0 <= t <= 1.0
        return _one_hit(at(seg, t), t)
    end
    return _NO_HITS
end

@inline function _circle_arc_hit(arc::CircleArc2F, pt::Point2F, d_theta::Float64)
    angle = atan(pt.y - arc.center.y, pt.x - arc.center.x)
    rel_angle = mod(angle - arc.theta1, 2pi)
    if rel_angle <= d_theta + 1e-9 || (2pi - rel_angle) <= 1e-9
        local_t = d_theta == 0.0 ? 0.0 : clamp(rel_angle / d_theta, 0.0, 1.0)
        return true, at(arc, local_t), local_t
    end
    return false, Point2F(0.0, 0.0), 0.0
end

@inline function _combine_hits(
    h1::Tuple{Bool,Point2F,Float64},
    h2::Tuple{Bool,Point2F,Float64},
)
    ok1, p1, t1 = h1
    ok2, p2, t2 = h2
    if ok1 && ok2
        return _LineCurveHits(2, p1, t1, p2, t2)
    elseif ok1
        return _one_hit(p1, t1)
    elseif ok2
        return _one_hit(p2, t2)
    end
    return _NO_HITS
end

@inline function _line_curve_hits(line::Line{2,Float64}, arc::CircleArc2F)
    u_norm_sq = dot(line.u, line.u)
    u_norm_sq == 0.0 && return _NO_HITS
    t_proj = dot(arc.center - line.p, line.u) / u_norm_sq
    proj = line.p + t_proj * line.u
    d_sq = dist_sq(arc.center, proj)
    r_sq = arc.radius^2
    d_sq > r_sq + 1e-12 && return _NO_HITS

    h = sqrt(max(0.0, r_sq - d_sq))
    u_unit = line.u / sqrt(u_norm_sq)
    d_theta = arc.theta2 - arc.theta1
    if d_theta < 0
        d_theta = mod(d_theta, 2pi)
    end

    if h < 1e-12
        ok, p1, t1 = _circle_arc_hit(arc, proj, d_theta)
        return ok ? _one_hit(p1, t1) : _NO_HITS
    end
    step = h * u_unit
    h1 = _circle_arc_hit(arc, proj + step, d_theta)
    h2 = _circle_arc_hit(arc, proj - step, d_theta)
    return _combine_hits(h1, h2)
end

@inline function _line_curve_hits_fallback(line::Line{2,Float64}, curve)
    pts_t = intersect_line_curve(line, curve)
    n = length(pts_t)
    n == 0 && return _NO_HITS
    p1, t1 = pts_t[1]
    n == 1 && return _one_hit(Point2F(p1), Float64(t1))
    p2, t2 = pts_t[2]
    return _LineCurveHits(2, Point2F(p1), Float64(t1), Point2F(p2), Float64(t2))
end

# Explicit branches avoid 5-way Union{Curve2D} dynamic dispatch inside tight loops.
@inline function _dispatch_line_curve_hits(line::Line{2,Float64}, c)
    if c isa Segment2F
        return _line_curve_hits(line, c)
    elseif c isa CircleArc2F
        return _line_curve_hits(line, c)
    elseif c isa ParabolaArc
        return _line_curve_hits_fallback(line, c)
    elseif c isa Parabola
        return _line_curve_hits_fallback(line, c)
    else
        return _line_curve_hits_fallback(line, c)
    end
end

@inline function _curve_length(c)
    if c isa Segment2F
        return geom_length(c)
    elseif c isa CircleArc2F
        return geom_length(c)
    elseif c isa ParabolaArc
        return geom_length(c)
    elseif c isa Parabola
        return geom_length(c)
    else
        return geom_length(c)
    end
end

"""
    intersect_line_curve(line::Line{2, Float64}, curve::Curve2D)

Compute all intersection points between infinite 2D line and `curve` (Segment2F or CircleArc2F).
Returns a vector of `(point::Point2F, t_curve::Float64)` where `t_curve ∈ [0, 1]` is the local curve parameter.
"""
intersect_line_curve(line::Line{2,Float64}, seg::Segment2F) =
    _hits_to_vector(_line_curve_hits(line, seg))

intersect_line_curve(line::Line{2,Float64}, arc::CircleArc2F) =
    _hits_to_vector(_line_curve_hits(line, arc))

"""
    tangent_intersections_cp(CP::CurvePolygon2D, Q, alpha::Real)

Compute the tangent line to `Q` at angle `alpha`, tangency point `q`, and the two intersection points `u` (right turn from O->q->u) and `v` (left turn from O->q->v) on `CP`.
Returns `(q, u, t_u, v, t_v)`.
"""
function tangent_intersections_cp(CP::CurvePolygon2D, Q, alpha::Real)
    line, q = tangent_line(Q, alpha)
    u_norm_sq = dot(line.u, line.u)

    # Along tangent line L(t) = q + t*line.u:
    # - t > 0 -> left turn (v, after_tangent_to_polygon)
    # - t < 0 -> right turn (u, before_tangent_to_polygon)
    u_pt, u_pos = q, 0.0
    v_pt, v_pos = q, 0.0
    min_t_pos = Inf
    max_t_neg = -Inf

    acc_len = 0.0
    for curve in CP.curves
        curve_len = Float64(_curve_length(curve))
        hits = _dispatch_line_curve_hits(line, curve)
        for k in 1:hits.n
            pt, local_t = k == 1 ? (hits.p1, hits.t1) : (hits.p2, hits.t2)
            t_line = dot(pt - q, line.u) / u_norm_sq
            if t_line < -1e-9 && t_line > max_t_neg
                max_t_neg = t_line
                u_pt = pt
                u_pos = acc_len + local_t * curve_len
            elseif t_line > 1e-9 && t_line < min_t_pos
                min_t_pos = t_line
                v_pt = pt
                v_pos = acc_len + local_t * curve_len
            end
        end
        acc_len += curve_len
    end

    u_t = isfinite(max_t_neg) ? u_pos / acc_len : 0.0
    v_t = isfinite(min_t_pos) ? v_pos / acc_len : 0.0
    return q, u_pt, u_t, v_pt, v_t
end

"""
    before_tangent_to_polygon(CP::CurvePolygon2D, Q, alpha::Real)

Return the parametric location t in [0, 1] of intersection point u on `CP`, 
where u is the intersection point on the tangent line to `Q` at angle `alpha` forming a right turn O -> q -> u.
"""
function before_tangent_to_polygon(CP::CurvePolygon2D, Q, alpha::Real)
    _, _, u_t, _, _ = tangent_intersections_cp(CP, Q, alpha)
    return u_t
end

"""
    after_tangent_to_polygon(CP::CurvePolygon2D, Q, alpha::Real)

Return the parametric location t in [0, 1] of intersection point v on `CP`, 
where v is the intersection point on the tangent line to `Q` at angle `alpha` forming a left turn O -> q -> v.
"""
function after_tangent_to_polygon(CP::CurvePolygon2D, Q, alpha::Real)
    _, _, _, _, v_t = tangent_intersections_cp(CP, Q, alpha)
    return v_t
end

"""
    is_inside(q::Point{2,S}, cp::CurvePolygon2D{T})

Test whether query point `q` lies strictly inside the 2D curved polygon `cp` by shooting a ray downward.
"""
function is_inside(q::Point{2,S}, cp::CurvePolygon2D{T}) where {S,T}
    ray = Line(Point{2,Float64}(q), Point{2,Float64}(0.0, -1.0))
    crossings = 0
    qy = Float64(q.y)
    for c in cp.curves
        hits = _dispatch_line_curve_hits(ray, c)
        for k in 1:hits.n
            pt, s = k == 1 ? (hits.p1, hits.t1) : (hits.p2, hits.t2)
            if pt.y < qy - 1e-11
                if s < 1.0 - 1e-9
                    crossings += 1
                end
            elseif abs(pt.y - qy) <= 1e-11
                return true # Point is on the boundary
            end
        end
    end
    return isodd(crossings)
end
