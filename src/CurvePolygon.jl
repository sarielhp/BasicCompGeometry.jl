# Curve type and CurvePolygon2D

"""
    Curve2D{T}

Type alias for a 2D curve primitive, which can be a `Segment{2, T}`, `CircleArc{T}`, or `Hyperbola{T}`.
"""
const Curve2D{T} = Union{Segment{2,T}, CircleArc{T}, Hyperbola{T}, Parabola{T}, ParabolaArc{T}}
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
    for i = 1:n
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

"""
    point_on(pnt_seq::AbsPntSeq, t::Real)

Return the point on `pnt_seq` at relative distance `t` in `[0, 1]`, traversing along the sequence (or closed polygon).
If treated as a closed polygon, pass `point_on(pnt_seq, t; closed=true)`.
"""
function point_on(pnt_seq::AbsPntSeq{D,T}, t::Real; closed::Bool=false) where {D,T}
    n = length(pnt_seq)
    n == 0 && throw(ArgumentError("Empty point sequence"))
    n == 1 && return pnt_seq.pnts[1]

    if !closed
        return at(pnt_seq, t)
    else
        pnts = Points(pnt_seq)
        closed_pnts = pnts[1] == pnts[end] ? pnts : vcat(pnts, [pnts[1]])
        m = length(closed_pnts)
        lens = zeros(m)
        for i = 2:m
            lens[i] = lens[i-1] + dist(closed_pnts[i-1], closed_pnts[i])
        end
        total_len = lens[end]
        target = clamp(Float64(t), 0.0, 1.0) * total_len

        if target <= 0.0
            return closed_pnts[1]
        elseif target >= total_len
            return closed_pnts[end]
        end

        i = searchsortedfirst(lens, target)
        i <= 1 && return closed_pnts[1]

        seg_len = lens[i] - lens[i-1]
        local_t = (target - lens[i-1]) / seg_len
        return convex_comb(closed_pnts[i-1], closed_pnts[i], local_t)
    end
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

"""
    intersect_line_curve(line::Line{2, Float64}, curve::Curve2D)

Compute all intersection points between infinite 2D line and `curve` (Segment2F or CircleArc2F).
Returns a vector of `(point::Point2F, t_curve::Float64)` where `t_curve ∈ [0, 1]` is the local curve parameter.
"""
function intersect_line_curve(line::Line{2,Float64}, seg::Segment2F)
    res = Tuple{Point2F,Float64}[]
    # Line 1: line.p + s * line.u
    # Line 2 (segment): seg.p + t * (seg.q - seg.p)
    dp = seg.q - seg.p
    det_m = line.u[1] * (-dp[2]) - line.u[2] * (-dp[1])
    if abs(det_m) < 1e-12
        return res
    end
    diff = seg.p - line.p
    s = (diff[1] * (-dp[2]) - diff[2] * (-dp[1])) / det_m
    t = (line.u[1] * diff[2] - line.u[2] * diff[1]) / det_m
    if 0.0 <= t <= 1.0
        pt = at(seg, t)
        push!(res, (pt, t))
    end
    return res
end

function intersect_line_curve(line::Line{2,Float64}, arc::CircleArc2F)
    res = Tuple{Point2F,Float64}[]
    # Distance from arc.center to line
    # Projection of center onto line:
    u_norm_sq = dot(line.u, line.u)
    u_norm_sq == 0.0 && return res
    t_proj = dot(arc.center - line.p, line.u) / u_norm_sq
    proj = line.p + t_proj * line.u
    d_sq = dist_sq(arc.center, proj)
    r_sq = arc.radius^2

    if d_sq > r_sq + 1e-12
        return res
    end

    h_sq = max(0.0, r_sq - d_sq)
    h = sqrt(h_sq)
    u_unit = line.u / sqrt(u_norm_sq)

    candidate_pts = Point2F[]
    if h < 1e-12
        push!(candidate_pts, proj)
    else
        push!(candidate_pts, proj + h * u_unit)
        push!(candidate_pts, proj - h * u_unit)
    end

    d_theta = arc.theta2 - arc.theta1
    if d_theta < 0
        d_theta = mod(d_theta, 2pi)
    end

    for pt in candidate_pts
        angle = atan(pt[2] - arc.center[2], pt[1] - arc.center[1])
        rel_angle = mod(angle - arc.theta1, 2pi)
        if rel_angle <= d_theta + 1e-9 || (2pi - rel_angle) <= 1e-9
            local_t = clamp(rel_angle / d_theta, 0.0, 1.0)
            # Verify actual point matches arc point
            actual_pt = at(arc, local_t)
            push!(res, (actual_pt, local_t))
        end
    end
    return res
end

"""
    tangent_intersections_cp(CP::CurvePolygon2D, Q, alpha::Real)

Compute the tangent line to `Q` at angle `alpha`, tangency point `q`, and the two intersection points `u` (right turn from O->q->u) and `v` (left turn from O->q->v) on `CP`.
Returns `(q, u, t_u, v, t_v)`.
"""
function tangent_intersections_cp(CP::CurvePolygon2D, Q, alpha::Real)
    line, q = tangent_line(Q, alpha)
    plens = prefix_lengths(CP)
    total_len = plens[end]

    intersections = Tuple{Point2F,Float64}[] # (point, parametric_location_on_CP)
    n_curves = length(CP.curves)
    for i = 1:n_curves
        curve = CP.curves[i]
        pts_t = intersect_line_curve(line, curve)
        for (pt, local_t) in pts_t
            cp_t = (plens[i] + local_t * (plens[i+1] - plens[i])) / total_len
            push!(intersections, (pt, cp_t))
        end
    end

    # Classify into right turn (u) and left turn (v)
    # Along tangent line L(t) = q + t*line.u:
    # line.u is orthogonal to outer normal v = direction(alpha), pointing counter-clockwise (left).
    # Since q is on boundary of Q, O -> q -> (q + t*line.u) turns:
    # - Left for t > 0 (v, after_tangent_to_polygon)
    # - Right for t < 0 (u, before_tangent_to_polygon)
    u_pt, u_t = q, 0.0
    v_pt, v_t = q, 0.0

    min_t_pos = Inf  # smallest t > 0 (closest intersection in +line.u direction -> v)
    max_t_neg = -Inf # largest t < 0 (closest intersection in -line.u direction -> u)

    for (pt, cp_t_val) in intersections
        # Express pt - q as t_line * line.u
        t_line = dot(pt - q, line.u) / dot(line.u, line.u)
        if t_line < -1e-9 # Right turn -> u
            if t_line > max_t_neg
                max_t_neg = t_line
                u_pt = pt
                u_t = cp_t_val
            end
        elseif t_line > 1e-9 # Left turn -> v
            if t_line < min_t_pos
                min_t_pos = t_line
                v_pt = pt
                v_t = cp_t_val
            end
        end
    end

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

Test whether query point `q` lies strictly inside the 2D curved polygon `cp` using ray casting.
"""
function is_inside(q::Point{2,S}, cp::CurvePolygon2D{T}) where {S,T}
    angles = (0.0, 0.3141592653589793, 0.7853981633974483, 1.2345678901234567, 2.1543210987654321)
    for ang in angles
        dir = Point{2,Float64}(cos(ang), sin(ang))
        ray = Line(Point{2,Float64}(q), dir)
        valid = true
        crossings = 0
        for c in cp.curves
            intersections = intersect_line_curve(ray, c)
            for (pt, s) in intersections
                dist_along_ray = dot(pt - Point{2,Float64}(q), dir)
                if abs(dist_along_ray) < 1e-11
                    return true # Point is on the boundary
                end
                if dist_along_ray > 0.0
                    if abs(s) < 1e-8 || abs(s - 1.0) < 1e-8
                        valid = false
                        break
                    end
                    crossings += 1
                end
            end
            if !valid
                break
            end
        end
        if valid
            return isodd(crossings)
        end
    end
    # Fallback to winding number
    return _is_inside_winding(q, cp)
end

function _is_inside_winding(q::Point{2,S}, cp::CurvePolygon2D{T}, n_per_curve::Int=64) where {S,T}
    total_angle = 0.0
    for c in cp.curves
        p_prev = at(c, 0.0)
        v_prev = p_prev - q
        for k in 1:n_per_curve
            s = k / n_per_curve
            p_curr = at(c, s)
            v_curr = p_curr - q
            cp_cross = v_prev[1] * v_curr[2] - v_prev[2] * v_curr[1]
            dp_dot = v_prev[1] * v_curr[1] + v_prev[2] * v_curr[2]
            total_angle += atan(cp_cross, dp_dot)
            v_prev = v_curr
        end
    end
    return abs(total_angle) > pi
end


