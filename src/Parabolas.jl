###############################################
### Parabola and ParabolaArc types

"""
    Parabola{T}

A parabola in 2D defined by its `focus` point and `directrix` line.
The parabola is the set of points `x` such that `dist(x, focus) = dist(x, directrix)`.
"""
struct Parabola{T} <: AbsCurve2D{T}
    focus::Point{2,T}
    directrix::Line{2,T}
end

"""
    Parabola(point::Point{2,T}, line::Line{2,T})

Construct a parabola from a focus point and a directrix line.
"""

"""
    bisector(line::Line{2,T}, point::Point{2,T})

Construct the bisector parabola of a line and a point.
This is the set of points equidistant from the line and the point.
Returns a `Parabola(point, line)`.
"""
function bisector(line::Line{2,T}, point::Point{2,T}) where {T}
    return Parabola(point, line)
end

"""
    vertex(h::Parabola)

Return the vertex of the parabola (midpoint between focus and its projection onto the directrix).
"""
function vertex(h::Parabola{T}) where {T}
    n = Point{2,T}(-h.directrix.u.y, h.directrix.u.x)
    n = n / norm(n)
    d_line = dot(h.directrix.p, n)
    d_focus = dot(h.focus, n)
    p = abs(d_focus - d_line) / 2
    if d_focus >= d_line
        return h.focus - p * n
    else
        return h.focus + p * n
    end
end

"""
    p_coeff(h::Parabola)

Return the focal parameter `p = dist(focus, directrix) / 2`.
"""
function p_coeff(h::Parabola{T}) where {T}
    n = Point{2,T}(-h.directrix.u.y, h.directrix.u.x)
    n = n / norm(n)
    d_line = dot(h.directrix.p, n)
    d_focus = dot(h.focus, n)
    return abs(d_focus - d_line) / 2
end

"""
    axis_direction(h::Parabola)

Return the unit vector from the directrix toward the focus (the axis direction).
"""
function axis_direction(h::Parabola{T}) where {T}
    n = Point{2,T}(-h.directrix.u.y, h.directrix.u.x)
    n = n / norm(n)
    d_line = dot(h.directrix.p, n)
    d_focus = dot(h.focus, n)
    return d_focus >= d_line ? n : -n
end

"""
    at(h::Parabola, t::Real)

Return a point on the parabola at parameter `t` (t ∈ ℝ).
Uses the parametric form: `x = V + t²/(4p)·n + t·u`,
where `V` is the vertex, `n` is the axis direction, and `u` is a perpendicular unit vector.
"""
function at(h::Parabola{T}, t::Real) where {T}
    V = vertex(h)
    n = axis_direction(h)
    u = Point{2,T}(-n.y, n.x)
    p = p_coeff(h)
    return V + (Float64(t)^2 / (4 * p)) * n + Float64(t) * u
end

Base.eltype(::Parabola{T}) where {T} = T

geom_length(::Parabola) = Inf

function Base.show(io::IO, h::Parabola{T}) where {T}
    print(io, "Parabola(focus=", h.focus, ", directrix=", h.directrix, ")")
end

"""
    intersect_line_curve(line::Line{2,Float64}, h::Parabola)

Compute all intersection points between an infinite 2D line and a parabola.
Returns a vector of `(point::Point2F, t_curve::Float64)`.
"""
function intersect_line_curve(line::Line{2,Float64}, h::Parabola)
    res = Tuple{Point2F,Float64}[]
    V = vertex(h)
    n = axis_direction(h)
    u = Point{2,Float64}(-n.y, n.x)
    p = p_coeff(h)

    lp = Point{2,Float64}(line.p)
    lu = Point{2,Float64}(line.u)

    A = lu.x * n.x + lu.y * n.y
    B = lu.x * u.x + lu.y * u.y
    C = (lp.x - V.x) * n.x + (lp.y - V.y) * n.y
    D = (lp.x - V.x) * u.x + (lp.y - V.y) * u.y

    if abs(B) < 1e-12
        if abs(A) < 1e-12
            return res
        end
        t_curve = D
        s = (D^2 / (4 * p) - C) / A
        pt = lp + s * lu
        push!(res, (pt, t_curve))
        return res
    end

    a = 1.0
    b = -4.0 * p * A / B
    c = 4.0 * p * (A * D - B * C) / B

    disc = b^2 - 4 * a * c
    if disc < -1e-12
        return res
    end
    disc = max(0.0, disc)
    sqrt_disc = sqrt(disc)

    for sgn in (1.0, -1.0)
        t_curve = (-b + sgn * sqrt_disc) / (2 * a)
        s = (t_curve - D) / B
        pt = lp + s * lu
        push!(res, (pt, t_curve))
    end

    return res
end

"""
    distance(x::Point{2,T}, h::Parabola{T})

Compute the minimum distance from point `x` to the parabola `h`.
Uses numerical sampling.
"""
function distance(x::Point{2,T}, h::Parabola{T}) where {T}
    best = Inf
    for t in range(-5.0, 5.0, length=500)
        p = at(h, t)
        d = dist(x, p)
        if d < best
            best = d
        end
    end
    return best
end

###############################################
### ParabolaArc type

"""
    ParabolaArc{T}

A bounded segment of a 2D parabola between parameters `t1` and `t2`.
"""
struct ParabolaArc{T} <: AbsCurve2D{T}
    parabola::Parabola{T}
    t1::T
    t2::T
end

const ParabolaArc2F = ParabolaArc{Float64}

"""
    ParabolaArc(parabola::Parabola{T}, t1::Real, t2::Real)

Construct a `ParabolaArc` from a parabola and parameter bounds `t1`, `t2`.
"""
function ParabolaArc(parabola::Parabola{T}, t1::Real, t2::Real) where {T}
    return ParabolaArc{T}(parabola, T(t1), T(t2))
end

"""
    ParabolaArc(parabola::Parabola{T}, p1::Point{2,T}, p2::Point{2,T})

Construct a `ParabolaArc` from a parabola and two endpoints `p1`, `p2` on the parabola.
"""
function ParabolaArc(parabola::Parabola{T}, p1::Point{2,T}, p2::Point{2,T}) where {T}
    V = vertex(parabola)
    n = axis_direction(parabola)
    u = Point{2,T}(-n.y, n.x)
    t1 = dot(p1 - V, u)
    t2 = dot(p2 - V, u)
    return ParabolaArc{T}(parabola, t1, t2)
end

Base.eltype(::ParabolaArc{T}) where {T} = T

"""
    at(arc::ParabolaArc, s::Real)

Return a point on `arc` at relative parameter `s` in `[0, 1]`.
`s=0` corresponds to `t1` and `s=1` corresponds to `t2`.
"""
function at(arc::ParabolaArc{T}, s::Real) where {T}
    t = arc.t1 + T(s) * (arc.t2 - arc.t1)
    return at(arc.parabola, t)
end

"""
    geom_length(arc::ParabolaArc)

Compute the exact arc length of the parabolic arc between parameters `t1` and `t2`.
"""
function geom_length(arc::ParabolaArc{T}) where {T}
    p = p_coeff(arc.parabola)
    u1 = arc.t1 / (2 * p)
    u2 = arc.t2 / (2 * p)
    F(u) = p * (u * sqrt(1 + u^2) + asinh(u))
    return abs(F(u2) - F(u1))
end

function Base.show(io::IO, arc::ParabolaArc{T}) where {T}
    print(io, "ParabolaArc(", arc.parabola, ", t1=", arc.t1, ", t2=", arc.t2, ")")
end

"""
    intersect_line_curve(line::Line{2,Float64}, arc::ParabolaArc)

Compute all intersection points between an infinite 2D line and a parabolic arc.
Returns a vector of `(point::Point2F, s::Float64)` where `s ∈ [0, 1]` is the local arc parameter.
"""
function intersect_line_curve(line::Line{2,Float64}, arc::ParabolaArc)
    res = Tuple{Point2F,Float64}[]
    pts = intersect_line_curve(line, arc.parabola)
    min_t = min(arc.t1, arc.t2)
    max_t = max(arc.t1, arc.t2)
    dt = arc.t2 - arc.t1
    for (pt, t_curve) in pts
        if min_t - 1e-11 <= t_curve <= max_t + 1e-11
            s = abs(dt) < 1e-12 ? 0.0 : clamp((t_curve - arc.t1) / dt, 0.0, 1.0)
            push!(res, (pt, s))
        end
    end
    return res
end

export Parabola, ParabolaArc, ParabolaArc2F, vertex, p_coeff, axis_direction
