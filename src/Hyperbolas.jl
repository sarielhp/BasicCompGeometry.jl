###############################################
### Hyperbola type

"""
    Hyperbola{T}

A hyperbola in 2D defined by its two foci `f1`, `f2` and constant distance difference `d`.
The hyperbola is the set of points `x` such that `|dist(x, f1) - dist(x, f2)| = d`.

The transverse axis is along the direction from `f1` to `f2`.
The center is the midpoint of the two foci.
"""
struct Hyperbola{T} <: AbsCurve2D{T}
    f1::Point{2,T}
    f2::Point{2,T}
    d::T
end

"""
    Hyperbola(center, a, b, angle)

Construct a hyperbola from its `center`, semi-transverse axis `a`, semi-conjugate axis `b`,
and rotation `angle` of the transverse axis from the x-axis.
"""
function Hyperbola(center::Point{2,T}, a::Real, b::Real, angle::Real) where {T}
    c = sqrt(a^2 + b^2)
    f_offset = Point{2,T}(c * cos(angle), c * sin(angle))
    f1 = center - f_offset
    f2 = center + f_offset
    d = 2 * a
    return Hyperbola{T}(f1, f2, T(d))
end

function Hyperbola(center::AbstractVector, a::Real, b::Real, angle::Real)
    T = promote_type(eltype(center), typeof(a), typeof(b), typeof(angle))
    c = Point{2,T}(center...)
    return Hyperbola(c, T(a), T(b), T(angle))
end

"""
    Hyperbola(f1, f2, d)

Construct a hyperbola from two foci `f1`, `f2` and constant distance difference `d`.
Requires `d < dist(f1, f2)`.
"""
function Hyperbola(f1::Point{2,T}, f2::Point{2,T}, d::Real) where {T}
    @assert d < dist(f1, f2) "Constant difference d must be less than distance between foci"
    return Hyperbola{T}(f1, f2, T(d))
end

"""
    Hyperbola(line::Line{2,T}, point::Point{2,T})

Construct the bisector of a line and a point.
The bisector is the set of points `x` such that `dist(x, line) = dist(x, point)`.
This is a parabola (focus-directrix definition), not a hyperbola.

Deprecated: use `Parabola(point, line)` instead.
"""
function Hyperbola(line::Line{2,T}, point::Point{2,T}) where {T}
    Base.depwarn("`Hyperbola(line, point)` is deprecated, use `Parabola(point, line)` instead.", :Hyperbola)
    return Parabola(point, line)
end

"""
    center(h::Hyperbola)

Return the center of the hyperbola (midpoint of the two foci).
"""
center(h::Hyperbola{T}) where {T} = (h.f1 + h.f2) / 2

"""
    a_coeff(h::Hyperbola)

Return the semi-transverse axis length `a = d/2`.
"""
a_coeff(h::Hyperbola{T}) where {T} = h.d / 2

"""
    c_coeff(h::Hyperbola)

Return the focal distance `c = dist(f1, f2)/2`.
"""
c_coeff(h::Hyperbola{T}) where {T} = dist(h.f1, h.f2) / 2

"""
    b_coeff(h::Hyperbola)

Return the semi-conjugate axis length `b = sqrt(c^2 - a^2)`.
"""
b_coeff(h::Hyperbola{T}) where {T} = sqrt(max(0.0, c_coeff(h)^2 - a_coeff(h)^2))

"""
    rotation_angle(h::Hyperbola)

Return the rotation angle of the transverse axis from the x-axis.
"""
rotation_angle(h::Hyperbola{T}) where {T} = atan(h.f2.y - h.f1.y, h.f2.x - h.f1.x)

"""
    eltype(::Hyperbola{T})

Return the element type `T` of the hyperbola.
"""
Base.eltype(::Hyperbola{T}) where {T} = T

"""
    at(h::Hyperbola, t::Real; branch::Int=1)

Return a point on the hyperbola at parameter `t`.
`branch=1` returns points on the branch closer to `f1`, `branch=2` returns points on the branch closer to `f2`.
The parameter `t` ranges over all real numbers.
Uses the parametric form: `x = center + a*cosh(t)*u + b*sinh(t)*v` for branch 1,
and `x = center - a*cosh(t)*u + b*sinh(t)*v` for branch 2,
where `u` is the unit vector along the transverse axis and `v` is the unit vector along the conjugate axis.
"""
function at(h::Hyperbola{T}, t::Real; branch::Int=1) where {T}
    c = center(h)
    a = a_coeff(h)
    b = b_coeff(h)
    ang = rotation_angle(h)
    u = Point{2,T}(cos(ang), sin(ang))
    v = Point{2,T}(-sin(ang), cos(ang))
    ct = cosh(Float64(t))
    st = sinh(Float64(t))
    if b < 1e-12
        if branch == 1
            return c + a * ct * u + a * st * v
        else
            return c - a * ct * u + a * st * v
        end
    else
        if branch == 1
            return c + a * ct * u + b * st * v
        else
            return c - a * ct * u + b * st * v
        end
    end
end

"""
    geom_length(h::Hyperbola)

Return `Inf` since a hyperbola is an unbounded curve.
"""
geom_length(::Hyperbola) = Inf

function Base.show(io::IO, h::Hyperbola{T}) where {T}
    print(io, "Hyperbola(f1=", h.f1, ", f2=", h.f2, ", d=", h.d, ")")
end

"""
    intersect_line_curve(line::Line{2,Float64}, h::Hyperbola)

Compute all intersection points between an infinite 2D line and a hyperbola.
Returns a vector of `(point::Point2F, t_curve::Float64)` where `t_curve` is the hyperbola parameter.
"""
function intersect_line_curve(line::Line{2,Float64}, h::Hyperbola)
    res = Tuple{Point2F,Float64}[]
    c = center(h)
    a = a_coeff(h)
    b = b_coeff(h)
    ang = rotation_angle(h)
    u = Point{2,Float64}(cos(ang), sin(ang))
    v = Point{2,Float64}(-sin(ang), cos(ang))

    lp = Point{2,Float64}(line.p)
    lu = Point{2,Float64}(line.u)

    A = lu.x * v.x + lu.y * v.y
    B = lu.x * u.x + lu.y * u.y
    C = (lp.x - c.x) * v.x + (lp.y - c.y) * v.y
    D = (lp.x - c.x) * u.x + (lp.y - c.y) * u.y

    if b < 1e-12
        A2 = 0.0
    else
        A2 = A^2 / b^2
    end
    B2 = B^2 / a^2
    denom = A2 - B2
    if abs(denom) < 1e-12
        return res
    end

    if b < 1e-12
        alpha = -(B * D) / a^2
        beta = -D^2 / a^2 + 1.0
    else
        alpha = (A * C) / b^2 - (B * D) / a^2
        beta = C^2 / b^2 - D^2 / a^2 + 1.0
    end

    disc = alpha^2 - denom * beta
    if disc < -1e-12
        return res
    end
    disc = max(0.0, disc)
    sqrt_disc = sqrt(disc)

    for sgn in (1.0, -1.0)
        t_line = (-alpha + sgn * sqrt_disc) / denom
        pt = lp + t_line * lu
        dx = (pt.x - c.x) * u.x + (pt.y - c.y) * u.y
        dy = (pt.x - c.x) * v.x + (pt.y - c.y) * v.y
        if b < 1e-12
            t_curve = dx / a
        else
            t_curve = asinh(dy / b)
        end
        push!(res, (pt, t_curve))
    end

    return res
end

"""
    intersect_hyperbolas(h1::Hyperbola, h2::Hyperbola)

Compute intersection points of two hyperbolas by solving the system of equations
`|dist(x, h1.f1) - dist(x, h1.f2)| = h1.d` and `|dist(x, h2.f1) - dist(x, h2.f2)| = h2.d`.
Uses numerical root-finding along a search grid.
Returns a vector of intersection points.
"""
function intersect_hyperbolas(h1::Hyperbola{T}, h2::Hyperbola{T}) where {T}
    pts = Point{2,T}[]
    c1 = center(h1)
    c2 = center(h2)
    r = max(dist(c1, c2) + max(c_coeff(h1), c_coeff(h2)), 10.0)
    for branch1 in (1, 2)
        for branch2 in (1, 2)
            for t in range(-3.0, 3.0, length=200)
                p1 = at(h1, t; branch=branch1)
                d1 = abs(dist(p1, h1.f1) - dist(p1, h1.f2)) - h1.d
                d2 = abs(dist(p1, h2.f1) - dist(p1, h2.f2)) - h2.d
                if abs(d1) < 0.01 && abs(d2) < 0.01
                    push!(pts, p1)
                end
            end
        end
    end
    unique!(pts)
    return filter!(p -> !any(isnan, p), pts)
end

"""
    distance(x::Point{2,T}, h::Hyperbola{T})

Compute the minimum distance from point `x` to the hyperbola `h`.
Uses numerical optimization by sampling points on both branches.
"""
function distance(x::Point{2,T}, h::Hyperbola{T}) where {T}
    best = Inf
    for branch in (1, 2)
        for t in range(-5.0, 5.0, length=500)
            p = at(h, t; branch=branch)
            d = dist(x, p)
            if d < best
                best = d
            end
        end
    end
    return best
end

export Hyperbola, center, a_coeff, b_coeff, c_coeff, rotation_angle, intersect_hyperbolas