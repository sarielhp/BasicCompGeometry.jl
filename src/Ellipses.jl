# Ellipses and Elliptic Arcs in 2D

"""
    Ellipse{T} <: AbsCurve2D{T}

A 2D ellipse defined by its `center` (`Point{2, T}`), semi-major axis `r_major` (`T`),
semi-minor axis `r_minor` (`T`), and counter-clockwise orientation `angle` in radians.
"""
struct Ellipse{T} <: AbsCurve2D{T}
    center::Point{2,T}
    r_major::T
    r_minor::T
    angle::T

    function Ellipse{T}(center::Point{2,T}, r_major::T, r_minor::T, angle::T) where {T}
        @assert r_major >= zero(T) "Semi-major axis must be non-negative"
        @assert r_minor >= zero(T) "Semi-minor axis must be non-negative"
        return new{T}(center, r_major, r_minor, angle)
    end
end

"""
    Ellipse2F

Type alias for an ellipse with `Float64` coordinates.
"""
const Ellipse2F = Ellipse{Float64}

"""
    Ellipse(center, r_major, r_minor, angle=0.0)

Construct an `Ellipse` with given parameters, promoting numeric types.
"""
function Ellipse(center::Point{2,T}, r_major::Real, r_minor::Real, angle::Real = 0.0) where {T}
    S = promote_type(T, typeof(r_major), typeof(r_minor), typeof(angle))
    c = Point{2,S}(S(center[1]), S(center[2]))
    return Ellipse{S}(c, S(r_major), S(r_minor), S(angle))
end

function Ellipse(center::NTuple{2,Real}, r_major::Real, r_minor::Real, angle::Real = 0.0)
    S = promote_type(typeof(center[1]), typeof(center[2]), typeof(r_major), typeof(r_minor), typeof(angle))
    return Ellipse{S}(Point{2,S}(center[1], center[2]), S(r_major), S(r_minor), S(angle))
end

center(e::Ellipse) = e.center
r_major(e::Ellipse) = e.r_major
r_minor(e::Ellipse) = e.r_minor
Base.angle(e::Ellipse) = e.angle
rotation_angle(e::Ellipse) = e.angle
area(e::Ellipse) = π * e.r_major * e.r_minor

"""
    point_on(e::Ellipse, t::Real)

Parametric point on the ellipse at parameter `t` in radians.
"""
function point_on(e::Ellipse{T}, t::Real) where {T}
    ct, st = cos(T(t)), sin(T(t))
    ca, sa = cos(e.angle), sin(e.angle)
    x = e.center[1] + e.r_major * ct * ca - e.r_minor * st * sa
    y = e.center[2] + e.r_major * ct * sa + e.r_minor * st * ca
    return Point{2,T}(x, y)
end

"""
    is_inside(p::Point{2}, e::Ellipse)

Test whether point `p` is strictly inside or on the boundary of `e`.
"""
function is_inside(p::Point{2}, e::Ellipse{T}) where {T}
    dx = p[1] - e.center[1]
    dy = p[2] - e.center[2]
    ca, sa = cos(e.angle), sin(e.angle)
    u = (dx * ca + dy * sa) / e.r_major
    v = (-dx * sa + dy * ca) / e.r_minor
    return (u * u + v * v) <= one(T) + 1e-12
end

"""
    bbox(e::Ellipse)

Compute the exact axis-aligned bounding box of a rotated ellipse.
"""
function bbox(e::Ellipse{T}) where {T}
    ca, sa = cos(e.angle), sin(e.angle)
    half_w = sqrt((e.r_major * ca)^2 + (e.r_minor * sa)^2)
    half_h = sqrt((e.r_major * sa)^2 + (e.r_minor * ca)^2)
    bl = Point{2,T}(e.center[1] - half_w, e.center[2] - half_h)
    tr = Point{2,T}(e.center[1] + half_w, e.center[2] + half_h)
    return BBox(bl, tr)
end

BBox(e::Ellipse) = bbox(e)

"""
    sample_uniformly(e::Ellipse, n::Int)

Sample `n` points along the perimeter of the ellipse, returning a `PntSeq{2, T}`.
"""
function sample_uniformly(e::Ellipse{T}, n::Int) where {T}
    @assert n >= 3 "Sample count must be at least 3"
    pts = [point_on(e, 2π * i / n) for i in 0:(n-1)]
    return PntSeq(pts)
end

function Base.show(io::IO, e::Ellipse{T}) where {T}
    print(io, "Ellipse(", e.center, ", a=", e.r_major, ", b=", e.r_minor, ", θ=", round(e.angle, digits=3), ")")
end

# -----------------------------------------------------------------------------
# Elliptic Arc
# -----------------------------------------------------------------------------

"""
    EllipticArc{T} <: AbsCurve2D{T}

A 2D elliptic arc defined by an underlying `ellipse` (`Ellipse{T}`) and counter-clockwise
angular span from `alpha1` to `alpha2` (in radians).
"""
struct EllipticArc{T} <: AbsCurve2D{T}
    ellipse::Ellipse{T}
    alpha1::T
    alpha2::T
end

"""
    EllipticArc2F

Type alias for an elliptic arc with `Float64` coordinates.
"""
const EllipticArc2F = EllipticArc{Float64}

"""
    EllipticArc(ellipse, alpha1, alpha2)
    EllipticArc(center, r_major, r_minor, angle, alpha1, alpha2)

Construct an `EllipticArc`.
"""
function EllipticArc(e::Ellipse{T}, alpha1::Real, alpha2::Real) where {T}
    S = promote_type(T, typeof(alpha1), typeof(alpha2))
    e_s = Ellipse{S}(Point{2,S}(e.center[1], e.center[2]), S(e.r_major), S(e.r_minor), S(e.angle))
    return EllipticArc{S}(e_s, S(alpha1), S(alpha2))
end

function EllipticArc(center::Point{2}, r_major::Real, r_minor::Real, angle::Real, alpha1::Real, alpha2::Real)
    e = Ellipse(center, r_major, r_minor, angle)
    return EllipticArc(e, alpha1, alpha2)
end

center(arc::EllipticArc) = arc.ellipse.center
r_major(arc::EllipticArc) = arc.ellipse.r_major
r_minor(arc::EllipticArc) = arc.ellipse.r_minor
angle(arc::EllipticArc) = arc.ellipse.angle
alpha1(arc::EllipticArc) = arc.alpha1
alpha2(arc::EllipticArc) = arc.alpha2

"""
    point_on(arc::EllipticArc, t::Real)

Evaluate point on the arc for parameter `t ∈ [0, 1]`.
"""
function point_on(arc::EllipticArc{T}, t::Real) where {T}
    span = arc.alpha2 >= arc.alpha1 ? (arc.alpha2 - arc.alpha1) : (arc.alpha2 - arc.alpha1 + 2π)
    angle_param = arc.alpha1 + T(t) * span
    return point_on(arc.ellipse, angle_param)
end

"""
    start_point(arc::EllipticArc)
    end_point(arc::EllipticArc)

End points of the elliptic arc.
"""
start_point(arc::EllipticArc) = point_on(arc.ellipse, arc.alpha1)
end_point(arc::EllipticArc) = point_on(arc.ellipse, arc.alpha2)

function _angle_in_ccw_span(a::Real, a1::Real, a2::Real)
    two_pi = 2 * π
    da = mod(a - a1, two_pi)
    dspan = mod(a2 - a1, two_pi)
    dspan == 0 && return false
    return da <= dspan
end

"""
    bbox(arc::EllipticArc)

Compute exact axis-aligned bounding box of an elliptic arc.
"""
function bbox(arc::EllipticArc{T}) where {T}
    pts = Point{2,T}[start_point(arc), end_point(arc)]
    e = arc.ellipse
    ca, sa = cos(e.angle), sin(e.angle)

    # Critical angles for x: dx/dt = 0 => tan(t) = -(b*sin(θ))/(a*cos(θ))
    denom_x = e.r_major * ca
    if abs(denom_x) > 1e-12
        t_x = atan(-e.r_minor * sa, denom_x)
        for t in (t_x, t_x + π)
            _angle_in_ccw_span(t, arc.alpha1, arc.alpha2) && push!(pts, point_on(e, t))
        end
    end

    # Critical angles for y: dy/dt = 0 => tan(t) = (b*cos(θ))/(a*sin(θ))
    denom_y = e.r_major * sa
    if abs(denom_y) > 1e-12
        t_y = atan(e.r_minor * ca, denom_y)
        for t in (t_y, t_y + π)
            _angle_in_ccw_span(t, arc.alpha1, arc.alpha2) && push!(pts, point_on(e, t))
        end
    end

    min_x = minimum(p[1] for p in pts)
    max_x = maximum(p[1] for p in pts)
    min_y = minimum(p[2] for p in pts)
    max_y = maximum(p[2] for p in pts)
    return BBox(Point{2,T}(min_x, min_y), Point{2,T}(max_x, max_y))
end

"""
    sample_uniformly(arc::EllipticArc, n::Int)

Sample `n` points along the arc, returning a `PntSeq{2, T}`.
"""
function sample_uniformly(arc::EllipticArc{T}, n::Int) where {T}
    @assert n >= 2 "Sample count must be at least 2"
    pts = [point_on(arc, T(i) / (n - 1)) for i in 0:(n-1)]
    return PntSeq(pts)
end

function Base.show(io::IO, arc::EllipticArc{T}) where {T}
    print(io, "EllipticArc(", arc.ellipse, ", α1=", round(arc.alpha1, digits=3), ", α2=", round(arc.alpha2, digits=3), ")")
end

BBox(arc::EllipticArc) = bbox(arc)

export Ellipse, Ellipse2F, EllipticArc, EllipticArc2F
export r_major, r_minor, alpha1, alpha2, start_point, end_point, area, bbox
