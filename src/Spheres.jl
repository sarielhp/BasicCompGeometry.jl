###############################################
### Sphere type
###############################################

"""
    Sphere{D, T}

A sphere in `D` dimensions defined by its `center` (a `Point{D, T}`) and scalar `radius` of type `T`.
"""
struct Sphere{D,T}
    center::Point{D,T}
    radius::T
end

"""
    Sphere(center::Point{D, T}, radius::Real)

Construct a `Sphere` with center `center` and radius converted to type `T`.
"""
function Sphere(center::Point{D,T}, radius::Real) where {D,T}
    return Sphere{D,T}(center, T(radius))
end

"""
    Sphere(center::AbstractVector{S}, radius::Real)

Construct a `Sphere` from a vector coordinate center.
"""
function Sphere(center::AbstractVector{S}, radius::Real) where {S}
    T = promote_type(S, typeof(radius))
    p = Point{length(center),T}(center...)
    return Sphere(p, T(radius))
end

"""
    Sphere(center::NTuple{D, S}, radius::Real)

Construct a `Sphere` from a tuple coordinate center.
"""
function Sphere(center::NTuple{D,S}, radius::Real) where {D,S}
    T = promote_type(S, typeof(radius))
    p = Point{D,T}(center...)
    return Sphere(p, T(radius))
end

"""
    Circle{T}

Type alias for a 2D sphere (`Sphere{2, T}`).
"""
const Circle{T} = Sphere{2,T}

"""
    Circle2F

Type alias for a 2D circle with `Float64` coordinates and radius.
"""
const Circle2F = Circle{Float64}

"""
    Circle2I

Type alias for a 2D circle with `Int64` coordinates and radius.
"""
const Circle2I = Circle{Int64}

"""
    Sphere2F

Type alias for a 2D sphere with `Float64` coordinates and radius.
"""
const Sphere2F = Sphere{2,Float64}

"""
    Sphere3F

Type alias for a 3D sphere with `Float64` coordinates and radius.
"""
const Sphere3F = Sphere{3,Float64}

"""
    Circle(center, radius)

Construct a 2D circle (`Sphere{2, T}`).
"""
Circle(center::Point{2,T}, radius::Real) where {T} = Sphere(center, radius)
Circle(center::AbstractVector, radius::Real) = Sphere(center, radius)
Circle(center::NTuple{2,Real}, radius::Real) = Sphere(center, radius)

function Base.show(io::IO, s::Sphere{D,T}) where {D,T}
    if D == 2
        print(io, "Circle(", s.center, ", r=", s.radius, ")")
    else
        print(io, "Sphere{", D, "}(", s.center, ", r=", s.radius, ")")
    end
end

"""
    is_inside(p::Point{D, T1}, s::Sphere{D, T2})

Return `true` if point `p` is inside or on the boundary of sphere `s`.
"""
@inline function is_inside(p::Point{D,T1}, s::Sphere{D,T2}) where {D,T1,T2}
    return dist_sq(p, s.center) <= s.radius^2
end

"""
    dist(p::Point{D, T1}, s::Sphere{D, T2})

Calculate the Euclidean distance from point `p` to the surface of sphere `s`.
If `p` is inside `s`, returns 0.0.
"""
@inline function dist(p::Point{D,T1}, s::Sphere{D,T2}) where {D,T1,T2}
    d_center = dist(p, s.center)
    return max(0.0, d_center - Float64(s.radius))
end

###############################################
### Circular Inversion Functions
###############################################

"""
    invert(p::Point{D, T}, ref::Sphere{D, S} = Sphere(zeros(Point{D, Float64}), 1.0))

Invert point `p` with respect to reference sphere `ref` centered at `O` with radius `R`.
The inverse point is `P' = O + (R^2 / ||p - O||^2) * (p - O)`.
"""
function invert(
    p::Point{D,T},
    ref::Sphere{D,S} = Sphere(zeros(Point{D,Float64}), 1.0),
) where {D,T,S}
    O = Point{D,Float64}(ref.center)
    R2 = Float64(ref.radius)^2
    diff = Point{D,Float64}(p) - O
    d2 = sum(abs2, diff)
    d2 == 0.0 && error("Cannot invert the center of the reference sphere (maps to point at infinity).")
    return O + (R2 / d2) * diff
end

"""
    invert(l::Line{2, T}, ref::Sphere{2, S} = Sphere(zeros(Point{2, Float64}), 1.0); tol=1e-12)

Invert a 2D line `l` with respect to reference circle `ref`.
Returns either a `Line{2, Float64}` (if `l` passes through the center of `ref`)
or a `Sphere{2, Float64}` (Circle) passing through the center of `ref`.
"""
function invert(
    l::Line{2,T},
    ref::Sphere{2,S} = Sphere(zeros(Point{2,Float64}), 1.0);
    tol::Real = 1e-12,
) where {T,S}
    O = Point2F(ref.center)
    R = Float64(ref.radius)
    R2 = R^2

    # Find orthogonal projection P0 of O onto line l
    u = Point2F(l.u)
    u_norm_sq = sum(abs2, u)
    u_norm_sq == 0.0 && error("Line direction vector u cannot be zero.")

    t0 = dot(O - Point2F(l.p), u) / u_norm_sq
    P0 = Point2F(l.p) + t0 * u
    d = dist(O, P0)

    if d <= tol * R
        # Line passes through center of inversion O -> inverts to itself
        return Line{2,Float64}(Point2F(l.p), u)
    else
        # Line does not pass through O -> inverts to a circle passing through O.
        # P0 inverts to P0', and OP0' is a diameter of the inverted circle.
        P0_inv = O + (R2 / (d^2)) * (P0 - O)
        center_inv = (O + P0_inv) / 2.0
        radius_inv = R2 / (2.0 * d)
        return Sphere{2,Float64}(center_inv, radius_inv)
    end
end

"""
    invert(c::Sphere{2, T}, ref::Sphere{2, S} = Sphere(zeros(Point{2, Float64}), 1.0); tol=1e-12)

Invert a 2D circle `c` with respect to reference circle `ref`.
Returns either a `Line{2, Float64}` (if circle `c` passes through the center of `ref`)
or a `Sphere{2, Float64}` (Circle).
"""
function invert(
    c::Sphere{2,T},
    ref::Sphere{2,S} = Sphere(zeros(Point{2,Float64}), 1.0);
    tol::Real = 1e-12,
) where {T,S}
    O = Point2F(ref.center)
    R = Float64(ref.radius)
    R2 = R^2

    A = Point2F(c.center)
    r = Float64(c.radius)

    d_O = dist(A, O)

    if abs(d_O - r) <= tol * max(R, r)
        # Circle c passes through center of inversion O -> inverts to a line.
        diff = A - O
        B_inv = O + (R2 / (2.0 * r^2)) * diff
        dir = Point2F(-diff[2], diff[1]) # Vector perpendicular to diff
        return Line{2,Float64}(B_inv, dir)
    else
        # Circle c does not pass through O -> inverts to another circle.
        denom = d_O^2 - r^2
        center_inv = O + (R2 / denom) * (A - O)
        radius_inv = (R2 * r) / abs(denom)
        return Sphere{2,Float64}(center_inv, radius_inv)
    end
end

export Sphere, Circle, Circle2F, Circle2I, Sphere2F, Sphere3F
export CircleArc, CircleArc2F
export invert

###############################################
### CircleArc type
###############################################

"""
    CircleArc{T}

A 2D circular arc defined by a `center` (`Point{2, T}`), scalar `radius` of type `T`, 
and counter-clockwise angular span from `theta1` to `theta2` (in radians).
"""
struct CircleArc{T}
    center::Point{2,T}
    radius::T
    theta1::T
    theta2::T
end

"""
    CircleArc2F

Type alias for a 2D circular arc with `Float64` coordinates and angles.
"""
const CircleArc2F = CircleArc{Float64}

"""
    CircleArc(center, radius, theta1, theta2)

Construct a `CircleArc`.
"""
function CircleArc(center::Point{2,T}, radius::Real, theta1::Real, theta2::Real) where {T}
    return CircleArc{T}(center, T(radius), T(theta1), T(theta2))
end

function CircleArc(center::NTuple{2,Real}, radius::Real, theta1::Real, theta2::Real)
    T = promote_type(typeof(center[1]), typeof(center[2]), typeof(radius), typeof(theta1), typeof(theta2))
    return CircleArc{T}(Point{2,T}(center), T(radius), T(theta1), T(theta2))
end

function Base.show(io::IO, arc::CircleArc{T}) where {T}
    print(io, "CircleArc(", arc.center, ", r=", arc.radius, ", θ1=", arc.theta1, ", θ2=", arc.theta2, ")")
end

"""
    invert(seg::Segment{2, T}, ref::Sphere{2, S} = Sphere(zeros(Point{2, Float64}), 1.0); tol=1e-12)

Invert a 2D line segment `seg` with respect to reference circle `ref`.
Returns either a `CircleArc{Float64}` (if the supporting line does not pass through the center of `ref`)
or a `Segment{2, Float64}` (if the supporting line passes through the center of `ref`).
"""
function invert(
    seg::Segment{2,T},
    ref::Sphere{2,S} = Sphere(zeros(Point{2,Float64}), 1.0);
    tol::Real = 1e-12,
) where {T,S}
    O = Point2F(ref.center)
    P1 = Point2F(seg.p)
    P2 = Point2F(seg.q)

    dir = P2 - P1
    sum(abs2, dir) == 0.0 && return Segment(invert(P1, ref), invert(P2, ref))
    l = Line(P1, dir)

    inv_l = invert(l, ref; tol=tol)

    if inv_l isa Sphere{2,Float64}
        A_prime = inv_l.center
        r_prime = inv_l.radius

        P1_prime = invert(P1, ref)
        P2_prime = invert(P2, ref)

        th1 = atan(P1_prime[2] - A_prime[2], P1_prime[1] - A_prime[1])
        th2 = atan(P2_prime[2] - A_prime[2], P2_prime[1] - A_prime[1])
        thO = atan(O[2] - A_prime[2], O[1] - A_prime[1])

        d_th2 = mod(th2 - th1, 2pi)
        d_thO = mod(thO - th1, 2pi)

        if 0.0 < d_thO < d_th2
            th_start = th2
            th_end = th1 < th2 ? th1 + 2pi : th1
        else
            th_start = th1
            th_end = th2 < th1 ? th2 + 2pi : th2
        end

        return CircleArc{Float64}(A_prime, r_prime, th_start, th_end)
    else
        return Segment(invert(P1, ref), invert(P2, ref))
    end
end
