###############################################
### Plane type

"""
    Plane{D, T}

A (d-1)-dimensional hyperplane (plane) in D dimensions,
orthogonal to the normal direction `n` and passing through point `p`.
The plane consists of all points `x` satisfying `dot(n, x - p) = 0`,
equivalently `dot(n, x) = dot(n, p)`.

In 2D, a Plane is the same as a Line.
"""
struct Plane{D,T}
    p::Point{D,T}
    n::Point{D,T}
end

"""
    Plane2F

A 2-dimensional plane (line) with `Float64` coordinates.
"""
const Plane2F = Plane{2,Float64}

"""
    polar(p)

Return the polar hyperplane (or line) of point `p` with respect to the unit sphere/circle.
In 2D, this returns a `Line{2, T}`. In D ≥ 3, this returns a `Plane{D, T}`.
"""
function polar(p::Point{D,T}) where {D,T}
    d = dot(p, p)
    pl = Plane{D,T}(p / d, p)
    return D == 2 ? Line(pl) : pl
end

"""
    polar(pl)

Return the polar point of plane `pl` with respect to the unit sphere.
The plane is defined by `dot(n, x) = dot(n, p)`, and the polar point is `n / dot(n, p)`.
"""
function polar(pl::Plane{D,T}) where {D,T}
    return pl.n / dot(pl.n, pl.p)
end

"""
    Plane(l::Line{2,T})

Construct a 2D `Plane` from a `Line`.
The line's direction vector `u` is converted to the plane's normal vector `n = (u_y, -u_x)`.
"""
Plane(l::Line{2,T}) where {T} = Plane{2,T}(l.p, Point{2,T}(l.u[2], -l.u[1]))

"""
    Line(pl::Plane{2,T})

Construct a 2D `Line` from a `Plane`.
The plane's normal vector `n` is converted to the line's direction vector `u = (-n_y, n_x)`.
"""
Line(pl::Plane{2,T}) where {T} = Line{2,T}(pl.p, Point{2,T}(-pl.n[2], pl.n[1]))

"""
    polar(l)

Return the polar point of the 2D line `l` with respect to the unit sphere.
This is the same as `polar(::Plane)` since in 2D a line is a (d-1)-dimensional subspace.
"""
function polar(l::Line{2,T}) where {T}
    return polar(Plane(l))
end

"""
    polar(::Type{Line}, p::Point{2,T})

Return the polar line of a 2D point `p` with respect to the unit circle directly as a `Line{2, T}`.
"""
function polar(::Type{Line}, p::Point{2,T}) where {T}
    res = polar(p)
    return res isa Line ? res : Line(res)
end

"""
    polar_line(p::Point{2,T})

Return the polar line of a 2D point `p` with respect to the unit circle directly as a `Line{2, T}`.
"""
polar_line(p::Point{2,T}) where {T} = polar(Line, p)

export Plane, Plane2F, polar, polar_line