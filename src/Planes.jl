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

Return the polar plane of point `p` with respect to the unit sphere.
The polar plane has normal `p` and passes through `p / dot(p, p)`.
"""
function polar(p::Point{D,T}) where {D,T}
    d = dot(p, p)
    return Plane{D,T}(p / d, p)
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
In 2D, a plane and a line are the same (d-1 dimensional subspace),
so this simply maps the line's direction to the plane's normal.
"""
Plane(l::Line{2,T}) where {T} = Plane{2,T}(l.p, l.u)

"""
    Line(pl::Plane{2,T})

Construct a 2D `Line` from a `Plane`.
In 2D, a line and a plane are the same (d-1 dimensional subspace),
so this simply maps the plane's normal to the line's direction.
"""
Line(pl::Plane{2,T}) where {T} = Line{2,T}(pl.p, pl.n)

"""
    polar(l)

Return the polar point of the 2D line `l` with respect to the unit sphere.
This is the same as `polar(::Plane)` since in 2D a line is a (d-1)-dimensional subspace.
"""
function polar(l::Line{2,T}) where {T}
    return polar(Plane(l))
end

export Plane, Plane2F, polar