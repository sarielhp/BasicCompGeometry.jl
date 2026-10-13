# Kernel-first methods are extension points; keyword wrappers are the user API.
@inline _kernel_point(kernel, p::Point) = map(x -> _kernel_scalar(kernel, x), p)

@inline function turn_sign(kernel::BuiltinKernel, p::Point{2}, q::Point{2}, r::Point{2})
    a, b, c = _kernel_point(kernel, p), _kernel_point(kernel, q), _kernel_point(kernel, r)
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
end

"""
    compare_xy(p, q; kernel=FloatingKernel())

Compare two 2D points lexicographically, returning -1, 0, or 1.
"""
@inline compare_xy(p::Point{2}, q::Point{2}; kernel::AbstractKernel=FloatingKernel()) =
    compare_xy(kernel, p, q)

@inline function compare_xy(kernel::BuiltinKernel, p::Point{2}, q::Point{2})
    a, b = _kernel_point(kernel, p), _kernel_point(kernel, q)
    x_order = _predicate_sign(a.x - b.x)
    return x_order == 0 ? _predicate_sign(a.y - b.y) : x_order
end

@inline function _kernel_distance_sq(kernel, p::Point{D}, q::Point{D}) where {D}
    return sum(abs2, _kernel_point(kernel, p) - _kernel_point(kernel, q))
end

"""
    compare_squared_distance(p, q, squared; kernel=FloatingKernel())

Compare the squared distance between points with `squared`, returning -1, 0,
or 1. Supply an exact squared threshold when using `ExactKernel`; the kernel
cannot undo rounding or overflow in a threshold computed by the caller.
"""
@inline compare_squared_distance(p::Point{D}, q::Point{D}, squared::Real;
    kernel::AbstractKernel=FloatingKernel()) where {D} =
    compare_squared_distance(kernel, p, q, squared)

@inline function compare_squared_distance(kernel::BuiltinKernel, p::Point{D},
    q::Point{D}, squared::Real) where {D}
    return _predicate_sign(_kernel_distance_sq(kernel, p, q) - _kernel_scalar(kernel, squared))
end

"""
    compare_distance(p, q, radius; kernel=FloatingKernel())

Compare the distance between points with a nonnegative radius, returning -1,
0, or 1. Squaring occurs inside the selected kernel, without taking a square
root. Equality includes tangency when this predicate is used for disk tests.
"""
@inline compare_distance(p::Point{D}, q::Point{D}, radius::Real;
    kernel::AbstractKernel=FloatingKernel()) where {D} =
    compare_distance(kernel, p, q, radius)

@inline function _kernel_radius(kernel, radius)
    value = _kernel_scalar(kernel, radius)
    value >= 0 || throw(DomainError(radius, "radius must be nonnegative"))
    return value
end

@inline function compare_distance(kernel::BuiltinKernel, p::Point{D}, q::Point{D},
    radius::Real) where {D}
    r = _kernel_radius(kernel, radius)
    return _predicate_sign(_kernel_distance_sq(kernel, p, q) - r * r)
end

# Sign of d + sqrt(a) - sqrt(b), with nonnegative radicands.
# For d > 0, isolation gives 2d sqrt(a) versus b-a-d^2. Check the
# latter's sign before squaring so equality and tangency remain exact.
function _radical_sign(kernel::ExactKernel, d, a, b)
    d < 0 && return -_radical_sign(kernel, -d, b, a)
    d == 0 && return _predicate_sign(a - b)
    t = b - a - d * d
    t < 0 && return 1
    return _predicate_sign(4 * d * d * a - t * t)
end

@inline _radical_sign(::FloatingKernel, d, a, b) = _predicate_sign(d + sqrt(a) - sqrt(b))

@inline function _arc_radicand(kernel, center::Point{2}, height, radius)
    offset = height - _kernel_scalar(kernel, center.y)
    value = radius * radius - offset * offset
    value >= 0 || throw(DomainError(height, "height lies outside the circle's arc domain"))
    return value
end

"""
    compare_right_arcs(p, q, height, radius; kernel=FloatingKernel())

Compare the right circular boundaries `p.x + sqrt(radius^2-(height-p.y)^2)`
and the corresponding expression for `q`, returning -1, 0, or 1. Both circles
have the same nonnegative radius. The height must be in both arc domains.
The exact kernel compares radicals algebraically, without constructing them.
"""
@inline compare_right_arcs(p::Point{2}, q::Point{2}, height::Real, radius::Real;
    kernel::AbstractKernel=FloatingKernel()) = compare_right_arcs(kernel, p, q, height, radius)

function compare_right_arcs(kernel::BuiltinKernel, p::Point{2}, q::Point{2},
    height::Real, radius::Real)
    r, y = _kernel_radius(kernel, radius), _kernel_scalar(kernel, height)
    a, b = _arc_radicand(kernel, p, y, r), _arc_radicand(kernel, q, y, r)
    d = _kernel_scalar(kernel, p.x) - _kernel_scalar(kernel, q.x)
    return _radical_sign(kernel, d, a, b)
end

"""
    compare_right_arc_to_point(center, q, radius; kernel=FloatingKernel())

Compare the right circular boundary at `q.y` with `q.x`, returning -1, 0, or 1.
The height must be in the circle's arc domain; the radius must be nonnegative.
"""
@inline compare_right_arc_to_point(center::Point{2}, q::Point{2}, radius::Real;
    kernel::AbstractKernel=FloatingKernel()) = compare_right_arc_to_point(kernel, center, q, radius)

function compare_right_arc_to_point(kernel::BuiltinKernel, center::Point{2},
    q::Point{2}, radius::Real)
    r, y = _kernel_radius(kernel, radius), _kernel_scalar(kernel, q.y)
    a = _arc_radicand(kernel, center, y, r)
    d = _kernel_scalar(kernel, center.x) - _kernel_scalar(kernel, q.x)
    return _radical_sign(kernel, d, a, zero(a))
end

export compare_xy, compare_squared_distance, compare_distance
export compare_right_arcs, compare_right_arc_to_point
