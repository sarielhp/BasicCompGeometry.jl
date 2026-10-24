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
@inline compare_xy(p::Point{2}, q::Point{2}; kernel::AbstractKernel = FloatingKernel()) =
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
@inline compare_squared_distance(
    p::Point{D},
    q::Point{D},
    squared::Real;
    kernel::AbstractKernel = FloatingKernel(),
) where {D} = compare_squared_distance(kernel, p, q, squared)

@inline function compare_squared_distance(
    kernel::BuiltinKernel,
    p::Point{D},
    q::Point{D},
    squared::Real,
) where {D}
    return _predicate_sign(
        _kernel_distance_sq(kernel, p, q) - _kernel_scalar(kernel, squared),
    )
end

"""
    compare_distance(p, q, radius; kernel=FloatingKernel())

Compare the distance between points with a nonnegative radius, returning -1,
0, or 1. Squaring occurs inside the selected kernel, without taking a square
root. Equality includes tangency when this predicate is used for disk tests.
"""
@inline compare_distance(
    p::Point{D},
    q::Point{D},
    radius::Real;
    kernel::AbstractKernel = FloatingKernel(),
) where {D} = compare_distance(kernel, p, q, radius)

@inline function _kernel_radius(kernel, radius)
    value = _kernel_scalar(kernel, radius)
    value >= 0 || throw(DomainError(radius, "radius must be nonnegative"))
    return value
end

@inline function compare_distance(
    kernel::BuiltinKernel,
    p::Point{D},
    q::Point{D},
    radius::Real,
) where {D}
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
@inline compare_right_arcs(
    p::Point{2},
    q::Point{2},
    height::Real,
    radius::Real;
    kernel::AbstractKernel = FloatingKernel(),
) = compare_right_arcs(kernel, p, q, height, radius)

function compare_right_arcs(
    kernel::BuiltinKernel,
    p::Point{2},
    q::Point{2},
    height::Real,
    radius::Real,
)
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
@inline compare_right_arc_to_point(
    center::Point{2},
    q::Point{2},
    radius::Real;
    kernel::AbstractKernel = FloatingKernel(),
) = compare_right_arc_to_point(kernel, center, q, radius)

function compare_right_arc_to_point(
    kernel::BuiltinKernel,
    center::Point{2},
    q::Point{2},
    radius::Real,
)
    r, y = _kernel_radius(kernel, radius), _kernel_scalar(kernel, q.y)
    a = _arc_radicand(kernel, center, y, r)
    d = _kernel_scalar(kernel, center.x) - _kernel_scalar(kernel, q.x)
    return _radical_sign(kernel, d, a, zero(a))
end

"""
    orient(a, b, c; kernel=FloatingKernel())
    orient(kernel, a, b, c)

Calculate the exact orientation determinant of the sequence of 2D points `a -> b -> c`:
`(b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)`.
Returns the signed determinant value (not normalized to sign). Under `ExactKernel()`,
this evaluates without rounding using exact rational/integer arithmetic.
"""
@inline orient(
    a::Point{2},
    b::Point{2},
    c::Point{2};
    kernel::AbstractKernel = FloatingKernel(),
) = orient(kernel, a, b, c)

@inline orient(kernel::AbstractKernel, a::Point{2}, b::Point{2}, c::Point{2}) =
    turn_sign(kernel, a, b, c)

"""
    orient_sign(a, b, c; kernel=FloatingKernel())
    orient_sign(kernel, a, b, c)

Return the normalized orientation sign of the sequence of 2D points `a -> b -> c`:
`1` for a counterclockwise (left) turn, `-1` for a clockwise (right) turn, and `0` for collinear.
Preserves exact zero under `ExactKernel()`.
"""
@inline orient_sign(
    a::Point{2},
    b::Point{2},
    c::Point{2};
    kernel::AbstractKernel = FloatingKernel(),
) = orient_sign(kernel, a, b, c)

@inline orient_sign(kernel::AbstractKernel, a::Point{2}, b::Point{2}, c::Point{2}) =
    _predicate_sign(turn_sign(kernel, a, b, c))

"""
    on_segment(p, a, b; kernel=FloatingKernel())

Determine whether point `p` lies on the closed segment `[a, b]`, including endpoints.
Collinearity and coordinate bounds are evaluated under `kernel`.
"""
@inline on_segment(
    p::Point{2},
    a::Point{2},
    b::Point{2};
    kernel::AbstractKernel = FloatingKernel(),
) = on_segment(kernel, p, a, b)

@inline function on_segment(
    kernel::AbstractKernel,
    p::Point{2},
    a::Point{2},
    b::Point{2},
)
    turn_sign(kernel, a, b, p) == 0 || return false
    kp, ka, kb = _kernel_point(kernel, p), _kernel_point(kernel, a), _kernel_point(kernel, b)
    return min(ka.x, kb.x) <= kp.x <= max(ka.x, kb.x) &&
           min(ka.y, kb.y) <= kp.y <= max(ka.y, kb.y)
end

"""
    proper_crossing(a, b, c, d; kernel=FloatingKernel())

Return `true` if and only if segment `[a, b]` and segment `[c, d]` intersect at a single point
interior to both segments (excluding endpoints and collinear overlaps).
"""
@inline proper_crossing(
    a::Point{2},
    b::Point{2},
    c::Point{2},
    d::Point{2};
    kernel::AbstractKernel = FloatingKernel(),
) = proper_crossing(kernel, a, b, c, d)

@inline function proper_crossing(
    kernel::AbstractKernel,
    a::Point{2},
    b::Point{2},
    c::Point{2},
    d::Point{2},
)
    s_abc = orient_sign(kernel, a, b, c)
    s_abd = orient_sign(kernel, a, b, d)
    s_abc * s_abd < 0 || return false
    s_cda = orient_sign(kernel, c, d, a)
    s_cdb = orient_sign(kernel, c, d, b)
    return s_cda * s_cdb < 0
end

"""
    compare_angle(a, b; kernel=FloatingKernel())
    compare_angle(kernel, a, b)

Compare two nonzero 2D direction vectors `a` and `b` in counterclockwise angular order
starting from the positive x-axis ([0, 2π)).
Returns `-1` if `a` strictly precedes `b`, `1` if `b` strictly precedes `a`, and `0` if `a` and `b`
point in the exact same direction.

# Preconditions
- Vectors `a` and `b` must both be nonzero. Throws `ArgumentError` if either is zero.
"""
@inline compare_angle(
    a::Point{2},
    b::Point{2};
    kernel::AbstractKernel = FloatingKernel(),
) = compare_angle(kernel, a, b)

@inline function compare_angle(
    kernel::AbstractKernel,
    a::Point{2},
    b::Point{2},
)
    (iszero(a.x) && iszero(a.y)) &&
        throw(ArgumentError("angular comparison requires nonzero vectors; got a = $a"))
    (iszero(b.x) && iszero(b.y)) &&
        throw(ArgumentError("angular comparison requires nonzero vectors; got b = $b"))
    ka, kb = _kernel_point(kernel, a), _kernel_point(kernel, b)
    ha = (ka.y > 0 || (ka.y == 0 && ka.x > 0)) ? 0 : 1
    hb = (kb.y > 0 || (kb.y == 0 && kb.x > 0)) ? 0 : 1
    if ha != hb
        return ha < hb ? -1 : 1
    end
    cp = turn_sign(kernel, point(0, 0), a, b)
    cs = _predicate_sign(cp)
    return cs > 0 ? -1 : (cs < 0 ? 1 : 0)
end

"""
    angle_less(a, b; kernel=FloatingKernel())
    angle_less(kernel, a, b)

Return `true` if nonzero vector `a` strictly precedes nonzero vector `b` in counterclockwise
angular order starting from the positive x-axis ([0, 2π)). Vectors pointing in the same direction
are equivalent and return `false`.

# Preconditions
- Vectors `a` and `b` must both be nonzero. Throws `ArgumentError` if either is zero.
"""
@inline angle_less(
    a::Point{2},
    b::Point{2};
    kernel::AbstractKernel = FloatingKernel(),
) = angle_less(kernel, a, b)

@inline angle_less(kernel::AbstractKernel, a::Point{2}, b::Point{2}) =
    compare_angle(kernel, a, b) < 0

"""
    id_pair(a::Integer, b::Integer)

Return the sorted pair `(min(a, b), max(a, b))` as a canonical key for an undirected edge by vertex indices.
"""
@inline id_pair(a::Integer, b::Integer) = a < b ? (a, b) : (b, a)

"""
    twice_area(vertices; kernel=FloatingKernel())

Compute twice the signed area of the polygon formed by `vertices`.
For `ExactKernel()`, returns the exact signed area without floating-point rounding.
"""
function twice_area(
    vertices;
    kernel::AbstractKernel = FloatingKernel(),
)
    n = length(vertices)
    n < 3 && return zero(_kernel_scalar(kernel, 0))
    s = orient(kernel, point(0, 0), vertices[1], vertices[2])
    for i in 2:n
        s += orient(kernel, point(0, 0), vertices[i], vertices[mod1(i + 1, n)])
    end
    return s
end

export compare_xy, compare_squared_distance, compare_distance
export orient, orient_sign, on_segment, proper_crossing
export compare_angle, angle_less, id_pair, twice_area
