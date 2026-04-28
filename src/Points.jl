###################################################################
### Point type

"""
    Point{D, T}

Type alias for `SVector{D, T}` from `StaticArrays`. 
Provides a concise way to work with fixed-size coordinate vectors.
It inherits all efficient operations from StaticArrays.
"""
const Point = SVector

"""
    Point2I

A 2-dimensional point with `Int64` coordinates.
"""
const Point2I = Point{2,Int64}

"""
    Point2F

A 2-dimensional point with `Float64` coordinates.
"""
const Point2F = Point{2,Float64}

"""
    Point3F

A 3-dimensional point with `Float64` coordinates.
"""
const Point3F = Point{3,Float64}

"""
    dist_sq(p, q)

Calculate the square of the Euclidean distance between two points `p` and `q`.
This is more efficient than `dist(p, q)` if only relative distances are needed.
"""
@inline dist_sq(p::Point{D,T}, q::Point{D,T}) where {D,T} = sum(abs2, p - q)

"""
    dist(p, q)

Calculate the Euclidean distance between two points `p` and `q`.
"""
@inline dist(p::Point{D,T}, q::Point{D,T}) where {D,T} = sqrt(dist_sq(p, q))

"""
    convex_comb(p, q, t)

Return the convex combination of two points: `(1-t)p + tq`.
When `t=0`, returns `p`. When `t=1`, returns `q`.
"""
@inline function convex_comb(p::Point{D,T}, q::Point{D,T}, t::Real) where {D,T}
    return p + (q - p) * t
end

"""
    turn_sign(p, q, r)

Return the sign of the turn formed by the sequence of 2D points `p -> q -> r`.
A positive value indicates a counter-clockwise (left) turn, 
a negative value indicates a clockwise (right) turn, 
and zero indicates the points are collinear.
"""
@inline function turn_sign(p::Point{2,T}, q::Point{2,T}, r::Point{2,T}) where {T}
    return (q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1])
end

"""
    is_left_turn(p, q, r)

Return `true` if the sequence of 2D points `p -> q -> r` performs a counter-clockwise (left) turn.
"""
@inline is_left_turn(p::Point{2,T}, q::Point{2,T}, r::Point{2,T}) where {T} =
    turn_sign(p, q, r) > zero(T)

"""
    is_right_turn(p, q, r)

Return `true` if the sequence of 2D points `p -> q -> r` performs a clockwise (right) turn.
"""
@inline is_right_turn(p::Point{2,T}, q::Point{2,T}, r::Point{2,T}) where {T} =
    turn_sign(p, q, r) < zero(T)

"""
    is_left_eq_turn(p, q, r)

Return `true` if the sequence of 2D points `p -> q -> r` performs a counter-clockwise (left) turn 
or if the points are collinear.
"""
@inline is_left_eq_turn(p::Point{2,T}, q::Point{2,T}, r::Point{2,T}) where {T} =
    turn_sign(p, q, r) >= zero(T)

"""
    is_right_eq_turn(p, q, r)

Return `true` if the sequence of 2D points `p -> q -> r` performs a clockwise (right) turn 
or if the points are collinear.
"""
@inline is_right_eq_turn(p::Point{2,T}, q::Point{2,T}, r::Point{2,T}) where {T} =
    turn_sign(p, q, r) <= zero(T)

"""
    is_collinear(p, q, r)

Return `true` if the three 2D points `p`, `q`, and `r` are collinear.
"""
@inline is_collinear(p::Point{2,T}, q::Point{2,T}, r::Point{2,T}) where {T} =
    turn_sign(p, q, r) == zero(T)

"""
    point(args...)

Construct a `Point` (SVector) from a list of coordinates or a single collection.
Example: `point(1.0, 2.0)` or `point([1, 2, 3])`.
"""
@inline point(args...) = Point(args...)

"""
    rand_gaussian(D, T=Float64)

Return a random `Point{D, T}` where each coordinate is sampled from a standard normal distribution.
"""
function rand_gaussian(D::Int, T::Type{<:Real} = Float64)
    return Point{D,T}(randn(T, D))
end

"""
    rand_point(D, T=Float64)

Return a random `Point{D, T}` where each coordinate is sampled uniformly from the interval `[0, 1)`.
"""
rand_point(D::Int, T::Type{<:Real} = Float64) = rand(Point{D,T})

"""
    max(p, q)

Return a point where each coordinate is the maximum of the corresponding coordinates of `p` and `q`.
"""
@inline Base.max(p::Point{D,T}, q::Point{D,T}) where {D,T} = Point{D,T}(max.(p, q))

"""
    min(p, q)

Return a point where each coordinate is the minimum of the corresponding coordinates of `p` and `q`.
"""
@inline Base.min(p::Point{D,T}, q::Point{D,T}) where {D,T} = Point{D,T}(min.(p, q))

export Point, Point2F, Point2I, Point3F
export dist, dist_sq
export rand_point, rand_gaussian
export convex_comb
export turn_sign,
    is_left_turn, is_right_turn, is_left_eq_turn, is_right_eq_turn, is_collinear
export point
