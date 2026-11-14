###############################################
### Line type

"""
    Line{D, T}

Infinite line in D dimensions defined by a point `p` and a direction vector `u`.
The line is parametrised as `p + u * t` for all real `t`.
`u` represents the orientation and does not need to be normalized.
"""
struct Line{D,T}
    p::Point{D,T}
    u::Point{D,T}
end

Line(l::Line) = l
function Line(p::Point{D,T1}, u::Point{D,T2}) where {D,T1,T2}
    T = promote_type(T1, T2)
    return Line{D,T}(Point{D,T}(p), Point{D,T}(u))
end

###############################################
### Segment type

"""
    Segment{D, T}

A directed line segment in D dimensions defined by its endpoints `p` (start) and `q` (end).
"""
struct Segment{D,T} <: AbsCurve2D{T}
    p::Point{D,T}
    q::Point{D,T}
end

Base.iterate(s::Segment) = (s.p, 1)
Base.iterate(s::Segment, state) = state == 1 ? (s.q, 2) : nothing
Base.length(::Segment) = 2
Base.firstindex(::Segment) = 1
Base.lastindex(::Segment) = 2
@inline Base.getindex(s::Segment, i::Int) = i == 1 ? s.p : (i == 2 ? s.q : throw(BoundsError(s, i)))

"""
    Segment{D, T}()

Default constructor creating a zero-length segment at the origin.
"""
function Segment{D,T}() where {D,T}
    return Segment{D,T}(zeros(Point{D,T}), zeros(Point{D,T}))
end

"""
    at(seg, t)

Return the point on the segment at parameter `t`.
`t=0` returns the start point `p`, and `t=1` returns the end point `q`.
The function uses a convex combination of the endpoints.
"""
@inline at(seg::Segment{D,T}, t::Real) where {D,T} = convex_comb(seg.p, seg.q, t)

"""
    midpoint(seg::Segment)

Return the exact midpoint of segment `seg`.
"""
@inline midpoint(seg::Segment) = midpoint(seg.p, seg.q)

"""
    affine_interpolate(seg::Segment{D}, t::Real)

Compute the affine interpolation on `seg` at parameter `t`.
"""
@inline affine_interpolate(seg::Segment{D}, t::Real) where {D} = convex_comb(seg.p, seg.q, t)

function Base.show(io::IO, s::Segment{D,T}) where {D,T}
    print(io, "[[", s.p, " .. ", s.q, "]]")
end

"""
    geom_length(seg)

Return the Euclidean length of the segment.
"""
@inline geom_length(seg::Segment{D,T}) where {D,T} = dist(seg.p, seg.q)

"""
    dist(seg, point)

Return the minimum Euclidean distance between the segment `seg` and a query point `qr`.
This calculates the distance to the nearest point on the segment.
"""
@inline function dist(s::Segment{D,T}, qr::Point{D,T}) where {D,T}
    return dist(nn_point(s, qr), qr)
end

"""
    dist(seg1, seg2)

Return the minimum Euclidean distance between two segments `a` and `b` in D dimensions.
This uses an optimization approach to find the closest points on both segments.
"""
@inline function dist(a::Segment{D,T}, b::Segment{D,T}) where {D,T}
    return dist_segment_segment(a.p, a.q, b.p, b.q)
end

@inline dist(p::Point{D,T}, s::Segment{D,T}) where {D,T} = dist(s, p)

"""
    convex_coef(seg, point)

Return the parameter `t` in `[0, 1]` such that `at(seg, t)` is the nearest point on the segment to `qr`.
"""
@inline function convex_coef(s::Segment{D,T}, qr::Point{D,T}) where {D,T}
    len = geom_length(s)
    return len == 0.0 ? 0.0 : dist(s.p, qr) / len
end

"""
    nn_point(p, q, qr)

Return the nearest point on the segment defined by endpoints `p` and `q` to the query point `qr`.
"""
function nn_point(s_p::Point{D,T}, s_q::Point{D,T}, qr::Point{D,T}) where {D,T}
    R = typeof(one(T) / one(T))
    p = Point{D,R}(s_p)
    q = Point{D,R}(s_q)
    r = Point{D,R}(qr)

    delta = q - p
    d_sq = dot(delta, delta)
    iszero(d_sq) && return p

    # Projection parameter calculation: t = <qr - p, q - p> / ||q - p||^2
    t = dot(r - p, delta) / d_sq
    return convex_comb(p, q, clamp(t, zero(t), one(t)))
end

"""
    nn_point(seg, qr)

Return the nearest point on segment `seg` to query point `qr`.
"""
@inline nn_point(s::Segment{D,T}, qr::Point{D,T}) where {D,T} = nn_point(s.p, s.q, qr)

@inline function _min_endpoint_segment_dist(
    a_p::Point{D,R},
    a_q::Point{D,R},
    b_p::Point{D,R},
    b_q::Point{D,R},
) where {D,R}
    return min(
        dist_point_segment(b_p, a_p, a_q),
        dist_point_segment(b_q, a_p, a_q),
        dist_point_segment(a_p, b_p, b_q),
        dist_point_segment(a_q, b_p, b_q),
    )
end

"""
    dist_segment_segment(a_p, a_q, b_p, b_q)

Calculate the minimum distance between two segments defined by `[a_p, a_q]` and `[b_p, b_q]`.
For non-parallel segments, it solves for the parameters `s, t` that minimize the distance.
If segments are parallel, it falls back to point-segment distance calculations.
"""
function dist_segment_segment(
    a_p::Point{D,T},
    a_q::Point{D,T},
    b_p::Point{D,T},
    b_q::Point{D,T},
) where {D,T}
    R = typeof(one(T) / one(T))
    ap, aq = Point{D,R}(a_p), Point{D,R}(a_q)
    bp, bq = Point{D,R}(b_p), Point{D,R}(b_q)

    v_1 = ap - bp
    v_2 = aq - ap
    v_3 = bp - bq

    # Minimizing ||v_1 + s*v_2 + t*v_3||^2 via normal equations [A B; B C] [s; t] = [r1; r2]
    A = dot(v_2, v_2)
    B = dot(v_2, v_3)
    C = dot(v_3, v_3)
    r1 = -dot(v_1, v_2)
    r2 = -dot(v_1, v_3)
    Δ = A * C - B * B

    edge_dist = _min_endpoint_segment_dist(ap, aq, bp, bq)
    tol = R <: AbstractFloat ? max(R(1e-12), 4 * eps(R)) : zero(R)
    if Δ <= tol * (A * C)
        return edge_dist
    end

    s = clamp((r1 * C - B * r2) / Δ, zero(R), one(R))
    t = clamp((A * r2 - B * r1) / Δ, zero(R), one(R))
    d = dist(convex_comb(ap, aq, s), convex_comb(bp, bq, t))
    return min(d, edge_dist)
end

"""
    dist_point_segment(qr, p, q)

Calculate the minimum Euclidean distance from query point `qr` to the segment `[p, q]`.
"""
@inline function dist_point_segment(
    qr::Point{D,T},
    s_p::Point{D,T},
    s_q::Point{D,T},
) where {D,T}
    return dist(nn_point(s_p, s_q, qr), qr)
end

"""
    bisection_point(seg, p, q)

Determine if the segment `seg` intersects the bisector plane of points `p` and `q`.
Return `(intersects::Bool, t::Float64, point::Point)`.
`t` is the parameter on the segment where intersection occurs.
"""
function bisection_point(seg::Segment{D,T}, p::Point{D,T}, q::Point{D,T}) where {D,T}
    R = typeof(one(T) / one(T))
    sp, sq = Point{D,R}(seg.p), Point{D,R}(seg.q)
    pr, qr = Point{D,R}(p), Point{D,R}(q)

    dir = qr - pr
    mid = (qr + pr) / R(2)
    pos = dot(dir, mid)
    denom = dot(dir, sq - sp)

    is_parallel = R <: AbstractFloat ? abs(denom) < R(1e-12) : iszero(denom)
    if is_parallel
        return false, zero(R), sp
    end

    t = (pos - dot(dir, sp)) / denom
    return zero(R) <= t <= one(R), t, convex_comb(sp, sq, t)
end

"""
    segment_at_x(a::Point{2}, b::Point{2}, x::Real)
    segment_at_x(seg::Segment{2}, x::Real)

Compute the exact intersection point of the nonvertical segment `[a, b]` (or `seg`)
with the vertical line at coordinate `x`.
Returns `Point(x, y)` if `x` lies within the closed projection `[min(a.x, b.x), max(a.x, b.x)]`,
or `nothing` if `x` is strictly outside.

# Preconditions
- The segment must be nonvertical (`a.x != b.x`). Throws `ArgumentError` if `a.x == b.x`.
"""
function segment_at_x(a::Point{2}, b::Point{2}, x::Real)
    a.x == b.x && throw(ArgumentError("segment is vertical: cannot evaluate nonvertical segment at x = $x"))
    min_x, max_x = min(a.x, b.x), max(a.x, b.x)
    (x < min_x || x > max_x) && return nothing
    x == a.x && return a
    x == b.x && return b
    dx = b.x - a.x
    dy = b.y - a.y
    y = a.y + (x - a.x) * dy / dx
    return Point(x, y)
end

@inline segment_at_x(seg::Segment{2}, x::Real) = segment_at_x(seg.p, seg.q, x)

@inline on_segment(
    p::Point{2},
    seg::Segment{2};
    kernel::AbstractKernel = FloatingKernel(),
) = on_segment(kernel, p, seg.p, seg.q)

@inline on_segment(
    kernel::AbstractKernel,
    p::Point{2},
    seg::Segment{2},
) = on_segment(kernel, p, seg.p, seg.q)

@inline proper_crossing(
    seg1::Segment{2},
    seg2::Segment{2};
    kernel::AbstractKernel = FloatingKernel(),
) = proper_crossing(kernel, seg1.p, seg1.q, seg2.p, seg2.q)

@inline proper_crossing(
    kernel::AbstractKernel,
    seg1::Segment{2},
    seg2::Segment{2},
) = proper_crossing(kernel, seg1.p, seg1.q, seg2.p, seg2.q)

"""
    segment_key(a, b; kernel=FloatingKernel())
    segment_key(seg; kernel=FloatingKernel())

Return the canonical undirected segment key with endpoints ordered lexicographically
(first by x-coordinate, then by y-coordinate).
"""
@inline segment_key(
    a::Point{2},
    b::Point{2};
    kernel::AbstractKernel = FloatingKernel(),
) = segment_key(kernel, a, b)

@inline segment_key(
    seg::Segment{2};
    kernel::AbstractKernel = FloatingKernel(),
) = segment_key(kernel, seg.p, seg.q)

@inline segment_key(kernel::AbstractKernel, seg::Segment{2}) =
    segment_key(kernel, seg.p, seg.q)

@inline function segment_key(
    kernel::AbstractKernel,
    a::Point{2},
    b::Point{2},
)
    return compare_xy(kernel, a, b) <= 0 ? Segment(a, b) : Segment(b, a)
end

@inline function segment_key(
    seg::Tuple{Point{2},Point{2}};
    kernel::AbstractKernel = FloatingKernel(),
)
    return compare_xy(kernel, seg[1], seg[2]) <= 0 ? (seg[1], seg[2]) : (seg[2], seg[1])
end

"""
    polygon_segments(vertices)

Return the sequence of directed segments `[Segment(v[i], v[i+1])]` connecting successive vertices in cyclic order.
"""
function polygon_segments(vertices)
    n = length(vertices)
    return [Segment(vertices[i], vertices[mod1(i + 1, n)]) for i in 1:n]
end

"""
    Segment2F

A 2-dimensional directed segment with `Float64` coordinates.
"""
const Segment2F = Segment{2,Float64}

export Line, Segment, Segment2F
export nn_point, geom_length, midpoint, affine_interpolate, segment_at_x
export on_segment, proper_crossing, segment_key, polygon_segments
