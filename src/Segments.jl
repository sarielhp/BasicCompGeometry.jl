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
    d_sq = dist_sq(s_p, s_q)
    if d_sq == 0.0
        return s_p
    end
    # Projection parameter calculation: t = <qr - p, q - p> / ||q - p||^2
    t = dot(qr - s_p, s_q - s_p) / d_sq
    return convex_comb(s_p, s_q, clamp(t, 0.0, 1.0))
end

"""
    nn_point(seg, qr)

Return the nearest point on segment `seg` to query point `qr`.
"""
@inline nn_point(s::Segment{D,T}, qr::Point{D,T}) where {D,T} = nn_point(s.p, s.q, qr)

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
    v_1 = a_p - b_p
    v_2 = a_q - a_p
    v_3 = b_p - b_q

    # Minimizing ||v_1 + s*v_2 + t*v_3||^2
    # Normal equations lead to a 2x2 system
    c = [-dot(v_1, v_2), -dot(v_1, v_3)]
    m = [
        dot(v_2, v_2) dot(v_2, v_3);
        dot(v_2, v_3) dot(v_3, v_3)
    ]

    # Degenerate case check (parallel segments)
    if det(m) < 1e-12 * norm(m)
        return min(
            dist_point_segment(b_p, a_p, a_q),
            dist_point_segment(b_q, a_p, a_q),
            dist_point_segment(a_p, b_p, b_q),
            dist_point_segment(a_q, b_p, b_q),
        )
    end

    sol = m \ c
    s = clamp(sol[1], 0.0, 1.0)
    t = clamp(sol[2], 0.0, 1.0)

    d = dist(convex_comb(a_p, a_q, s), convex_comb(b_p, b_q, t))

    # Boundary check for global minimum
    return min(
        d,
        dist_point_segment(b_p, a_p, a_q),
        dist_point_segment(b_q, a_p, a_q),
        dist_point_segment(a_p, b_p, b_q),
        dist_point_segment(a_q, b_p, b_q),
    )
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
    dir = q - p
    mid = (q + p) / 2
    pos = dot(dir, mid)

    seg_vec = seg.q - seg.p
    denom = dot(dir, seg_vec)

    if abs(denom) < 1e-12 # Segment parallel to bisector plane
        return false, 0.0, seg.p
    end

    t = (pos - dot(dir, seg.p)) / denom
    on_seg = (0.0 <= t <= 1.0)
    return on_seg, t, at(seg, t)
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
