# Cubic Bézier Curves and Splines

"""
    CubicBezier{D, T}

A cubic Bézier curve in `D` dimensions defined by four control points:
`p0`, `p1`, `p2`, and `p3` of type `Point{D, T}`.
"""
struct CubicBezier{D,T}
    p0::Point{D,T}
    p1::Point{D,T}
    p2::Point{D,T}
    p3::Point{D,T}
end

const CubicBezier2F = CubicBezier{2,Float64}
const CubicBezier3F = CubicBezier{3,Float64}

function CubicBezier(
    p0::Point{D,T},
    p1::Point{D,S},
    p2::Point{D,U},
    p3::Point{D,V},
) where {D,T,S,U,V}
    R = promote_type(T, S, U, V)
    return CubicBezier{D,R}(Point{D,R}(p0), Point{D,R}(p1), Point{D,R}(p2), Point{D,R}(p3))
end

function (b::CubicBezier{D,T})(t::Real) where {D,T}
    u = one(T) - T(t)
    tt = T(t)
    # de Casteljau evaluation
    q0 = u * b.p0 + tt * b.p1
    q1 = u * b.p1 + tt * b.p2
    q2 = u * b.p2 + tt * b.p3
    r0 = u * q0 + tt * q1
    r1 = u * q1 + tt * q2
    return u * r0 + tt * r1
end

"""
    derivative(b::CubicBezier, t::Real)

Compute the tangent vector B'(t) at parameter `t`.
"""
function derivative(b::CubicBezier{D,T}, t::Real) where {D,T}
    u = one(T) - T(t)
    tt = T(t)
    return 3 * u * u * (b.p1 - b.p0) +
           6 * u * tt * (b.p2 - b.p1) +
           3 * tt * tt * (b.p3 - b.p2)
end

import Base: split

"""
    split(b::CubicBezier, t::Real = 0.5)
    subdivide(b::CubicBezier, t::Real = 0.5)

Subdivide a cubic Bézier curve at parameter `t` using de Casteljau's algorithm.
Returns a tuple `(left, right)` of two `CubicBezier` curves.
"""
function Base.split(b::CubicBezier{D,T}, t::Real = 0.5) where {D,T}
    u = one(T) - T(t)
    tt = T(t)
    q0 = u * b.p0 + tt * b.p1
    q1 = u * b.p1 + tt * b.p2
    q2 = u * b.p2 + tt * b.p3
    r0 = u * q0 + tt * q1
    r1 = u * q1 + tt * q2
    mid = u * r0 + tt * r1
    left = CubicBezier{D,T}(b.p0, q0, r0, mid)
    right = CubicBezier{D,T}(mid, r1, q2, b.p3)
    return (left, right)
end

subdivide(b::CubicBezier, t::Real = 0.5) = Base.split(b, t)

@inline function _eval_1d_bezier(v0::T, v1::T, v2::T, v3::T, t::T) where {T}
    u = one(T) - t
    return u * (u * (u * v0 + t * v1) + t * (u * v1 + t * v2)) +
           t * (u * (u * v1 + t * v2) + t * (u * v2 + t * v3))
end

@inline function _update_coord_extrema(
    lo::T,
    hi::T,
    v0::T,
    v1::T,
    v2::T,
    v3::T,
    t::T,
) where {T}
    if zero(T) < t < one(T)
        val = _eval_1d_bezier(v0, v1, v2, v3, t)
        return min(lo, val), max(hi, val)
    end
    return lo, hi
end

function _coord_minmax(v0::T, v1::T, v2::T, v3::T) where {T}
    lo, hi = minmax(v0, v3)
    # Derivative is A*t^2 + B*t + C = 0
    A = 3 * (-v0 + 3 * v1 - 3 * v2 + v3)
    B = 6 * (v0 - 2 * v1 + v2)
    C = 3 * (v1 - v0)
    if abs(A) < 1e-12
        abs(B) > 1e-12 && return _update_coord_extrema(lo, hi, v0, v1, v2, v3, -C / B)
        return lo, hi
    end
    disc = B * B - 4 * A * C
    disc < zero(T) && return lo, hi
    s = sqrt(disc)
    lo, hi = _update_coord_extrema(lo, hi, v0, v1, v2, v3, (-B - s) / (2 * A))
    return _update_coord_extrema(lo, hi, v0, v1, v2, v3, (-B + s) / (2 * A))
end

@inline function _bezier_extrema(b::CubicBezier{D,T}) where {D,T}
    extrema_tup = ntuple(k -> _coord_minmax(b.p0[k], b.p1[k], b.p2[k], b.p3[k]), Val(D))
    mini = Point{D,T}(ntuple(k -> extrema_tup[k][1], Val(D)))
    maxi = Point{D,T}(ntuple(k -> extrema_tup[k][2], Val(D)))
    return mini, maxi
end

"""
    bbox(b::CubicBezier{D, T})

Compute exact axis-aligned bounding box of a cubic Bézier curve.
"""
function bbox(b::CubicBezier{D,T}) where {D,T}
    mini, maxi = _bezier_extrema(b)
    return BBox(mini, maxi)
end

"""
    flatten(b::CubicBezier; tol=1e-3, max_depth=10)

Adaptively flatten a cubic Bézier curve into a sequence of points within distance tolerance `tol`.
"""
function flatten(b::CubicBezier{D,T}; tol::Real = 1e-3, max_depth::Int = 10) where {D,T}
    pts = Point{D,T}[b.p0]
    _flatten_rec!(pts, b, T(tol), max_depth)
    return PntSeq(pts)
end

function _flatten_rec!(
    pts::Vector{Point{D,T}},
    b::CubicBezier{D,T},
    tol::T,
    depth::Int,
) where {D,T}
    seg = Segment(b.p0, b.p3)
    d1 = dist(b.p1, seg)
    d2 = dist(b.p2, seg)
    if (d1 <= tol && d2 <= tol) || depth <= 0
        push!(pts, b.p3)
    else
        left, right = split(b, 0.5)
        _flatten_rec!(pts, left, tol, depth - 1)
        _flatten_rec!(pts, right, tol, depth - 1)
    end
end

Base.show(io::IO, b::CubicBezier{D,T}) where {D,T} =
    print(io, "CubicBezier(", b.p0, ", ", b.p1, ", ", b.p2, ", ", b.p3, ")")


# -----------------------------------------------------------------------------
# Cubic Spline
# -----------------------------------------------------------------------------

"""
    CubicSpline{D, T}

A composite curve consisting of a chain of `CubicBezier{D, T}` segments.
"""
struct CubicSpline{D,T}
    segments::Vector{CubicBezier{D,T}}
    is_closed::Bool
end

CubicSpline(segs::Vector{CubicBezier{D,T}}; closed::Bool = false) where {D,T} =
    CubicSpline{D,T}(segs, closed)

Base.length(s::CubicSpline) = length(s.segments)

function (s::CubicSpline{D,T})(t::Real) where {D,T}
    n = length(s.segments)
    n == 0 && error("Empty spline")
    tc = clamp(T(t), zero(T), T(n))
    idx = min(floor(Int, tc) + 1, n)
    local_t = tc - (idx - 1)
    return s.segments[idx](local_t)
end

function bbox(s::CubicSpline{D,T}) where {D,T}
    isempty(s.segments) && error("Empty spline has no bbox")
    bb = bbox(s.segments[1])
    for i = 2:length(s.segments)
        mini, maxi = _bezier_extrema(s.segments[i])
        bound!(bb, mini)
        bound!(bb, maxi)
    end
    return bb
end

function flatten(s::CubicSpline{D,T}; tol::Real = 1e-3, max_depth::Int = 10) where {D,T}
    isempty(s.segments) && return PntSeq{D,T}(Point{D,T}[])
    all_pts = Point{D,T}[s.segments[1].p0]
    sizehint!(all_pts, length(s.segments) + 1)
    tol_T = T(tol)
    for seg in s.segments
        _flatten_rec!(all_pts, seg, tol_T, max_depth)
    end
    return PntSeq(all_pts)
end

"""
    interpolate_catmull_rom(pts::AbsPntSeq{D, T}; tension=0.5, closed=false)

Construct a smooth C^1 cubic spline interpolating the given points via Catmull-Rom.
"""
function interpolate_catmull_rom(
    pts::AbsPntSeq{D,T};
    tension::Real = 0.5,
    closed::Bool = false,
) where {D,T}
    n = cardin(pts)
    n >= 2 || error("Catmull-Rom interpolation requires at least 2 points")
    segs = CubicBezier{D,T}[]
    scale = (one(T) - T(tension)) / 3

    if closed
        _catmull_rom_closed!(segs, pts, scale)
    else
        _catmull_rom_open!(segs, pts, scale)
    end
    return CubicSpline(segs; closed = closed)
end

# Appends the segments of the closed Catmull-Rom spline through `pts`.
function _catmull_rom_closed!(segs, pts, scale)
    n = cardin(pts)
    for i = 1:n
        im1 = i == 1 ? n : i - 1
        ip1 = i == n ? 1 : i + 1
        ip2 = ip1 == n ? 1 : ip1 + 1
        p_prev = pts[im1]
        p_curr = pts[i]
        p_next = pts[ip1]
        p_next2 = pts[ip2]
        v_curr = scale * (p_next - p_prev)
        v_next = scale * (p_next2 - p_curr)
        push!(segs, CubicBezier(p_curr, p_curr + v_curr, p_next - v_next, p_next))
    end
    return segs
end

# Appends the segments of the open Catmull-Rom spline through `pts`; the end
# tangents are one-sided.
function _catmull_rom_open!(segs, pts, scale)
    n = cardin(pts)
    for i = 1:(n-1)
        p_curr = pts[i]
        p_next = pts[i+1]
        v_curr = i == 1 ? scale * (p_next - p_curr) * 2 : scale * (p_next - pts[i-1])
        v_next = i + 1 == n ? scale * (p_next - p_curr) * 2 : scale * (pts[i+2] - p_curr)
        push!(segs, CubicBezier(p_curr, p_curr + v_curr, p_next - v_next, p_next))
    end
    return segs
end

function _solve_natural_tridiag!(
    p1k::Vector{T},
    b::Vector{T},
    r::Vector{T},
    pts::AbsPntSeq{D,T},
    k::Int,
    m::Int,
) where {D,T}
    fill!(b, 4 * one(T))
    b[1] = 2 * one(T)
    b[m] = 7 * one(T)
    r[1] = pts[1][k] + 2 * pts[2][k]
    for i = 2:(m-1)
        r[i] = 4 * pts[i][k] + 2 * pts[i+1][k]
    end
    r[m] = 8 * pts[m][k] + pts[m+1][k]

    # Thomas forward elimination (superdiagonal is identically 1)
    for i = 2:m
        a_i = i == m ? 2 * one(T) : one(T)
        w = a_i / b[i-1]
        b[i] -= w
        r[i] -= w * r[i-1]
    end
    p1k[m] = r[m] / b[m]
    for i = (m-1):-1:1
        p1k[i] = (r[i] - p1k[i+1]) / b[i]
    end
    return p1k
end

"""
    interpolate_natural_spline(pts::AbsPntSeq{D, T}; closed=false)

Construct a smooth C^2 natural cubic spline interpolating the given points.
"""
function interpolate_natural_spline(pts::AbsPntSeq{D,T}; closed::Bool = false) where {D,T}
    n = cardin(pts)
    n >= 3 || return interpolate_catmull_rom(pts; closed = closed)
    m = n - 1
    P1_coords = ntuple(_ -> Vector{T}(undef, m), Val(D))
    b = Vector{T}(undef, m)
    r = Vector{T}(undef, m)
    for k = 1:D
        _solve_natural_tridiag!(P1_coords[k], b, r, pts, k, m)
    end

    segs = Vector{CubicBezier{D,T}}(undef, m)
    for i = 1:m
        p0 = pts[i]
        p3 = pts[i+1]
        p1 = Point{D,T}(ntuple(k -> P1_coords[k][i], Val(D)))
        p2 = if i < m
            Point{D,T}(ntuple(k -> 2 * p3[k] - P1_coords[k][i+1], Val(D)))
        else
            Point{D,T}(ntuple(k -> (p3[k] + P1_coords[k][m]) / 2, Val(D)))
        end
        segs[i] = CubicBezier{D,T}(p0, p1, p2, p3)
    end
    return CubicSpline(segs; closed = false)
end

BBox(b::CubicBezier) = bbox(b)
BBox(s::CubicSpline) = bbox(s)

export CubicBezier, CubicBezier2F, CubicBezier3F, CubicSpline
export subdivide, flatten, interpolate_catmull_rom, interpolate_natural_spline, bbox
