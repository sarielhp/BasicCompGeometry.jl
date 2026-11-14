# Minkowski Sum of 2D Polygons and Shapes

@inline _minkowski_coord_type(::ExactKernel, ::Type{T}) where {T} =
    T <: Rational ? T : Rational{BigInt}
@inline _minkowski_coord_type(::AbstractKernel, ::Type{T}) where {T} = T
@inline _minkowski_coord_type(
    kernel::AbstractKernel,
    ::Type{T1},
    ::Type{T2},
) where {T1,T2} = _minkowski_coord_type(kernel, promote_type(T1, T2))

@inline _pts_vec(s::Segment{2,T}) where {T} = Point{2,T}[s.p, s.q]
@inline _pts_vec(p::Point{2,T}) where {T} = Point{2,T}[p]

struct _TrianglePieces{R} <: AbstractVector{SVector{3,Point{2,R}}}
    pts::Vector{Point{2,R}}
    tris::Vector{NTuple{3,Int}}
end

Base.size(tp::_TrianglePieces) = (length(tp.tris),)
@inline function Base.getindex(tp::_TrianglePieces{R}, k::Int) where {R}
    i, j, m = tp.tris[k]
    return SVector{3,Point{2,R}}(tp.pts[i], tp.pts[j], tp.pts[m])
end

function _has_all_ccw_turns(
    pts::Vector{Point{2,R}},
    kernel::AbstractKernel,
    tol::Float64,
) where {R}
    k = length(pts)
    for i in 1:k
        p = pts[i == 1 ? k : i - 1]
        q = pts[i]
        r = pts[i == k ? 1 : i + 1]
        ts = orient_sign(kernel, p, q, r)
        if kernel isa ExactKernel
            ts < 0 && return false
        elseif ts < 0
            d1 = dist(p, q)
            d2 = dist(q, r)
            cross = (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x)
            cross < -tol * max(d1 * d2, 1.0) && return false
        end
    end
    return true
end

function _lowest_left_index(pts::AbstractVector{Point{2,R}}) where {R}
    idx_min = 1
    for i in 2:length(pts)
        p = pts[i]
        m = pts[idx_min]
        if p.y < m.y || (p.y == m.y && p.x < m.x)
            idx_min = i
        end
    end
    return idx_min
end

function _edges_angle_monotone(
    pts::Vector{Point{2,R}},
    idx_min::Int,
    kernel::AbstractKernel,
) where {R}
    k = length(pts)
    prev_edge = Point{2,R}(zero(R), zero(R))
    has_prev = false
    nz_count = 0
    idx = idx_min
    for _ in 1:k
        nxt = idx == k ? 1 : idx + 1
        e = pts[nxt] - pts[idx]
        idx = nxt
        (iszero(e.x) && iszero(e.y)) && continue
        nz_count += 1
        if has_prev && compare_angle(kernel, prev_edge, e) > 0
            return false
        end
        prev_edge = e
        has_prev = true
    end
    return nz_count >= 3
end

"""
    is_convex(poly; kernel=FloatingKernel())

Determine whether a 2D polygon `poly` is convex.

Returns `true` if the polygon's vertices form a valid non-degenerate convex polygon (all consecutive
turns are either all left turns or all right turns, and the boundary makes exactly one complete
rotation around the polygon without self-intersection). Collinear intermediate vertices along edges
are permitted. Single points and segments return `true`.

# Arguments
- `poly`: vertices of the polygon, supplied as an `AbstractVector{<:Point{2}}` or [`AbsPntSeq{2}`](@ref).
- `kernel`: computational geometry kernel ([`FloatingKernel()`](@ref) or [`ExactKernel()`](@ref)).
"""
function is_convex(
    poly::AbstractVector{<:Point{2}};
    kernel::AbstractKernel = FloatingKernel(),
)
    length(poly) <= 2 && return true
    R = _minkowski_coord_type(kernel, _scalar_type(poly))
    return _is_convex_impl(poly, kernel, R)
end

function _is_convex_impl(poly, kernel::AbstractKernel, ::Type{R}) where {R}
    deduped = Point{2,R}[]
    sizehint!(deduped, length(poly))
    _collect_deduped_points!(deduped, poly, kernel, 1e-11)

    length(deduped) <= 2 && return true
    area2 = twice_area(deduped; kernel = kernel)
    iszero(area2) && return false
    area2 < 0 && reverse!(deduped)

    _has_all_ccw_turns(deduped, kernel, 1e-11) || return false
    idx_min = _lowest_left_index(deduped)
    return _edges_angle_monotone(deduped, idx_min, kernel)
end

is_convex(poly::AbsPntSeq{2}; kernel::AbstractKernel = FloatingKernel()) =
    is_convex(Points(poly); kernel = kernel)

is_convex(::Point{2}; kernel::AbstractKernel = FloatingKernel()) = true
is_convex(::Segment{2}; kernel::AbstractKernel = FloatingKernel()) = true

is_convex(poly; kernel::AbstractKernel = FloatingKernel()) =
    is_convex(collect(poly); kernel = kernel)

"""
    minkowski_sum(poly1, poly2; kernel=FloatingKernel(), hull=false)

Compute the Minkowski sum ``P \\oplus Q = \\{ p + q \\mid p \\in P, q \\in Q \\}`` of two 2D polygons or geometric shapes.

For two convex polygons with ``n`` and ``m`` vertices, the Minkowski sum is computed in optimal ``O(n + m)``
time by merging the cyclically ordered edge lists by polar angle without quadratic pairwise expansion.
The vertices of the resulting convex polygon are returned in counterclockwise order.

If `hull=true`, computes the convex hull of the Minkowski sum ``\\mathrm{CH}(P \\oplus Q) = \\mathrm{CH}(P) \\oplus \\mathrm{CH}(Q)``
in ``O(n \\log n + m \\log m)`` time, valid for any input polygons (including non-convex polygons).
If `hull=false` (default), both polygons must be convex; otherwise an `ArgumentError` is thrown.
For non-convex polygons without taking the convex hull, see [`minkowski_sum_decomposition`](@ref).

# Arguments
- `poly1`, `poly2`: convex polygons (as `AbstractVector{<:Point{2}}` or [`AbsPntSeq{2}`](@ref)),
  segments ([`Segment{2}`](@ref)), or points ([`Point{2}`](@ref)).
- `kernel`: computational geometry kernel ([`FloatingKernel()`](@ref) or [`ExactKernel()`](@ref)).
- `hull`: whether to compute the convex hull of the Minkowski sum (default `false`).

# Returns
- `Vector{Point{2, R}}`: vertices of the Minkowski sum polygon in counterclockwise order.
  If both inputs are single points, returns a [`Point{2}`](@ref).
  If one input is a point and the other is a segment, returns a [`Segment{2}`](@ref).

# Examples
```julia
sq1 = [point(0.0, 0.0), point(1.0, 0.0), point(1.0, 1.0), point(0.0, 1.0)]
sq2 = [point(0.0, 0.0), point(2.0, 0.0), point(2.0, 2.0), point(0.0, 2.0)]
minkowski_sum(sq1, sq2) # 4 vertices from [0, 0] to [3, 3]
```
"""
function minkowski_sum(
    poly1::AbstractVector{<:Point{2}},
    poly2::AbstractVector{<:Point{2}};
    kernel::AbstractKernel = FloatingKernel(),
    hull::Bool = false,
)
    if hull
        h1 = convex_hull(poly1; kernel = kernel)
        h2 = convex_hull(poly2; kernel = kernel)
        return _minkowski_convex(Points(h1), Points(h2), kernel)
    end

    if !is_convex(poly1; kernel = kernel)
        throw(
            ArgumentError(
                "minkowski_sum: poly1 is not convex. Use minkowski_sum_decomposition for non-convex polygons, or set hull=true to compute the convex hull of the sum.",
            ),
        )
    end
    if !is_convex(poly2; kernel = kernel)
        throw(
            ArgumentError(
                "minkowski_sum: poly2 is not convex. Use minkowski_sum_decomposition for non-convex polygons, or set hull=true to compute the convex hull of the sum.",
            ),
        )
    end

    return _minkowski_convex(poly1, poly2, kernel)
end

function _minkowski_convex(poly1, poly2, kernel::AbstractKernel)
    R = _minkowski_coord_type(kernel, _scalar_type(poly1), _scalar_type(poly2))
    return _minkowski_convex_impl(poly1, poly2, kernel, R)
end

function _minkowski_convex_impl(poly1, poly2, kernel::AbstractKernel, ::Type{R}) where {R}
    p1 = Point{2,R}[]
    sizehint!(p1, length(poly1))
    _collect_deduped_points!(p1, poly1, kernel, 1e-11)

    p2 = Point{2,R}[]
    sizehint!(p2, length(poly2))
    _collect_deduped_points!(p2, poly2, kernel, 1e-11)

    n, m = length(p1), length(p2)
    (n == 0 || m == 0) && return Point{2,R}[]
    n == 1 && m == 1 && return Point{2,R}[p1[1]+p2[1]]
    n == 1 && return Point{2,R}[q + p1[1] for q in p2]
    m == 1 && return Point{2,R}[p + p2[1] for p in p1]

    n >= 3 && twice_area(p1; kernel = kernel) < 0 && reverse!(p1)
    m >= 3 && twice_area(p2; kernel = kernel) < 0 && reverse!(p2)

    scratch_raw = Point{2,R}[]
    scratch_dedup = Point{2,R}[]
    return _minkowski_convex_ccw(p1, p2, kernel, scratch_raw, scratch_dedup)
end

function _minkowski_merge_ccw!(
    raw_out::Vector{Point{2,R}},
    p1::AbstractVector{Point{2,R}},
    p2::AbstractVector{Point{2,R}},
    kernel::AbstractKernel,
) where {R}
    n = length(p1)
    m = length(p2)
    empty!(raw_out)
    sizehint!(raw_out, n + m)

    idx1 = _lowest_left_index(p1)
    idx2 = _lowest_left_index(p2)
    i = 1
    j = 1
    while i <= n || j <= m
        push!(raw_out, p1[idx1] + p2[idx2])
        if i > n
            j += 1
            idx2 = idx2 == m ? 1 : idx2 + 1
        elseif j > m
            i += 1
            idx1 = idx1 == n ? 1 : idx1 + 1
        else
            nxt1 = idx1 == n ? 1 : idx1 + 1
            nxt2 = idx2 == m ? 1 : idx2 + 1
            cmp = compare_angle(kernel, p1[nxt1] - p1[idx1], p2[nxt2] - p2[idx2])
            if cmp <= 0
                i += 1
                idx1 = nxt1
            end
            if cmp >= 0
                j += 1
                idx2 = nxt2
            end
        end
    end
    return raw_out
end

function _minkowski_convex_ccw(
    p1::AbstractVector{Point{2,R}},
    p2::AbstractVector{Point{2,R}},
    kernel::AbstractKernel,
    scratch_raw::Vector{Point{2,R}},
    scratch_dedup::Vector{Point{2,R}},
) where {R}
    n, m = length(p1), length(p2)
    (n == 0 || m == 0) && return Point{2,R}[]
    n == 1 && m == 1 && return Point{2,R}[p1[1]+p2[1]]
    n == 1 && return Point{2,R}[q + p1[1] for q in p2]
    m == 1 && return Point{2,R}[p + p2[1] for p in p1]

    _minkowski_merge_ccw!(scratch_raw, p1, p2, kernel)
    _collect_deduped_points!(scratch_dedup, scratch_raw, kernel, 1e-11)
    length(scratch_dedup) < 3 && return copy(scratch_dedup)

    final_poly = _prune_collinear(scratch_dedup, kernel, 1e-11)
    return isempty(final_poly) ? copy(scratch_dedup) : final_poly
end

# Specialized Point and Segment overloads
minkowski_sum(
    p1::Point{2,T1},
    p2::Point{2,T2};
    kernel::AbstractKernel = FloatingKernel(),
    kwargs...,
) where {T1,T2} = p1 + p2

minkowski_sum(
    p::Point{2,T1},
    s::Segment{2,T2};
    kernel::AbstractKernel = FloatingKernel(),
    kwargs...,
) where {T1,T2} = Segment(p + s.p, p + s.q)

minkowski_sum(
    s::Segment{2,T1},
    p::Point{2,T2};
    kernel::AbstractKernel = FloatingKernel(),
    kwargs...,
) where {T1,T2} = Segment(s.p + p, s.q + p)

minkowski_sum(
    s1::Segment{2},
    s2::Segment{2};
    kernel::AbstractKernel = FloatingKernel(),
    kwargs...,
) = minkowski_sum(_pts_vec(s1), _pts_vec(s2); kernel = kernel, kwargs...)

minkowski_sum(
    p::Point{2},
    poly::AbstractVector{<:Point{2}};
    kernel::AbstractKernel = FloatingKernel(),
    kwargs...,
) = minkowski_sum(_pts_vec(p), poly; kernel = kernel, kwargs...)

minkowski_sum(
    poly::AbstractVector{<:Point{2}},
    p::Point{2};
    kernel::AbstractKernel = FloatingKernel(),
    kwargs...,
) = minkowski_sum(poly, _pts_vec(p); kernel = kernel, kwargs...)

minkowski_sum(
    s::Segment{2},
    poly::AbstractVector{<:Point{2}};
    kernel::AbstractKernel = FloatingKernel(),
    kwargs...,
) = minkowski_sum(_pts_vec(s), poly; kernel = kernel, kwargs...)

minkowski_sum(
    poly::AbstractVector{<:Point{2}},
    s::Segment{2};
    kernel::AbstractKernel = FloatingKernel(),
    kwargs...,
) = minkowski_sum(poly, _pts_vec(s); kernel = kernel, kwargs...)

# PntSeq and generic collections
minkowski_sum(
    p1::AbsPntSeq{2},
    p2::AbsPntSeq{2};
    kernel::AbstractKernel = FloatingKernel(),
    kwargs...,
) = minkowski_sum(Points(p1), Points(p2); kernel = kernel, kwargs...)

minkowski_sum(p1::AbsPntSeq{2}, p2; kernel::AbstractKernel = FloatingKernel(), kwargs...) =
    minkowski_sum(Points(p1), p2; kernel = kernel, kwargs...)

minkowski_sum(p1, p2::AbsPntSeq{2}; kernel::AbstractKernel = FloatingKernel(), kwargs...) =
    minkowski_sum(p1, Points(p2); kernel = kernel, kwargs...)

minkowski_sum(poly1, poly2; kernel::AbstractKernel = FloatingKernel(), kwargs...) =
    minkowski_sum(collect(poly1), collect(poly2); kernel = kernel, kwargs...)

"""
    minkowski_sum_decomposition(poly1, poly2; kernel=FloatingKernel(), algorithm=:sweep)

Compute the Minkowski sum decomposition of two simple polygons into convex pieces.

If either polygon is non-convex, it is triangulated into convex triangles using [`triangulate_polygon`](@ref)
with the selected `algorithm` (`:sweep` by default for deterministic ``O(n \\log n + m \\log m)`` worst-case
triangulation; also supports `:seidel`, `:ear_cut_z_order`, `:ear_cut_kd_tree`, `:ear_cut_u_kd_tree`, and `:ear_cut`).
The Minkowski sum of each pair of convex pieces is computed in linear time.
Returns a `Vector{Vector{Point{2, R}}}` containing all convex pieces whose union is the Minkowski sum
``P \\oplus Q = \\bigcup_{i, j} (P_i \\oplus Q_j)``.

If both inputs are already convex, returns a 1-element vector containing the single convex Minkowski sum polygon.
"""
function minkowski_sum_decomposition(
    poly1::AbstractVector{<:Point{2}},
    poly2::AbstractVector{<:Point{2}};
    kernel::AbstractKernel = FloatingKernel(),
    algorithm::Symbol = :sweep,
)
    R = _minkowski_coord_type(kernel, _scalar_type(poly1), _scalar_type(poly2))
    return _minkowski_sum_decomposition_impl(poly1, poly2, kernel, algorithm, R)
end

function _normalize_convex_copy(p::Vector{Point{2,R}}, kernel::AbstractKernel) where {R}
    out = Point{2,R}[]
    sizehint!(out, length(p))
    _collect_deduped_points!(out, p, kernel, 1e-11)
    length(out) >= 3 && twice_area(out; kernel = kernel) < 0 && reverse!(out)
    return out
end

function _decompose_pieces(
    p::Vector{Point{2,R}},
    kernel::AbstractKernel,
    algorithm::Symbol,
) where {R}
    if length(p) <= 3 || is_convex(p; kernel = kernel)
        return (_normalize_convex_copy(p, kernel), NTuple{3,Int}[])
    end
    t = triangulate_polygon(p; kernel = kernel, algorithm = algorithm)
    # Triangulators may emit clockwise index triples (e.g. :ear_cut), but
    # _minkowski_convex_ccw requires counter-clockwise pieces.
    tris = Vector{NTuple{3,Int}}(undef, length(t.triangles))
    for (k, (i, j, m)) in enumerate(t.triangles)
        tris[k] = orient_sign(kernel, p[i], p[j], p[m]) < 0 ? (i, m, j) : (i, j, m)
    end
    return (p, tris)
end

function _minkowski_sum_decomposition_impl(
    poly1::AbstractVector{<:Point{2}},
    poly2::AbstractVector{<:Point{2}},
    kernel::AbstractKernel,
    algorithm::Symbol,
    ::Type{R},
) where {R}
    p1 = _to_point2_vec(poly1, R)
    p2 = _to_point2_vec(poly2, R)
    pts1, tris1 = _decompose_pieces(p1, kernel, algorithm)
    pts2, tris2 = _decompose_pieces(p2, kernel, algorithm)

    results = Vector{Point{2,R}}[]
    scratch_raw = Point{2,R}[]
    scratch_dedup = Point{2,R}[]
    return _dispatch_p1!(
        results,
        pts1,
        tris1,
        pts2,
        tris2,
        kernel,
        scratch_raw,
        scratch_dedup,
    )
end

function _dispatch_p1!(
    results,
    pts1,
    tris1,
    pts2,
    tris2,
    kernel,
    scratch_raw,
    scratch_dedup,
)
    if isempty(tris1)
        return _dispatch_p2!(
            results,
            (pts1,),
            pts2,
            tris2,
            kernel,
            scratch_raw,
            scratch_dedup,
        )
    end
    return _dispatch_p2!(
        results,
        _TrianglePieces(pts1, tris1),
        pts2,
        tris2,
        kernel,
        scratch_raw,
        scratch_dedup,
    )
end

function _dispatch_p2!(results, pieces1, pts2, tris2, kernel, scratch_raw, scratch_dedup)
    if isempty(tris2)
        return _minkowski_pair_loop!(
            results,
            pieces1,
            (pts2,),
            kernel,
            scratch_raw,
            scratch_dedup,
        )
    end
    return _minkowski_pair_loop!(
        results,
        pieces1,
        _TrianglePieces(pts2, tris2),
        kernel,
        scratch_raw,
        scratch_dedup,
    )
end

function _minkowski_pair_loop!(
    results::Vector{Vector{Point{2,R}}},
    pieces1,
    pieces2,
    kernel::AbstractKernel,
    scratch_raw::Vector{Point{2,R}},
    scratch_dedup::Vector{Point{2,R}},
) where {R}
    sizehint!(results, length(pieces1) * length(pieces2))
    for c1 in pieces1
        for c2 in pieces2
            push!(
                results,
                _minkowski_convex_ccw(c1, c2, kernel, scratch_raw, scratch_dedup),
            )
        end
    end
    return results
end

minkowski_sum_decomposition(
    p1::AbsPntSeq{2},
    p2::AbsPntSeq{2};
    kernel::AbstractKernel = FloatingKernel(),
    kwargs...,
) = minkowski_sum_decomposition(Points(p1), Points(p2); kernel = kernel, kwargs...)

minkowski_sum_decomposition(
    p1::AbsPntSeq{2},
    p2;
    kernel::AbstractKernel = FloatingKernel(),
    kwargs...,
) = minkowski_sum_decomposition(Points(p1), p2; kernel = kernel, kwargs...)

minkowski_sum_decomposition(
    p1,
    p2::AbsPntSeq{2};
    kernel::AbstractKernel = FloatingKernel(),
    kwargs...,
) = minkowski_sum_decomposition(p1, Points(p2); kernel = kernel, kwargs...)

minkowski_sum_decomposition(
    poly1,
    poly2;
    kernel::AbstractKernel = FloatingKernel(),
    kwargs...,
) = minkowski_sum_decomposition(collect(poly1), collect(poly2); kernel = kernel, kwargs...)
