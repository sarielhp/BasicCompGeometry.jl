# Minkowski Sum of 2D Polygons and Shapes

@inline _scalar_type(p::Point) = typeof(p.x)
@inline _scalar_type(s::Segment) = typeof(s.p.x)
function _scalar_type(poly)
    for p in poly
        return typeof(p.x)
    end
    return Float64
end

function _dedup_consecutive(pts::Vector{Point{2, R}}, kernel::AbstractKernel) where {R}
    n = length(pts)
    n <= 1 && return pts
    tol = 1e-11
    res = Point{2, R}[]
    for p in pts
        if isempty(res)
            push!(res, p)
        else
            is_dup = kernel isa ExactKernel ? (p == last(res)) : (dist(p, last(res)) < tol)
            if !is_dup
                push!(res, p)
            end
        end
    end
    if length(res) > 1
        is_dup = kernel isa ExactKernel ? (first(res) == last(res)) : (dist(first(res), last(res)) < tol)
        if is_dup
            pop!(res)
        end
    end
    return res
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
    n = length(poly)
    n <= 2 && return true

    R = _scalar_type(poly)
    if kernel isa ExactKernel && !(R <: Rational)
        R = Rational{BigInt}
    end
    p_conv = [point(R(p.x), R(p.y)) for p in poly]
    deduped = _dedup_consecutive(p_conv, kernel)

    k = length(deduped)
    k <= 2 && return true

    area2 = twice_area(deduped; kernel = kernel)
    if iszero(area2)
        return false
    end

    if area2 < 0
        deduped = reverse(deduped)
    end

    tol = 1e-11
    for i in 1:k
        p = deduped[mod1(i - 1, k)]
        q = deduped[i]
        r = deduped[mod1(i + 1, k)]

        ts = orient_sign(kernel, p, q, r)
        if kernel isa ExactKernel
            ts < 0 && return false
        else
            if ts < 0
                d1 = dist(p, q)
                d2 = dist(q, r)
                cross = (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x)
                if cross < -tol * max(d1 * d2, 1.0)
                    return false
                end
            end
        end
    end

    idx_min = 1
    for i in 2:k
        if deduped[i].y < deduped[idx_min].y || (deduped[i].y == deduped[idx_min].y && deduped[i].x < deduped[idx_min].x)
            idx_min = i
        end
    end

    edges = [deduped[mod1(idx_min + idx, k)] - deduped[mod1(idx_min + idx - 1, k)] for idx in 1:k]
    nz_edges = filter(e -> !(iszero(e.x) && iszero(e.y)), edges)
    length(nz_edges) < 3 && return false

    for idx in 1:(length(nz_edges) - 1)
        cmp = compare_angle(kernel, nz_edges[idx], nz_edges[idx + 1])
        if cmp > 0
            return false
        end
    end

    return true
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
        return _minkowski_convex(h1, h2, kernel)
    end

    if !is_convex(poly1; kernel = kernel)
        throw(ArgumentError("minkowski_sum: poly1 is not convex. Use minkowski_sum_decomposition for non-convex polygons, or set hull=true to compute the convex hull of the sum."))
    end
    if !is_convex(poly2; kernel = kernel)
        throw(ArgumentError("minkowski_sum: poly2 is not convex. Use minkowski_sum_decomposition for non-convex polygons, or set hull=true to compute the convex hull of the sum."))
    end

    return _minkowski_convex(poly1, poly2, kernel)
end

function _minkowski_convex(
    poly1,
    poly2,
    kernel::AbstractKernel,
)
    R = promote_type(_scalar_type(poly1), _scalar_type(poly2))
    if kernel isa ExactKernel && !(R <: Rational)
        R = Rational{BigInt}
    end

    p1 = Point{2, R}[point(R(p.x), R(p.y)) for p in poly1]
    p2 = Point{2, R}[point(R(p.x), R(p.y)) for p in poly2]

    p1 = _dedup_consecutive(p1, kernel)
    p2 = _dedup_consecutive(p2, kernel)

    n = length(p1)
    m = length(p2)

    if n == 0 || m == 0
        return Point{2, R}[]
    end

    if n == 1 && m == 1
        return Point{2, R}[p1[1] + p2[1]]
    elseif n == 1
        return Point{2, R}[q + p1[1] for q in p2]
    elseif m == 1
        return Point{2, R}[p + p2[1] for p in p1]
    end

    if n >= 3 && twice_area(p1; kernel = kernel) < 0
        p1 = reverse(p1)
    end
    if m >= 3 && twice_area(p2; kernel = kernel) < 0
        p2 = reverse(p2)
    end

    function find_start(pts)
        idx = 1
        for i in 2:length(pts)
            if pts[i].y < pts[idx].y || (pts[i].y == pts[idx].y && pts[i].x < pts[idx].x)
                idx = i
            end
        end
        return idx
    end

    i0 = find_start(p1)
    j0 = find_start(p2)

    v1 = [p1[mod1(i0 + k - 1, n)] for k in 1:n]
    v2 = [p2[mod1(j0 + k - 1, m)] for k in 1:m]

    edges1 = [v1[mod1(k + 1, n)] - v1[k] for k in 1:n]
    edges2 = [v2[mod1(k + 1, m)] - v2[k] for k in 1:m]

    result = Point{2, R}[]
    i = 1
    j = 1

    while i <= n || j <= m
        push!(result, v1[mod1(i, n)] + v2[mod1(j, m)])

        if i > n
            j += 1
        elseif j > m
            i += 1
        else
            e1 = edges1[i]
            e2 = edges2[j]

            cmp = compare_angle(kernel, e1, e2)
            if cmp < 0
                i += 1
            elseif cmp > 0
                j += 1
            else
                i += 1
                j += 1
            end
        end
    end

    deduped = _dedup_consecutive(result, kernel)
    k = length(deduped)
    if k < 3
        return deduped
    end

    tol = 1e-11
    final_poly = Point{2, R}[]
    for idx in 1:k
        pr = deduped[mod1(idx - 1, k)]
        cu = deduped[idx]
        nx = deduped[mod1(idx + 1, k)]
        if kernel isa ExactKernel
            if orient_sign(kernel, pr, cu, nx) != 0
                push!(final_poly, cu)
            end
        else
            cross = (cu.x - pr.x) * (nx.y - pr.y) - (cu.y - pr.y) * (nx.x - pr.x)
            d1 = dist(pr, cu)
            d2 = dist(cu, nx)
            if abs(cross) > tol * max(d1 * d2, 1.0)
                push!(final_poly, cu)
            end
        end
    end

    return isempty(final_poly) ? deduped : final_poly
end

# Specialized Point and Segment overloads
minkowski_sum(p1::Point{2, T1}, p2::Point{2, T2}; kernel::AbstractKernel = FloatingKernel(), kwargs...) where {T1, T2} =
    p1 + p2

minkowski_sum(p::Point{2, T1}, s::Segment{2, T2}; kernel::AbstractKernel = FloatingKernel(), kwargs...) where {T1, T2} =
    Segment(p + s.p, p + s.q)

minkowski_sum(s::Segment{2, T1}, p::Point{2, T2}; kernel::AbstractKernel = FloatingKernel(), kwargs...) where {T1, T2} =
    Segment(s.p + p, s.q + p)

minkowski_sum(s1::Segment{2}, s2::Segment{2}; kernel::AbstractKernel = FloatingKernel(), kwargs...) =
    minkowski_sum(Point{2}[s1.p, s1.q], Point{2}[s2.p, s2.q]; kernel = kernel, kwargs...)

minkowski_sum(p::Point{2}, poly::AbstractVector{<:Point{2}}; kernel::AbstractKernel = FloatingKernel(), kwargs...) =
    minkowski_sum(Point{2}[p], poly; kernel = kernel, kwargs...)

minkowski_sum(poly::AbstractVector{<:Point{2}}, p::Point{2}; kernel::AbstractKernel = FloatingKernel(), kwargs...) =
    minkowski_sum(poly, Point{2}[p]; kernel = kernel, kwargs...)

minkowski_sum(s::Segment{2}, poly::AbstractVector{<:Point{2}}; kernel::AbstractKernel = FloatingKernel(), kwargs...) =
    minkowski_sum(Point{2}[s.p, s.q], poly; kernel = kernel, kwargs...)

minkowski_sum(poly::AbstractVector{<:Point{2}}, s::Segment{2}; kernel::AbstractKernel = FloatingKernel(), kwargs...) =
    minkowski_sum(poly, Point{2}[s.p, s.q]; kernel = kernel, kwargs...)

# PntSeq and generic collections
minkowski_sum(p1::AbsPntSeq{2}, p2::AbsPntSeq{2}; kernel::AbstractKernel = FloatingKernel(), kwargs...) =
    minkowski_sum(Points(p1), Points(p2); kernel = kernel, kwargs...)

minkowski_sum(p1::AbsPntSeq{2}, p2; kernel::AbstractKernel = FloatingKernel(), kwargs...) =
    minkowski_sum(Points(p1), p2; kernel = kernel, kwargs...)

minkowski_sum(p1, p2::AbsPntSeq{2}; kernel::AbstractKernel = FloatingKernel(), kwargs...) =
    minkowski_sum(p1, Points(p2); kernel = kernel, kwargs...)

minkowski_sum(poly1, poly2; kernel::AbstractKernel = FloatingKernel(), kwargs...) =
    minkowski_sum(collect(poly1), collect(poly2); kernel = kernel, kwargs...)

"""
    minkowski_sum_decomposition(poly1, poly2; kernel=FloatingKernel())

Compute the Minkowski sum decomposition of two simple polygons into convex pieces.

If either polygon is non-convex, it is triangulated into convex triangles using [`triangulate_polygon`](@ref).
The Minkowski sum of each pair of convex pieces is computed in linear time.
Returns a `Vector{Vector{Point{2, R}}}` containing all convex pieces whose union is the Minkowski sum
``P \\oplus Q = \\bigcup_{i, j} (P_i \\oplus Q_j)``.

If both inputs are already convex, returns a 1-element vector containing the single convex Minkowski sum polygon.
"""
function minkowski_sum_decomposition(
    poly1::AbstractVector{<:Point{2}},
    poly2::AbstractVector{<:Point{2}};
    kernel::AbstractKernel = FloatingKernel(),
)
    R = promote_type(_scalar_type(poly1), _scalar_type(poly2))
    if kernel isa ExactKernel && !(R <: Rational)
        R = Rational{BigInt}
    end

    p1 = [point(R(p.x), R(p.y)) for p in poly1]
    p2 = [point(R(p.x), R(p.y)) for p in poly2]

    # Decompose p1
    pieces1 = if length(p1) <= 3 || is_convex(p1; kernel = kernel)
        Vector{Point{2, R}}[p1]
    else
        t1 = triangulate_polygon(p1; kernel = kernel)
        Vector{Point{2, R}}[[p1[i], p1[j], p1[m]] for (i, j, m) in t1.triangles]
    end

    # Decompose p2
    pieces2 = if length(p2) <= 3 || is_convex(p2; kernel = kernel)
        Vector{Point{2, R}}[p2]
    else
        t2 = triangulate_polygon(p2; kernel = kernel)
        Vector{Point{2, R}}[[p2[i], p2[j], p2[m]] for (i, j, m) in t2.triangles]
    end

    results = Vector{Point{2, R}}[]
    for c1 in pieces1
        for c2 in pieces2
            push!(results, _minkowski_convex(c1, c2, kernel))
        end
    end
    return results
end

minkowski_sum_decomposition(p1::AbsPntSeq{2}, p2::AbsPntSeq{2}; kernel::AbstractKernel = FloatingKernel()) =
    minkowski_sum_decomposition(Points(p1), Points(p2); kernel = kernel)

minkowski_sum_decomposition(p1::AbsPntSeq{2}, p2; kernel::AbstractKernel = FloatingKernel()) =
    minkowski_sum_decomposition(Points(p1), p2; kernel = kernel)

minkowski_sum_decomposition(p1, p2::AbsPntSeq{2}; kernel::AbstractKernel = FloatingKernel()) =
    minkowski_sum_decomposition(p1, Points(p2); kernel = kernel)

minkowski_sum_decomposition(poly1, poly2; kernel::AbstractKernel = FloatingKernel()) =
    minkowski_sum_decomposition(collect(poly1), collect(poly2); kernel = kernel)
