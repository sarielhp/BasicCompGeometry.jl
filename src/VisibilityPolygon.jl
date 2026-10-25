# Visibility Polygon from a Vantage Point inside a Simple Polygon

"""
    visibility_polygon(poly, q; kernel=FloatingKernel())

Compute the star-shaped visibility polygon ``\\mathrm{Vis}(q, P)`` consisting of all points inside
a counterclockwise simple polygon ``P`` that are visible from vantage point ``q \\in P``.

A point ``p \\in P`` is visible from ``q`` if the open line segment between ``q`` and ``p`` does not
intersect the exterior of ``P`` or any boundary edge of ``P``.

# Arguments
- `poly`: vertices of a counterclockwise simple polygon, supplied as an `AbstractVector{<:Point{2}}`
  or an [`AbsPntSeq{2}`](@ref).
- `q`: vantage point ([`Point{2}`](@ref) or 2-element coordinate vector) located in the interior or
  on the boundary of `poly`.
- `kernel`: computational geometry kernel ([`FloatingKernel()`](@ref) or [`ExactKernel()`](@ref)).
  Default is `FloatingKernel()`. Under `ExactKernel()`, rational coordinates and exact arithmetic
  are used throughout.

# Returns
- `Vector{Point{2, R}}`: vertices of the star-shaped visibility polygon in counterclockwise order,
  with collinear intermediate vertices pruned.

# Preconditions
- `poly` must have at least 3 vertices and be oriented counterclockwise with positive area
  (as verified by [`twice_area`](@ref)).
- `q` must lie strictly in the interior or on the boundary of `poly`. An [`ArgumentError`](@ref)
  is thrown if `q` is outside `poly`.

# Examples
```julia
using BasicCompGeometry

# U-shaped polygon
p_u = [
    point(0.0, 0.0), point(4.0, 0.0), point(4.0, 4.0), point(3.0, 4.0),
    point(3.0, 1.0), point(1.0, 1.0), point(1.0, 4.0), point(0.0, 4.0)
]

# Query from bottom corridor
vis = visibility_polygon(p_u, point(2.0, 0.5))
length(vis) # 6 vertices
twice_area(vis) # 9.0
```
"""
function visibility_polygon(
    poly::AbstractVector{<:Point{2}},
    q::Point{2};
    kernel::AbstractKernel = FloatingKernel(),
)
    n = length(poly)
    n >= 3 || throw(ArgumentError("polygon must have at least 3 vertices; got n = $n"))
    area2 = twice_area(poly; kernel = kernel)
    area2 > 0 || throw(ArgumentError("polygon must be oriented counterclockwise with positive area; got twice_area = $area2"))

    T_poly = eltype(eltype(poly))
    T_q = eltype(q)
    R = promote_type(T_poly, T_q)
    if kernel isa ExactKernel && !(R <: Rational)
        R = Rational{BigInt}
    end

    poly_conv = [point(R(p.x), R(p.y)) for p in poly]
    q_conv = point(R(q.x), R(q.y))

    loc_type, loc_idx = _point_in_polygon(kernel, q_conv, poly_conv)
    loc_type != :exterior || throw(ArgumentError("vantage point q = $q must lie inside or on the boundary of the polygon"))

    # Critical rays from q to all other polygon vertices
    directions = Point{2, R}[]
    for p in poly_conv
        if p != q_conv
            push!(directions, p - q_conv)
        end
    end

    sort!(directions, lt = (d1, d2) -> angle_less(kernel, d1, d2))

    unique_rays = Point{2, R}[]
    for d in directions
        if isempty(unique_rays) || compare_angle(kernel, d, last(unique_rays)) != 0
            push!(unique_rays, d)
        end
    end
    m = length(unique_rays)
    m >= 2 || throw(ArgumentError("polygon vertices must span at least 2 distinct directions from q"))

    edges = [(poly_conv[j], poly_conv[mod1(j + 1, n)]) for j in 1:n]

    # Determine sector range in CCW angular order
    local start_ray_idx::Int
    local end_ray_idx::Int

    if loc_type == :interior
        start_ray_idx = 1
        end_ray_idx = 1
    elseif loc_type == :vertex
        cu = poly_conv[loc_idx]
        nx = poly_conv[mod1(loc_idx + 1, n)]
        pr = poly_conv[mod1(loc_idx - 1, n)]
        u_start = nx - cu
        u_end = pr - cu
        start_ray_idx = findfirst(r -> compare_angle(kernel, r, u_start) == 0, unique_rays)
        end_ray_idx = findfirst(r -> compare_angle(kernel, r, u_end) == 0, unique_rays)
    elseif loc_type == :edge
        a = poly_conv[loc_idx]
        b = poly_conv[mod1(loc_idx + 1, n)]
        u_start = b - q_conv
        u_end = a - q_conv
        start_ray_idx = findfirst(r -> compare_angle(kernel, r, u_start) == 0, unique_rays)
        end_ray_idx = findfirst(r -> compare_angle(kernel, r, u_end) == 0, unique_rays)
    end

    raw_polygon = Point{2, R}[]
    curr_ray_idx = start_ray_idx

    while true
        next_ray_idx = mod1(curr_ray_idx + 1, m)
        u1 = unique_rays[curr_ray_idx]
        u2 = unique_rays[next_ray_idx]

        best_edge_idx = 0
        best_t1 = kernel isa ExactKernel ? typemax(Rational{BigInt}) : Inf
        best_t2 = kernel isa ExactKernel ? typemax(Rational{BigInt}) : Inf
        best_pt1 = q_conv
        best_pt2 = q_conv

        for (e_idx, (a, b)) in enumerate(edges)
            hit1, t1, pt1 = _ray_edge_intersect(kernel, q_conv, u1, a, b, R)
            hit2, t2, pt2 = _ray_edge_intersect(kernel, q_conv, u2, a, b, R)

            if hit1 && hit2
                is_closer = false
                if best_edge_idx == 0
                    is_closer = true
                else
                    diff1 = t1 - best_t1
                    diff2 = t2 - best_t2
                    if kernel isa ExactKernel
                        if (diff1 < 0 && diff2 <= 0) || (diff1 <= 0 && diff2 < 0)
                            is_closer = true
                        end
                    else
                        if (diff1 < -1e-9 && diff2 < 1e-9) || (diff2 < -1e-9 && diff1 < 1e-9)
                            is_closer = true
                        end
                    end
                end

                if is_closer
                    best_edge_idx = e_idx
                    best_t1 = t1
                    best_t2 = t2
                    best_pt1 = pt1
                    best_pt2 = pt2
                end
            end
        end

        if best_edge_idx != 0
            push!(raw_polygon, best_pt1)
            push!(raw_polygon, best_pt2)
        end

        if loc_type != :interior && next_ray_idx == end_ray_idx
            break
        end
        curr_ray_idx = next_ray_idx
        if loc_type == :interior && curr_ray_idx == start_ray_idx
            break
        end
    end

    if loc_type != :interior
        push!(raw_polygon, q_conv)
    end

    # Deduplicate consecutive identical points
    deduped = Point{2, R}[]
    for p in raw_polygon
        if isempty(deduped)
            push!(deduped, p)
        else
            is_dup = kernel isa ExactKernel ? (p == last(deduped)) : (dist(p, last(deduped)) < 1e-9)
            if !is_dup
                push!(deduped, p)
            end
        end
    end
    if length(deduped) > 1
        is_dup = kernel isa ExactKernel ? (first(deduped) == last(deduped)) : (dist(first(deduped), last(deduped)) < 1e-9)
        if is_dup
            pop!(deduped)
        end
    end

    # Remove collinear intermediate vertices along edges
    k = length(deduped)
    final_vertices = Point{2, R}[]
    for i in 1:k
        pr = deduped[mod1(i - 1, k)]
        cu = deduped[i]
        nx = deduped[mod1(i + 1, k)]
        if kernel isa ExactKernel
            s = orient_sign(kernel, pr, cu, nx)
            if s != 0
                push!(final_vertices, cu)
            end
        else
            cross = (cu.x - pr.x) * (nx.y - pr.y) - (cu.y - pr.y) * (nx.x - pr.x)
            d_prev = dist(pr, cu)
            d_next = dist(cu, nx)
            if abs(cross) > 1e-9 * max(d_prev * d_next, 1.0)
                push!(final_vertices, cu)
            end
        end
    end

    return final_vertices
end

function visibility_polygon(
    poly::AbsPntSeq{2},
    q::Point{2};
    kernel::AbstractKernel = FloatingKernel(),
)
    return visibility_polygon(Points(poly), q; kernel = kernel)
end

function visibility_polygon(
    poly,
    q::Point{2};
    kernel::AbstractKernel = FloatingKernel(),
)
    return visibility_polygon(collect(poly), q; kernel = kernel)
end

function visibility_polygon(
    poly,
    q::AbstractVector{<:Real};
    kernel::AbstractKernel = FloatingKernel(),
)
    length(q) == 2 || throw(ArgumentError("vantage point q must be 2-dimensional; got length = $(length(q))"))
    return visibility_polygon(poly, point(q[1], q[2]); kernel = kernel)
end

"""
    is_inside(q::Point{2}, poly; kernel=FloatingKernel())

Return `true` if query point `q` lies strictly inside or on the boundary of the simple polygon `poly`.
"""
function is_inside(
    q::Point{2},
    poly::AbstractVector{<:Point{2}};
    kernel::AbstractKernel = FloatingKernel(),
)
    loc_type, _ = _point_in_polygon(kernel, q, poly)
    return loc_type != :exterior
end

function is_inside(
    q::Point{2},
    poly::AbsPntSeq{2};
    kernel::AbstractKernel = FloatingKernel(),
)
    return is_inside(q, Points(poly); kernel = kernel)
end

"""
    polygon_kernel(poly; kernel=FloatingKernel())

Compute the kernel of a simple counterclockwise polygon ``P``, denoted ``\\mathrm{Ker}(P)``.

The kernel is the set of all points in ``P`` from which the entire polygon is visible:
``\\mathrm{Ker}(P) = \\{ q \\in P \\mid \\forall p \\in P, \\overline{qp} \\subset P \\}``.
It is the intersection of all closed inward halfplanes supported by the directed edges of ``P``.

# Arguments
- `poly`: vertices of a counterclockwise simple polygon, supplied as an `AbstractVector{<:Point{2}}`
  or an [`AbsPntSeq{2}`](@ref).
- `kernel`: computational geometry kernel ([`FloatingKernel()`](@ref) or [`ExactKernel()`](@ref)).
  Default is `FloatingKernel()`.

# Returns
- `Vector{Point{2, R}}`: vertices of the convex polygon forming the kernel in counterclockwise order,
  or an empty vector if the polygon has an empty kernel (i.e. is not star-shaped).

# Preconditions
- `poly` must have at least 3 vertices and be oriented counterclockwise with positive area.
"""
function polygon_kernel(
    poly::AbstractVector{<:Point{2}};
    kernel::AbstractKernel = FloatingKernel(),
)
    n = length(poly)
    n >= 3 || throw(ArgumentError("polygon must have at least 3 vertices; got n = $n"))
    area2 = twice_area(poly; kernel = kernel)
    area2 > 0 || throw(ArgumentError("polygon must be oriented counterclockwise with positive area; got twice_area = $area2"))

    T = eltype(eltype(poly))
    R = T
    if kernel isa ExactKernel && !(R <: Rational)
        R = Rational{BigInt}
    end
    poly_conv = [point(R(p.x), R(p.y)) for p in poly]

    min_x = minimum(p.x for p in poly_conv)
    max_x = maximum(p.x for p in poly_conv)
    min_y = minimum(p.y for p in poly_conv)
    max_y = maximum(p.y for p in poly_conv)

    curr_poly = Point{2, R}[
        point(min_x, min_y),
        point(max_x, min_y),
        point(max_x, max_y),
        point(min_x, max_y),
    ]

    tol = 1e-11

    for i in 1:n
        a = poly_conv[i]
        b = poly_conv[mod1(i + 1, n)]

        isempty(curr_poly) && break
        next_poly = Point{2, R}[]
        m = length(curr_poly)

        ab = b - a
        ab_len = kernel isa ExactKernel ? one(R) : max(hypot(ab.x, ab.y), 1.0)

        for j in 1:m
            p_prev = curr_poly[mod1(j - 1, m)]
            p_curr = curr_poly[j]

            ts_prev = turn_sign(kernel, a, b, p_prev)
            ts_curr = turn_sign(kernel, a, b, p_curr)

            if kernel isa ExactKernel
                prev_in = ts_prev >= 0
                curr_in = ts_curr >= 0
            else
                prev_in = ts_prev >= -tol * ab_len
                curr_in = ts_curr >= -tol * ab_len
            end

            if curr_in
                if !prev_in
                    d = p_curr - p_prev
                    denom = d.x * ab.y - d.y * ab.x
                    if !iszero(denom)
                        ap = a - p_prev
                        num = ap.x * ab.y - ap.y * ab.x
                        t = kernel isa ExactKernel ? (num // denom) : (num / denom)
                        t_clamped = clamp(t, zero(R), one(R))
                        pt = point(p_prev.x + t_clamped * d.x, p_prev.y + t_clamped * d.y)
                        push!(next_poly, pt)
                    end
                end
                push!(next_poly, p_curr)
            elseif prev_in
                d = p_curr - p_prev
                denom = d.x * ab.y - d.y * ab.x
                if !iszero(denom)
                    ap = a - p_prev
                    num = ap.x * ab.y - ap.y * ab.x
                    t = kernel isa ExactKernel ? (num // denom) : (num / denom)
                    t_clamped = clamp(t, zero(R), one(R))
                    pt = point(p_prev.x + t_clamped * d.x, p_prev.y + t_clamped * d.y)
                    push!(next_poly, pt)
                end
            end
        end

        # Deduplicate consecutive points
        deduped = Point{2, R}[]
        for p in next_poly
            if isempty(deduped)
                push!(deduped, p)
            else
                is_dup = kernel isa ExactKernel ? (p == last(deduped)) : (dist(p, last(deduped)) < tol)
                if !is_dup
                    push!(deduped, p)
                end
            end
        end
        if length(deduped) > 1
            is_dup = kernel isa ExactKernel ? (first(deduped) == last(deduped)) : (dist(first(deduped), last(deduped)) < tol)
            if is_dup
                pop!(deduped)
            end
        end
        curr_poly = deduped
    end

    k = length(curr_poly)
    if k < 3
        return curr_poly
    end

    # Remove collinear intermediate vertices
    final_poly = Point{2, R}[]
    for i in 1:k
        pr = curr_poly[mod1(i - 1, k)]
        cu = curr_poly[i]
        nx = curr_poly[mod1(i + 1, k)]
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

    return final_poly
end

function polygon_kernel(
    poly::AbsPntSeq{2};
    kernel::AbstractKernel = FloatingKernel(),
)
    return polygon_kernel(Points(poly); kernel = kernel)
end

function polygon_kernel(
    poly;
    kernel::AbstractKernel = FloatingKernel(),
)
    return polygon_kernel(collect(poly); kernel = kernel)
end

"""
    is_star_shaped(poly; kernel=FloatingKernel())

Return `true` if simple polygon `poly` is star-shaped (i.e. has a non-empty kernel).
"""
function is_star_shaped(
    poly;
    kernel::AbstractKernel = FloatingKernel(),
)
    k = polygon_kernel(poly; kernel = kernel)
    return !isempty(k) && (length(k) < 3 || twice_area(k; kernel = kernel) > 0)
end

# Internal helpers

@inline function _ray_edge_intersect(kernel::AbstractKernel, q::Point{2, R}, u::Point{2, R}, a::Point{2, R}, b::Point{2, R}, ::Type{R}) where {R}
    s_a = orient_sign(kernel, q, q + u, a)
    s_b = orient_sign(kernel, q, q + u, b)

    if (s_a > 0 && s_b > 0) || (s_a < 0 && s_b < 0)
        return (false, zero(R), q)
    end

    ba_x = b.x - a.x
    ba_y = b.y - a.y
    denom = u.x * ba_y - u.y * ba_x

    if iszero(denom)
        return (false, zero(R), q)
    end

    aq_x = a.x - q.x
    aq_y = a.y - q.y
    num_t = aq_x * ba_y - aq_y * ba_x
    num_s = u.x * aq_y - u.y * aq_x

    if kernel isa ExactKernel
        t = num_t // denom
        s = -num_s // denom
        if t > 0 && 0 <= s <= 1
            pt = point(q.x + t * u.x, q.y + t * u.y)
            return (true, t, pt)
        end
    else
        t = num_t / denom
        s = -num_s / denom
        if t > 1e-11 && -1e-11 <= s <= 1.0 + 1e-11
            s_clamped = clamp(s, 0.0, 1.0)
            pt = point(a.x + s_clamped * ba_x, a.y + s_clamped * ba_y)
            return (true, t, pt)
        end
    end
    return (false, zero(R), q)
end

@inline function _point_in_polygon(kernel::AbstractKernel, q::Point{2}, poly)
    n = length(poly)
    for i in 1:n
        if q == poly[i]
            return (:vertex, i)
        end
    end

    wn = 0
    for i in 1:n
        a = poly[i]
        b = poly[mod1(i + 1, n)]

        s = orient_sign(kernel, a, b, q)
        if s == 0
            min_x, max_x = min(a.x, b.x), max(a.x, b.x)
            min_y, max_y = min(a.y, b.y), max(a.y, b.y)
            if min_x <= q.x <= max_x && min_y <= q.y <= max_y
                return (:edge, i)
            end
        end

        if a.y <= q.y
            if b.y > q.y && orient_sign(kernel, a, b, q) > 0
                wn += 1
            end
        else
            if b.y <= q.y && orient_sign(kernel, a, b, q) < 0
                wn -= 1
            end
        end
    end

    return wn != 0 ? (:interior, 0) : (:exterior, 0)
end
