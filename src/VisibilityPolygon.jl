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
- `q` must lie strictly in the interior or on the boundary of `poly`. An `ArgumentError`
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
    area2 > 0 || throw(
        ArgumentError(
            "polygon must be oriented counterclockwise with positive area; got twice_area = $area2",
        ),
    )

    R = _vis_coord_type(kernel, _scalar_type(poly), eltype(q))
    return _visibility_polygon_impl(poly, q, kernel, R)
end

function _visibility_polygon_impl(
    poly::AbstractVector{<:Point{2}},
    q::Point{2},
    kernel::AbstractKernel,
    ::Type{R},
) where {R}
    poly_conv = _to_point2_vec(poly, R)
    q_conv = point(R(q.x), R(q.y))

    loc_type, loc_idx = _point_in_polygon(kernel, q_conv, poly_conv)
    loc_type != :exterior || throw(
        ArgumentError(
            "vantage point q = $q must lie inside or on the boundary of the polygon",
        ),
    )

    unique_rays = _sorted_unique_rays(poly_conv, q_conv, kernel, R)
    m = length(unique_rays)
    start_ray_idx, end_ray_idx =
        _sector_ray_range(kernel, loc_type, loc_idx, poly_conv, q_conv, unique_rays)

    raw_polygon = Point{2,R}[]
    sizehint!(raw_polygon, 2 * m + 1)
    curr_ray_idx = start_ray_idx
    while true
        next_ray_idx = curr_ray_idx == m ? 1 : curr_ray_idx + 1
        hit, best_pt1, best_pt2 = _closest_sector_hit(
            kernel,
            q_conv,
            unique_rays[curr_ray_idx],
            unique_rays[next_ray_idx],
            poly_conv,
            R,
        )
        if hit
            push!(raw_polygon, best_pt1, best_pt2)
        end
        if loc_type != :interior && next_ray_idx == end_ray_idx
            break
        end
        curr_ray_idx = next_ray_idx
        if loc_type == :interior && curr_ray_idx == start_ray_idx
            break
        end
    end

    loc_type != :interior && push!(raw_polygon, q_conv)

    deduped = Point{2,R}[]
    sizehint!(deduped, length(raw_polygon))
    _collect_deduped_points!(deduped, raw_polygon, kernel, 1e-9)
    return _prune_collinear(deduped, kernel, 1e-9)
end

function visibility_polygon(
    poly::AbsPntSeq{2},
    q::Point{2};
    kernel::AbstractKernel = FloatingKernel(),
)
    return visibility_polygon(Points(poly), q; kernel = kernel)
end

function visibility_polygon(poly, q::Point{2}; kernel::AbstractKernel = FloatingKernel())
    return visibility_polygon(collect(poly), q; kernel = kernel)
end

function visibility_polygon(
    poly,
    q::AbstractVector{<:Real};
    kernel::AbstractKernel = FloatingKernel(),
)
    length(q) == 2 || throw(
        ArgumentError("vantage point q must be 2-dimensional; got length = $(length(q))"),
    )
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
    area2 > 0 || throw(
        ArgumentError(
            "polygon must be oriented counterclockwise with positive area; got twice_area = $area2",
        ),
    )

    T = _scalar_type(poly)
    R = _vis_coord_type(kernel, T, T)
    return _polygon_kernel_impl(poly, kernel, R)
end

function _polygon_kernel_impl(
    poly::AbstractVector{<:Point{2}},
    kernel::AbstractKernel,
    ::Type{R},
) where {R}
    poly_conv = _to_point2_vec(poly, R)
    n = length(poly_conv)
    min_x, max_x, min_y, max_y = _bounding_rect(poly_conv)

    curr_poly = Point{2,R}[
        point(min_x, min_y),
        point(max_x, min_y),
        point(max_x, max_y),
        point(min_x, max_y),
    ]
    next_poly = Point{2,R}[]
    deduped = Point{2,R}[]
    sizehint!(curr_poly, 16)
    sizehint!(next_poly, 16)
    sizehint!(deduped, 16)

    tol = 1e-11
    for i = 1:n
        isempty(curr_poly) && break
        a = poly_conv[i]
        b = poly_conv[i == n ? 1 : i + 1]
        _clip_polygon_halfplane!(next_poly, curr_poly, a, b, kernel, R, tol)
        _collect_deduped_points!(deduped, next_poly, kernel, tol)
        curr_poly, deduped = deduped, curr_poly
    end

    length(curr_poly) < 3 && return curr_poly
    return _prune_collinear(curr_poly, kernel, tol)
end

function polygon_kernel(poly::AbsPntSeq{2}; kernel::AbstractKernel = FloatingKernel())
    return polygon_kernel(Points(poly); kernel = kernel)
end

function polygon_kernel(poly; kernel::AbstractKernel = FloatingKernel())
    return polygon_kernel(collect(poly); kernel = kernel)
end

"""
    is_star_shaped(poly; kernel=FloatingKernel())

Return `true` if simple polygon `poly` is star-shaped (i.e. has a non-empty kernel).
"""
function is_star_shaped(poly; kernel::AbstractKernel = FloatingKernel())
    k = polygon_kernel(poly; kernel = kernel)
    return !isempty(k) && (length(k) < 3 || twice_area(k; kernel = kernel) > 0)
end

# Internal helpers

@inline _scalar_type(::Point{2,T}) where {T} = T
@inline _scalar_type(::Segment{2,T}) where {T} = T
@inline _scalar_type(::AbsPntSeq{2,T}) where {T} = T
@inline _scalar_type(::AbstractVector{<:Point{2,T}}) where {T} = T
function _scalar_type(poly)
    isempty(poly) && return Float64
    T = typeof(first(poly).x)
    for p in poly
        T = promote_type(T, typeof(p.x))
    end
    return T
end

@inline _vis_coord_type(::ExactKernel, ::Type{T1}, ::Type{T2}) where {T1,T2} =
    promote_type(T1, T2) <: Rational ? promote_type(T1, T2) : Rational{BigInt}
@inline _vis_coord_type(::AbstractKernel, ::Type{T1}, ::Type{T2}) where {T1,T2} =
    float(promote_type(T1, T2))

@inline _to_point2_vec(poly::Vector{Point{2,R}}, ::Type{R}) where {R} = poly
@inline _to_point2_vec(poly, ::Type{R}) where {R} =
    Point{2,R}[point(R(p.x), R(p.y)) for p in poly]

function _collect_deduped_points!(
    out::Vector{Point{2,R}},
    poly,
    kernel::AbstractKernel,
    tol::Float64,
) where {R}
    empty!(out)
    for raw_p in poly
        p = point(R(raw_p.x), R(raw_p.y))
        if isempty(out)
            push!(out, p)
        else
            is_dup = kernel isa ExactKernel ? (p == last(out)) : (dist(p, last(out)) < tol)
            is_dup || push!(out, p)
        end
    end
    if length(out) > 1
        is_dup =
            kernel isa ExactKernel ? (first(out) == last(out)) :
            (dist(first(out), last(out)) < tol)
        is_dup && pop!(out)
    end
    return out
end

function _prune_collinear(
    pts::Vector{Point{2,R}},
    kernel::AbstractKernel,
    tol::Float64,
) where {R}
    k = length(pts)
    out = Point{2,R}[]
    sizehint!(out, k)
    for i = 1:k
        pr = pts[i == 1 ? k : i - 1]
        cu = pts[i]
        nx = pts[i == k ? 1 : i + 1]
        keep = if kernel isa ExactKernel
            orient_sign(kernel, pr, cu, nx) != 0
        else
            cross = (cu.x - pr.x) * (nx.y - pr.y) - (cu.y - pr.y) * (nx.x - pr.x)
            abs(cross) > tol * max(dist(pr, cu) * dist(cu, nx), 1.0)
        end
        keep && push!(out, cu)
    end
    return out
end

function _sorted_unique_rays(
    poly_conv::Vector{Point{2,R}},
    q_conv::Point{2,R},
    kernel::AbstractKernel,
    ::Type{R},
) where {R}
    directions = Point{2,R}[]
    sizehint!(directions, length(poly_conv))
    for p in poly_conv
        p != q_conv && push!(directions, p - q_conv)
    end
    sort!(directions, lt = (d1, d2) -> angle_less(kernel, d1, d2))

    unique_rays = Point{2,R}[]
    sizehint!(unique_rays, length(directions))
    for d in directions
        if isempty(unique_rays) || compare_angle(kernel, d, last(unique_rays)) != 0
            push!(unique_rays, d)
        end
    end
    length(unique_rays) >= 2 || throw(
        ArgumentError("polygon vertices must span at least 2 distinct directions from q"),
    )
    return unique_rays
end

@inline function _find_ray_index(
    kernel::AbstractKernel,
    unique_rays::Vector{Point{2,R}},
    u::Point{2,R},
) where {R}
    for idx in eachindex(unique_rays)
        compare_angle(kernel, unique_rays[idx], u) == 0 && return idx
    end
    return 0
end

function _sector_ray_range(
    kernel::AbstractKernel,
    loc_type::Symbol,
    loc_idx::Int,
    poly_conv::Vector{Point{2,R}},
    q_conv::Point{2,R},
    unique_rays::Vector{Point{2,R}},
) where {R}
    loc_type == :interior && return (1, 1)
    n = length(poly_conv)
    if loc_type == :vertex
        cu = poly_conv[loc_idx]
        nx = poly_conv[loc_idx == n ? 1 : loc_idx + 1]
        pr = poly_conv[loc_idx == 1 ? n : loc_idx - 1]
        return (
            _find_ray_index(kernel, unique_rays, nx - cu),
            _find_ray_index(kernel, unique_rays, pr - cu),
        )
    end
    a = poly_conv[loc_idx]
    b = poly_conv[loc_idx == n ? 1 : loc_idx + 1]
    return (
        _find_ray_index(kernel, unique_rays, b - q_conv),
        _find_ray_index(kernel, unique_rays, a - q_conv),
    )
end

@inline function _is_closer_sector_hit(
    kernel::AbstractKernel,
    best_edge_idx::Int,
    t1,
    t2,
    best_t1,
    best_t2,
)
    best_edge_idx == 0 && return true
    diff1 = t1 - best_t1
    diff2 = t2 - best_t2
    if kernel isa ExactKernel
        return (diff1 < 0 && diff2 <= 0) || (diff1 <= 0 && diff2 < 0)
    end
    return (diff1 < -1e-9 && diff2 < 1e-9) || (diff2 < -1e-9 && diff1 < 1e-9)
end

function _closest_sector_hit(
    kernel::AbstractKernel,
    q_conv::Point{2,R},
    u1::Point{2,R},
    u2::Point{2,R},
    poly_conv::Vector{Point{2,R}},
    ::Type{R},
) where {R}
    n = length(poly_conv)
    best_edge_idx = 0
    best_t1 = kernel isa ExactKernel ? typemax(R) : R(Inf)
    best_t2 = kernel isa ExactKernel ? typemax(R) : R(Inf)
    best_pt1 = q_conv
    best_pt2 = q_conv

    for e_idx = 1:n
        a = poly_conv[e_idx]
        b = poly_conv[e_idx == n ? 1 : e_idx + 1]
        hit1, t1, pt1 = _ray_edge_intersect(kernel, q_conv, u1, a, b, R)
        hit1 || continue
        hit2, t2, pt2 = _ray_edge_intersect(kernel, q_conv, u2, a, b, R)
        hit2 || continue

        if _is_closer_sector_hit(kernel, best_edge_idx, t1, t2, best_t1, best_t2)
            best_edge_idx = e_idx
            best_t1 = t1
            best_t2 = t2
            best_pt1 = pt1
            best_pt2 = pt2
        end
    end
    return (best_edge_idx != 0, best_pt1, best_pt2)
end

function _bounding_rect(pts::Vector{Point{2,R}}) where {R}
    p1 = pts[1]
    min_x = max_x = p1.x
    min_y = max_y = p1.y
    for idx = 2:length(pts)
        p = pts[idx]
        min_x = min(min_x, p.x)
        max_x = max(max_x, p.x)
        min_y = min(min_y, p.y)
        max_y = max(max_y, p.y)
    end
    return (min_x, max_x, min_y, max_y)
end

@inline _halfplane_thresh(::ExactKernel, ab, ::Float64) = 0.0
@inline _halfplane_thresh(::AbstractKernel, ab, tol::Float64) =
    -tol * max(Float64(hypot(ab.x, ab.y)), 1.0)

@inline _halfplane_inside(kernel::ExactKernel, a, b, p, ::Float64) =
    turn_sign(kernel, a, b, p) >= 0
@inline _halfplane_inside(kernel::AbstractKernel, a, b, p, thresh::Float64) =
    turn_sign(kernel, a, b, p) >= thresh

@inline function _halfplane_edge_intersect(
    a::Point{2,R},
    ab::Point{2,R},
    p_prev::Point{2,R},
    p_curr::Point{2,R},
    kernel::AbstractKernel,
    ::Type{R},
) where {R}
    d = p_curr - p_prev
    denom = d.x * ab.y - d.y * ab.x
    iszero(denom) && return (false, p_prev)
    ap = a - p_prev
    num = ap.x * ab.y - ap.y * ab.x
    t = kernel isa ExactKernel ? (num // denom) : (num / denom)
    t_clamped = clamp(t, zero(R), one(R))
    return (true, point(R(p_prev.x + t_clamped * d.x), R(p_prev.y + t_clamped * d.y)))
end

function _clip_polygon_halfplane!(
    next_poly::Vector{Point{2,R}},
    curr_poly::Vector{Point{2,R}},
    a::Point{2,R},
    b::Point{2,R},
    kernel::AbstractKernel,
    ::Type{R},
    tol::Float64,
) where {R}
    empty!(next_poly)
    m = length(curr_poly)
    ab = b - a
    thresh = _halfplane_thresh(kernel, ab, tol)

    p_prev = curr_poly[m]
    prev_in = _halfplane_inside(kernel, a, b, p_prev, thresh)
    for j = 1:m
        p_curr = curr_poly[j]
        curr_in = _halfplane_inside(kernel, a, b, p_curr, thresh)
        if curr_in != prev_in
            hit, pt = _halfplane_edge_intersect(a, ab, p_prev, p_curr, kernel, R)
            hit && push!(next_poly, pt)
        end
        curr_in && push!(next_poly, p_curr)
        p_prev = p_curr
        prev_in = curr_in
    end
    return next_poly
end

@inline function _ray_edge_intersect(
    kernel::AbstractKernel,
    q::Point{2,R},
    u::Point{2,R},
    a::Point{2,R},
    b::Point{2,R},
    ::Type{R},
) where {R}
    s_a = orient_sign(kernel, q, q + u, a)
    s_b = orient_sign(kernel, q, q + u, b)

    if (s_a > 0 && s_b > 0) || (s_a < 0 && s_b < 0)
        return (false, zero(R), q)
    end

    ba_x = b.x - a.x
    ba_y = b.y - a.y
    denom = u.x * ba_y - u.y * ba_x
    iszero(denom) && return (false, zero(R), q)

    aq_x = a.x - q.x
    aq_y = a.y - q.y
    num_t = aq_x * ba_y - aq_y * ba_x
    num_s = u.x * aq_y - u.y * aq_x

    if kernel isa ExactKernel
        t = num_t // denom
        s = -num_s // denom
        if t > 0 && 0 <= s <= 1
            pt = point(R(q.x + t * u.x), R(q.y + t * u.y))
            return (true, R(t), pt)
        end
    else
        t = num_t / denom
        s = -num_s / denom
        if t > 1e-11 && -1e-11 <= s <= 1.0 + 1e-11
            s_clamped = clamp(s, zero(t), one(t))
            pt = point(R(a.x + s_clamped * ba_x), R(a.y + s_clamped * ba_y))
            return (true, R(t), pt)
        end
    end
    return (false, zero(R), q)
end

@inline function _point_in_polygon(kernel::AbstractKernel, q::Point{2}, poly)
    n = length(poly)
    for i = 1:n
        q == poly[i] && return (:vertex, i)
    end

    wn = 0
    for i = 1:n
        a = poly[i]
        b = poly[i == n ? 1 : i + 1]

        s = orient_sign(kernel, a, b, q)
        if s == 0
            min_x, max_x = min(a.x, b.x), max(a.x, b.x)
            min_y, max_y = min(a.y, b.y), max(a.y, b.y)
            if min_x <= q.x <= max_x && min_y <= q.y <= max_y
                return (:edge, i)
            end
        end

        if a.y <= q.y
            b.y > q.y && s > 0 && (wn += 1)
        else
            b.y <= q.y && s < 0 && (wn -= 1)
        end
    end

    return wn != 0 ? (:interior, 0) : (:exterior, 0)
end
