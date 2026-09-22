###############################################
### Line and Halfplane Bundles, Vertices, and Default View

"""
    intersect(l1::Line{2}, l2::Line{2}; tol::Real = 1e-12)

Compute the unique intersection point of two 2D infinite lines `l1` and `l2`.
Returns `Point{2, Float64}` if lines intersect at a single point, or `nothing` if they are
parallel or collinear within relative tolerance `tol`.
"""
function Base.intersect(l1::Line{2}, l2::Line{2}; tol::Real = 1e-12)
    # Line 1: l1.p + s * l1.u
    # Line 2: l2.p + t * l2.u
    det_m = l1.u[1] * l2.u[2] - l1.u[2] * l2.u[1]
    u1_len = norm(l1.u)
    u2_len = norm(l2.u)
    if u1_len == 0 || u2_len == 0 || abs(det_m) <= tol * u1_len * u2_len
        return nothing
    end
    diff = l2.p - l1.p
    s = (diff[1] * l2.u[2] - diff[2] * l2.u[1]) / det_m
    return Point{2,Float64}(l1.p[1] + s * l1.u[1], l1.p[2] + s * l1.u[2])
end

"""
    intersect_lines(l1::Line{2}, l2::Line{2}; tol::Real = 1e-12)

Alias for `intersect(l1, l2)`.
"""
const intersect_lines = intersect

"""
    intersect(h1::Halfplane, h2::Halfplane; tol::Real = 1e-12)

Intersection point of the boundary lines of two halfplanes `h1` and `h2`.
"""
Base.intersect(h1::Halfplane, h2::Halfplane; tol::Real = 1e-12) =
    intersect(h1.boundary, h2.boundary; tol = tol)

"""
    pairwise_intersections(a, b; tol::Real = 1e-12)

Compute all intersection points between two geometric objects `a` and `b`.
Returns a `Vector{Point{2, Float64}}`.
"""
function pairwise_intersections(l1::Line{2}, l2::Line{2}; tol::Real = 1e-12)
    pt = intersect(l1, l2; tol = tol)
    return pt === nothing ? Point{2,Float64}[] : [pt]
end

pairwise_intersections(h1::Halfplane, h2::Halfplane; tol::Real = 1e-12) =
    pairwise_intersections(h1.boundary, h2.boundary; tol = tol)

function pairwise_intersections(l::Line{2}, c::AbsCurve2D; tol::Real = 1e-12)
    pts_t = intersect_line_curve(Line{2,Float64}(l.p, l.u), c)
    return [Point{2,Float64}(pt) for (pt, _) in pts_t]
end

pairwise_intersections(c::AbsCurve2D, l::Line{2}; tol::Real = 1e-12) =
    pairwise_intersections(l, c; tol = tol)

pairwise_intersections(h::Halfplane, c::AbsCurve2D; tol::Real = 1e-12) =
    pairwise_intersections(h.boundary, c; tol = tol)

pairwise_intersections(c::AbsCurve2D, h::Halfplane; tol::Real = 1e-12) =
    pairwise_intersections(h.boundary, c; tol = tol)

function pairwise_intersections(h1::Hyperbola, h2::Hyperbola; tol::Real = 1e-12)
    return intersect_hyperbolas(h1, h2)
end

function _unique_points(pts::Vector{Point{2,Float64}}; tol::Real = 1e-12)
    isempty(pts) && return pts
    res = Point{2,Float64}[]
    tol_sq = tol^2
    for p in pts
        if !any(q -> dist_sq(p, q) <= tol_sq, res)
            push!(res, p)
        end
    end
    return res
end

"""
    vertices(bundle; tol::Real = 1e-12, unique::Bool = false)

Compute the set of vertices (all pairwise intersection points) of a bundle of lines
or halfplanes (or curves).
All elements in the bundle must be of the same type.
Returns a `Vector{Point{2, Float64}}`.
If `unique=true`, duplicate intersection points within distance `tol` are merged.
"""
function vertices(lines::AbstractVector{<:Line{2}}; tol::Real = 1e-12, unique::Bool = false)
    n = length(lines)
    pts = Point{2,Float64}[]
    n < 2 && return pts
    for i = 1:n
        l1 = lines[i]
        for j = (i+1):n
            pt = intersect(l1, lines[j]; tol = tol)
            pt !== nothing && push!(pts, pt)
        end
    end
    return unique ? _unique_points(pts; tol = tol) : pts
end

function vertices(
    halfplanes::AbstractVector{<:Halfplane};
    tol::Real = 1e-12,
    unique::Bool = false,
)
    lines = [h.boundary for h in halfplanes]
    return vertices(lines; tol = tol, unique = unique)
end

function vertices(bundle::AbstractVector; tol::Real = 1e-12, unique::Bool = false)
    isempty(bundle) && return Point{2,Float64}[]
    if all(x -> x isa Line{2}, bundle)
        typed_lines = [Line{2,Float64}(l.p, l.u) for l in bundle]
        return vertices(typed_lines; tol = tol, unique = unique)
    elseif all(x -> x isa Halfplane, bundle)
        typed_hps = [Halfplane2F(h.boundary) for h in bundle]
        return vertices(typed_hps; tol = tol, unique = unique)
    elseif all(x -> x isa AbsCurve2D, bundle)
        pts = Point{2,Float64}[]
        n = length(bundle)
        for i = 1:n
            b1 = bundle[i]
            for j = (i+1):n
                append!(pts, pairwise_intersections(b1, bundle[j]; tol = tol))
            end
        end
        return unique ? _unique_points(pts; tol = tol) : pts
    else
        throw(
            ArgumentError(
                "Bundle elements must all be of the same type (all Line{2} or all Halfplane).",
            ),
        )
    end
end

"""
    BBox(bundle::AbstractVector{<:Line{2}}; tol::Real = 1e-12)
    BBox(bundle::AbstractVector{<:Halfplane}; tol::Real = 1e-12)

Compute the smallest axis-aligned bounding box enclosing all vertices (pairwise intersections)
of the line or halfplane bundle. Returns an uninitialized `BBox` if there are fewer than 2 vertices
or lines are parallel.
"""
function BBox(lines::AbstractVector{<:Line{2}}; tol::Real = 1e-12)
    verts = vertices(lines; tol = tol)
    bb = BBox{2,Float64}()
    isempty(verts) && return bb
    return bound!(bb, verts)
end

function BBox(halfplanes::AbstractVector{<:Halfplane}; tol::Real = 1e-12)
    verts = vertices(halfplanes; tol = tol)
    bb = BBox{2,Float64}()
    isempty(verts) && return bb
    return bound!(bb, verts)
end

"""
    bbox(bundle; kwargs...)

Alias for `BBox(bundle; kwargs...)`.
"""
@inline bbox(bundle; kwargs...) = BBox(bundle; kwargs...)

"""
    default_view(bundle; factor::Real = 1.1, tol::Real = 1e-12)
    default_view(bb::BBox; factor::Real = 1.1)

Compute the default view bounding box for a bundle of lines or halfplanes.
Computes the bounding box of the vertices of the bundle and rescales it by `factor`
(default `1.1`) around its center to provide a comfortable viewing margin.
If all vertices coincide or have zero width/height, padding is added to ensure a non-degenerate 2D view.
"""
function default_view(
    lines::AbstractVector{<:Line{2}};
    factor::Real = 1.1,
    tol::Real = 1e-12,
)
    bb = BBox(lines; tol = tol)
    return default_view(bb; factor = factor)
end

function default_view(
    halfplanes::AbstractVector{<:Halfplane};
    factor::Real = 1.1,
    tol::Real = 1e-12,
)
    bb = BBox(halfplanes; tol = tol)
    return default_view(bb; factor = factor)
end

function default_view(bundle::AbstractVector; factor::Real = 1.1, tol::Real = 1e-12)
    verts = vertices(bundle; tol = tol)
    bb = BBox{2,Float64}()
    if !isempty(verts)
        bound!(bb, verts)
    end
    return default_view(bb; factor = factor)
end

function default_view(bb::BBox{2,T}; factor::Real = 1.1) where {T}
    !bb.f_init && return bb
    res = deepcopy(bb)
    w = width(res)
    h = height(res)
    if w == 0.0 && h == 0.0
        expand_add!(res, 1.0)
    elseif w == 0.0
        res.mini[1] -= 0.5 * h
        res.maxi[1] += 0.5 * h
        expand!(res, factor)
    elseif h == 0.0
        res.mini[2] -= 0.5 * w
        res.maxi[2] += 0.5 * w
        expand!(res, factor)
    else
        expand!(res, factor)
    end
    return res
end

"""
    defaultview(bundle; factor::Real = 1.1)

Alias for `default_view`.
"""
const defaultview = default_view

export vertices, default_view, defaultview, bbox, intersect_lines, pairwise_intersections
