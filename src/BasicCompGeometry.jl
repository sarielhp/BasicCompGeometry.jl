"""
    BasicCompGeometry

A comprehensive library for basic computational geometry operations.
Supports high-dimensional points, segments, lines, point sequences, and axis-aligned bounding boxes.
Built for performance and ease of use using Julia's multiple dispatch system.
"""
module BasicCompGeometry

using StaticArrays
using LinearAlgebra
using DelimitedFiles
using Printf
using Random

"""
    AbsFMS

Abstract supertype for all finite metric spaces.
Points in a finite metric space are represented by integers 1:n.
"""
abstract type AbsFMS end

"""
    AbsCurve2D{T}

Abstract supertype for 2D curve primitives.
A curve may be bounded (e.g., `Segment{2, T}`, `CircleArc{T}`) or unbounded (e.g., `Hyperbola{T}`).
Concrete subtypes: `Segment{2, T}`, `CircleArc{T}`, `Hyperbola{T}`.
"""
abstract type AbsCurve2D{T} end

include("Kernels.jl")
include("Points.jl")
include("KernelPredicates.jl")
include("Segments.jl")
include("Planes.jl")
include("PntSeqs.jl")
include("BBoxes.jl")
include("Transforms2D.jl")
include("PntSeqHausdorff.jl")
include("ConvexHull.jl")
include("CairoDraw.jl")
include("ConvexHull3D.jl")
include("VirtArray.jl")
include("BBT.jl")
include("WSPD.jl")
include("Diameter.jl")
include("MetricSpace.jl")
include("MVBB.jl")
include("ReadWrite.jl")
include("LongestConvexSubset.jl")
include("PntSeqHierarchy.jl")

using .VirtArray
using .BBT
import .BBT: depth
import .WSPD
using .MVBB
using .MetricSpace
using .ReadWrite
using .LongestConvexSubset
using .ConvexHull3D

# Common high-level operations

"""
    distance_infty(P, Q)

Calculate the maximum (L_∞) distance between corresponding vertices of two point sequences.
Requires `length(P) == length(Q)`.
Useful for assessing the quality of point-to-point matching.
"""
function distance_infty(P::AbsPntSeq{D,T}, Q::AbsPntSeq{D,T}) where {D,T}
    n = length(P)
    @assert n == length(Q)
    n == 0 && return 0.0

    d = 0.0
    for i = 1:n
        d = max(d, dist(P[i], Q[i]))
    end
    return d
end

"""
    distance(y::Point, l::Line)

Calculate the minimum Euclidean distance between query point `y` and infinite line `l`.
"""
function distance(y::Point{D,T}, l::Line{D,T}) where {D,T}
    # projection parameter: t = <y - p, u> / <u, u>
    t = dot(y - l.p, l.u) / dot(l.u, l.u)
    x = l.p + t * l.u
    return dist(x, y)
end

"""
    distance(y::Point, pl::Plane)

Calculate the minimum Euclidean distance between query point `y` and hyperplane `pl`.
"""
function distance(y::Point{D,T}, pl::Plane{D,T}) where {D,T}
    return abs(dot(pl.n, y - pl.p)) / norm(pl.n)
end

"""
    centroid(P)

Calculate the centroid (arithmetic mean) of a collection of points `P`.
Works with point sequences or vectors of points.
"""
centroid(P) = sum(P) / length(P)

"""
    match_price(p_a, p_b, q_a, q_b)

Calculate a cost metric for matching the segment `[p_a, p_b]` to the segment `[q_a, q_b]`.
The price is defined as the product of the sum of the lengths of the two edges 
and the average distance between their corresponding endpoints.
"""
function match_price(
    p_a::Point{N,T},
    p_b::Point{N,T},
    q_a::Point{N,T},
    q_b::Point{N,T},
) where {N,T}
    l_avg_endpoint = (dist(p_a, q_a) + dist(p_b, q_b)) / 2.0
    l_sum_edges = dist(p_a, p_b) + dist(q_a, q_b)
    return l_sum_edges * l_avg_endpoint
end

"""
    VecPoint2I

Type alias for a vector of 2D integer points.
"""
const VecPoint2I = Vector{Point2I}

"""
    VecPntSeq2F

Type alias for a vector of 2D point sequences with Float64 coordinates.
"""
const VecPntSeq2F = Vector{PntSeq2F}

"""
    VecPolygon2F

Alias for backward compatibility.
"""
const VecPolygon2F = VecPntSeq2F

include("Spheres.jl")
include("Hyperbolas.jl")
include("Parabolas.jl")
include("CurvePolygon.jl")
include("Ellipses.jl")
include("Bezier.jl")
include("IpeDraw.jl")
include("Halfplanes.jl")
include("Bundles.jl")
include("Arrangements.jl")
include("MonotoneTriangulation.jl")
include("VisibilityPolygon.jl")
include("MinkowskiSum.jl")
include("RMQ.jl")
include("ChainIndex.jl")
include("LocalIncidence.jl")

# Core types and distinctive operations are exported. Generic geometry helpers
# remain available through qualified calls, explicit imports, or the opt-in
# Extras submodule; see the migration notes for v0.5.0 in CHANGELOG.md.
export Segment, BBox, BBox2F, Segment2F, Line, Plane, Plane2F
export Halfplane, Halfplane2F, boundary_line, rand_halfplane, in_interior, on_boundary
export default_view,
    bbox, union_bbox, intersect_lines, intersections, pairwise_intersections
export VerticalTrapezoid,
    VerticalTrapezoid2F,
    ArrVertDecomp,
    ArrVertDecomp2F,
    vertical_decomposition,
    trapezoids,
    view_box
export trapezoid_chains
export Sphere,
    Circle, Circle2F, Circle2I, Sphere2F, Sphere3F, CircleArc, CircleArc2F, invert
export Ellipse, Ellipse2F, EllipticArc, EllipticArc2F
export CubicBezier,
    CubicBezier2F,
    CubicBezier3F,
    CubicSpline,
    subdivide,
    flatten,
    interpolate_catmull_rom,
    interpolate_natural_spline
export Curve2D,
    Curve2DF,
    CurvePolygon2D,
    CurvePolygon2DF,
    AbsCurve2D,
    point_on,
    tangent_line,
    intersect_line_curve
export Hyperbola, Parabola, ParabolaArc, ParabolaArc2F, intersect_hyperbolas, bisector
export turn_sign,
    is_left_turn, is_right_turn, is_left_eq_turn, is_right_eq_turn, is_collinear
export dist, dist_sq
export exact_diameter, approx_diameter, approx_diameter_subspace, exact_diameter_subspace
export centroid, convex_comb, convex_hull
export orient,
    orient_sign,
    on_segment,
    proper_crossing,
    midpoint,
    affine_interpolate,
    segment_at_x,
    compare_angle,
    angle_less,
    segment_key,
    id_pair,
    twice_area,
    polygon_segments
export polar, polar_line
export bottom_left, corner, top_right, diam, max_dist, is_inside
export BBT, WSPD, MVBB, MetricSpace, LongestConvexSubset, ConvexHull3D, IpeDraw
export PolygonTriangulation,
    MonotoneTriangulation,
    triangulate_polygon,
    triangulate_monotone,
    triangulate_seidel,
    is_x_monotone,
    is_y_monotone,
    monotone_chains,
    monotone_chain_indices,
    merge_chains,
    stack_diagonals,
    polygon_diagonals,
    visibility_polygon,
    polygon_kernel,
    is_star_shaped,
    is_convex,
    minkowski_sum,
    minkowski_sum_decomposition,
    build_rmq,
    range_argmin,
    build_chain_index,
    chain_window,
    nedges,
    vertical_ray_segment,
    first_positive_vertical_contact,
    NoHit,
    PointHit,
    OverlapHit,
    open_sector_contains,
    incoming_sector_contains

"""
    BasicCompGeometry.Extras

Additional supported public API names omitted from the default export list.
Exports 91 helpers, accessors, compatibility aliases, auxiliary namespaces,
and Cairo/LaTeX utilities. All bindings refer to the original objects in
`BasicCompGeometry`; visualization implementations still require their
optional packages to be loaded.

Load both sets of public exports with:
```julia
import BasicCompGeometry as CG
using .CG, .CG.Extras
```
The explicit export list contains only supported public API names.
"""
module Extras

using ..BasicCompGeometry:
    AbsPolygon,
    Canvas,
    MatPolygon,
    Points,
    Polygon_random_sphere,
    ReadWrite,
    VecPntSeq2F,
    VecPoint2I,
    VecPolygon2F,
    VertTrapezoid,
    VirtArray,
    a_coeff,
    add_latex_macros!,
    add_latex_packages!,
    after_tangent_to_polygon,
    alpha1,
    alpha2,
    append_latex_preamble!,
    area,
    at,
    axis_direction,
    b_coeff,
    before_tangent_to_polygon,
    bisection_point,
    bound!,
    boundary,
    c_coeff,
    cairo_draw_halfplane,
    cairo_draw_halfplanes,
    cairo_draw_latex,
    cairo_draw_latex_page,
    cairo_draw_points,
    cairo_draw_polygon,
    cairo_draw_setup,
    cairo_set_line_width,
    cardin,
    center,
    clip,
    compare_right_arc_to_point,
    compare_right_arcs,
    complement,
    convex_coef,
    d_max,
    d_min,
    defaultview,
    depth,
    derivative,
    description,
    direction,
    dist_point_segment,
    dist_segment_segment,
    dist_subspace,
    distance,
    distance_infty,
    end_point,
    expand!,
    expand_add!,
    faces,
    get_file_path,
    get_latex_preamble,
    hausdorff_dist_subseg,
    height,
    init!,
    latex_to_pdf,
    lines,
    locate,
    match_price,
    max_dist_subspace,
    middle,
    open_canvas,
    p_coeff,
    pdf_merge,
    pnt_seq_random_sphere,
    prefix_lengths,
    push_distinct!,
    r_major,
    r_minor,
    rand_polygon,
    random_halfplane,
    read_halfplanes,
    read_latex_snippet,
    reset_latex_preamble!,
    rotation_angle,
    set_latex_preamble!,
    start_point,
    split_edges,
    vertex,
    vertices,
    width,
    write_halfplanes,
    write_plt

export AbsPolygon,
    Canvas,
    MatPolygon,
    Points,
    Polygon_random_sphere,
    ReadWrite,
    VecPntSeq2F,
    VecPoint2I,
    VecPolygon2F,
    VertTrapezoid,
    VirtArray,
    a_coeff,
    add_latex_macros!,
    add_latex_packages!,
    after_tangent_to_polygon,
    alpha1,
    alpha2,
    append_latex_preamble!,
    area,
    at,
    axis_direction,
    b_coeff,
    before_tangent_to_polygon,
    bisection_point,
    bound!,
    boundary,
    c_coeff,
    cairo_draw_halfplane,
    cairo_draw_halfplanes,
    cairo_draw_latex,
    cairo_draw_latex_page,
    cairo_draw_points,
    cairo_draw_polygon,
    cairo_draw_setup,
    cairo_set_line_width,
    cardin,
    center,
    clip,
    compare_right_arc_to_point,
    compare_right_arcs,
    complement,
    convex_coef,
    d_max,
    d_min,
    defaultview,
    depth,
    derivative,
    description,
    direction,
    dist_point_segment,
    dist_segment_segment,
    dist_subspace,
    distance,
    distance_infty,
    end_point,
    expand!,
    expand_add!,
    faces,
    get_file_path,
    get_latex_preamble,
    hausdorff_dist_subseg,
    height,
    init!,
    latex_to_pdf,
    lines,
    locate,
    match_price,
    max_dist_subspace,
    middle,
    open_canvas,
    p_coeff,
    pdf_merge,
    pnt_seq_random_sphere,
    prefix_lengths,
    push_distinct!,
    r_major,
    r_minor,
    rand_polygon,
    random_halfplane,
    read_halfplanes,
    read_latex_snippet,
    reset_latex_preamble!,
    rotation_angle,
    set_latex_preamble!,
    start_point,
    split_edges,
    vertex,
    vertices,
    width,
    write_halfplanes,
    write_plt

end # module Extras

end
