"""
    BasicCompGeometry

A comprehensive library for basic computational geometry operations.
Supports high-dimensional points, segments, lines, point sequences, and axis-aligned bounding boxes.
Built for performance and ease of use using Julia's multiple dispatch system.
"""
module BasicCompGeometry

using Parameters
using StaticArrays
using LinearAlgebra
using DelimitedFiles
using Distributions
using Random

"""
    AbsFMS

Abstract supertype for all finite metric spaces.
Points in a finite metric space are represented by integers 1:n.
"""
abstract type AbsFMS end

include("Points.jl")
include("Segments.jl")
include("PntSeqs.jl")
include("BBoxes.jl")
include("Transforms2D.jl")
include("PntSeqHausdorff.jl")
include("ConvexHull.jl")
include("ConvexHull3D.jl")
include("VirtArray.jl")
include("BBT.jl")
include("WSPD.jl")
include("Diameter.jl")
include("MetricSpace.jl")
include("MVBB.jl")
include("ReadWrite.jl")
include("LongestConvexSubset.jl")

using .VirtArray
using .BBT
using .WSPD
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

export Segment, BBox, BBox2F, Segment2F, Line
export turn_sign,
    is_left_turn, is_right_turn, is_left_eq_turn, is_right_eq_turn, is_collinear
export dist, dist_sq, distance_infty, distance
export exact_diameter, approx_diameter
export Points, centroid, convex_comb, convex_hull
export match_price, cardin, VecPntSeq2F, VecPolygon2F, VecPoint2I
export bottom_left, top_right, width, height, middle, diam, max_dist, is_inside
export VirtArray, BBT, WSPD, MVBB, MetricSpace, ReadWrite, LongestConvexSubset, ConvexHull3D

end
