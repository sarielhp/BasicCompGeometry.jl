"""
    MetricSpace

The `MetricSpace` module provides a framework for working with abstract finite metric spaces (FMS).
It is designed to decouple algorithms that operate on metric spaces (like clustering or greedy 
permutations) from the underlying representation of the points and the distance function.

# Core Concepts

- **Abstract Finite Metric Space (`AbsFMS`)**: The base type for all metric spaces in this module.
- **Indices as Points**: Points in a finite metric space are represented by integers `1:n`.
- **Interface**: Any subtype of `AbsFMS` must implement `Base.size` and `dist`.

# Supported Space Types

1. `PointsSpace`: Defined by a collection of points of any type that supports `dist(p1, p2)`.
2. `MPointsSpace`: Defined by a matrix where columns represent points. Uses Euclidean distance.
3. `PermutMetric`: A wrapper that allows reordering or subsetting an existing metric space via a permutation vector.
4. `SpherePMetric`: A specialized space where distances are defined as angles relative to a base point (useful for spherical geometry or directional data).

# Algorithms

- `greedy_permutation_naive`: Implements Gonzalez's algorithm for k-center clustering and generating 
  hierarchical orderings of points.
- `greedy_permutation_vanity`: A variation of the greedy permutation that incorporates a "vanity" 
  metric to choose between points that are equidistant from existing centers.
"""
module MetricSpace

using ..BasicCompGeometry
using Distances
using LinearAlgebra

import ..BasicCompGeometry: dist, AbsFMS

"""
    metric(space::AbsFMS, x::Int, y::Int)

Computes the distance between points `x` and `y` in the given `space`. 
This is a convenience wrapper around `dist(space, x, y)`.
"""
function metric(space::AbsFMS, x, y)
    return dist(space, x, y)
end

"""
    original(m::AbsFMS, i::Int)

Return the original index of the point at position `i`. For simple spaces, this is 
just `i`. For spaces that involve permutations (like `PermutMetric`), this 
returns the index in the underlying source space.
"""
function original(m::AbsFMS, i::Int)
    return i
end

"""
    size(space::AbsFMS)

Return the number of points in the finite metric space.
"""
function Base.size(space::AbsFMS)
    error("size not implemented for $(typeof(space))")
end

####################################################################
# PointsSpace: Metric space defined by a vector of points.
####################################################################

"""
    PointsSpace{PType} <: AbsFMS

A metric space defined by a vector of points of type `PType`.
The distance between points `i` and `j` is computed using `dist(P[i], P[j])`.

# Fields
- `P::Vector{PType}`: The collection of points.
"""
struct PointsSpace{PType} <: AbsFMS
    P::Vector{PType}
end

"""
    dist(P::PointsSpace, i::Int, j::Int)

Returns the distance between the `i`-th and `j`-th point in the vector.
"""
function dist(P::PointsSpace{PType}, x, y)::Float64 where {PType}
    if x == y
        return 0.0
    end
    return dist(P.P[x], P.P[y])
end

"""
    dist_real(P::PointsSpace, i::Int, y_real)

Returns the distance between the `i`-th point in the space and an external point `y_real`.
"""
function dist_real(P::PointsSpace{PType}, x, y_real)::Float64 where {PType}
    return dist(P.P[x], y_real)
end

function Base.size(P::PointsSpace)
    return length(P.P)
end

####################################################################
# MPointsSpace: Metric space defined by a matrix.
# Columns are the points.
####################################################################

"""
    MPointsSpace{PType} <: AbsFMS

A metric space where points are columns of a matrix. 
Distances are computed using the Euclidean metric.

# Fields
- `m::Matrix{PType}`: Matrix where each column is a point.
- `n::Int`: Number of points (columns) to consider.
"""
struct MPointsSpace{PType} <: AbsFMS
    m::Matrix{PType}
    n::Int     # Number of points

    function MPointsSpace(m::Matrix{PType}, n::Int = 0) where {PType}
        if n == 0
            n = size(m, 2)
        end
        @assert(0 < n <= size(m, 2))
        return new{PType}(m, n)
    end
end

"""
    dist(P::MPointsSpace, i::Int, j::Int)

Euclidean distance between column `i` and column `j`.
"""
function dist(P::MPointsSpace, x, y)::Float64
    if x == y
        return 0.0
    end
    return euclidean(@view(P.m[:, x]), @view(P.m[:, y]))
end

"""
    dist_real(P::MPointsSpace, i::Int, y_real)

Euclidean distance between column `i` and an external vector `y_real`.
"""
function dist_real(P::MPointsSpace, x, y_real)
    return euclidean(@view(P.m[:, x]), y_real)
end

function Base.size(P::MPointsSpace)
    return P.n
end

####################################################################
# PermutMetric: A permutation of a finite metric space.
####################################################################

"""
    PermutMetric{MetricType} <: AbsFMS

A decorator for an existing metric space that applies a permutation (or subsetting) 
to the point indices. This is useful for algorithms that reorder points.

# Fields
- `n::Int`: Number of points in the permuted space.
- `m::MetricType`: The underlying metric space.
- `I::Vector{Int}`: Mapping from the new index to the original index.
"""
struct PermutMetric{MetricType} <: AbsFMS
    n::Int  # Number of points
    m::MetricType
    I::Vector{Int}
end

function original(m::PermutMetric, i::Int)
    return m.I[i]
end

"""
    PermutMetric(m::AbsFMS)

Construct an identity permutation of the given metric space.
"""
function PermutMetric(m::MetricType) where {MetricType}
    n = size(m)
    return PermutMetric(n, m, collect(1:n))
end

"""
    PermutMetric(m::AbsFMS, I::Vector{Int})

Construct a permutation of `m` using the mapping `I`.
"""
function PermutMetric(m::MetricType, I::Vector{Int}) where {MetricType}
    return PermutMetric(size(m), m, copy(I))
end

function dist(P::PermutMetric, x::Int, y::Int)::Float64
    if x == y
        return 0.0
    end
    return dist(P.m, P.I[x], P.I[y])
end

function dist_real(P::PermutMetric, x, y_real)
    return dist_real(P.m, P.I[x], y_real)
end

function Base.size(P::PermutMetric)
    return P.n
end

"""
    swap!(P::PermutMetric, i::Int, j::Int)

Swap the points at positions `i` and `j` in the permutation vector.
"""
function swap!(P::PermutMetric, x, y)
    if x != y
        P.I[x], P.I[y] = P.I[y], P.I[x]
    end
    return nothing
end

####################################################################
# SpherePMetric: A permutation of a finite metric space, where the
#    distance is angular between the points.
####################################################################

"""
    SpherePMetric{MetricType} <: AbsFMS

A specialized metric space where distances are defined as the angle formed 
at a base point `b_point` between two points `i` and `j`.

Given base point B and points X, Y, the distance is the angle ∠XBY calculated 
via the Law of Cosines using the distances in the underlying metric space `m`.

# Fields
- `b_point::Int`: Index of the base point in the underlying space `m`.
- `m::MetricType`: The underlying metric space.
- `I::Vector{Int}`: Indices of points included in this spherical space.
- `D::Vector{Float64}`: Precomputed distances from each point in `I` to `b_point`.
"""
struct SpherePMetric{MetricType} <: AbsFMS
    b_point::Int
    m::MetricType
    I::Vector{Int}
    D::Vector{Float64}   # The distances of the points to the base point
end

function SpherePMetric(m::MetricType, b_point::Int) where {MetricType}
    return SpherePMetric(b_point, m, Int[], Float64[])
end

"""
    push!(P::SpherePMetric, pnt::Int, d::Float64)

Add a point with index `pnt` (from the underlying space) and its distance `d` 
from the base point to the spherical metric space.
"""
function Base.push!(P::SpherePMetric, pnt::Int, d::Float64)
    push!(P.I, pnt)
    push!(P.D, d)
    return P
end

function dist(P::SpherePMetric, x::Int, y::Int)::Float64
    if x == y
        return 0.0
    end

    ℓ_b_x = P.D[x]
    ℓ_b_y = P.D[y]
    ℓ_x_y = dist(P.m, P.I[x], P.I[y])

    return calculate_angle_facing_a(ℓ_x_y, ℓ_b_x, ℓ_b_y)
end

function Base.size(P::SpherePMetric)
    return length(P.I)
end

function swap!(P::SpherePMetric, x, y)
    if x != y
        P.I[x], P.I[y] = P.I[y], P.I[x]
        P.D[x], P.D[y] = P.D[y], P.D[x]
    end
    return nothing
end

######################################################################

function get_points_assigned_to_pos(M, I, D, N, pos, sz)
    J = Int[]
    x = I[pos]
    for i ∈ 1:sz
        y = I[i]
        d = dist(M, x, y)
        if d < D[i] || D[i] < 0.0
            push!(J, i)
        end
    end
    return J
end

function update_distances(M::AbsFMS, I, D, N, pos, sz)
    x = I[pos]
    for i ∈ (pos+1):sz
        y = I[i]
        d = dist(M, x, y)
        if d < D[i] || D[i] < 0.0
            D[i] = d
            N[i] = pos
        end
    end
end

"""
    greedy_permutation_naive(M::AbsFMS, n::Int)

Compute a greedy permutation of a metric space `M` using Gonzalez's algorithm.
The algorithm selects `n` points such that each selected point maximizes its 
minimum distance to the set of previously selected points.

# Returns
- `I`: Indices of the points in the greedy permutation.
- `D`: Distance of each point to its nearest selected center at the time it was selected.
- `N`: Index of the nearest center for each point.
"""
function greedy_permutation_naive(M::AbsFMS, n::Int)
    sz = size(M)
    I = collect(1:sz)
    N = fill(1, sz)
    D = fill(-1.0, sz)

    update_distances(M, I, D, N, 1, sz)

    for i ∈ 2:n
        pos = argmax(view(D, i:sz)) + i - 1
        D[i-1] = D[pos]
        I[i], I[pos] = I[pos], I[i]
        D[i], D[pos] = D[pos], D[i]
        update_distances(M, I, D, N, i, sz)
    end

    if n == sz
        D[sz] = 0.0
    end

    return I, D, N
end

"""
    greedy_permutation_vanity(M::AbsFMS, vanity::Vector{Float64}, n::Int)

A variant of the greedy permutation where ties (or near-ties) in distance 
are broken by minimizing a "vanity" score. 

The algorithm first identifies the point furthest from current centers, 
finds all points currently "assigned" to that furthest point (i.e., would 
be moved to the center set if it were chosen), and then picks the point 
within that set that has the minimum vanity score.
"""
function greedy_permutation_vanity(M::AbsFMS, vanity::Vector{Float64}, n::Int)
    sz = size(M)
    I = collect(1:sz)
    N = fill(1, sz)
    D = fill(-1.0, sz)

    update_distances(M, I, D, N, 1, sz)

    for i ∈ 2:n
        tpos = argmax(view(D, i:sz)) + i - 1
        J = get_points_assigned_to_pos(M, I, D, N, tpos, sz)
        j_min = argmin(vanity[j] for j ∈ J)
        pos = J[j_min]

        D[i-1] = D[pos]
        I[i], I[pos] = I[pos], I[i]
        D[i], D[pos] = D[pos], D[i]
        vanity[i], vanity[pos] = vanity[pos], vanity[i]
    end

    if n == sz
        D[sz] = 0.0
    end

    return I, D, N
end

"""
    calculate_angle_facing_a(a, b, c)

Calculates the angle (in radians) opposite side `a` in a triangle with sides `a`, `b`, and `c` 
using the Law of Cosines: `a² = b² + c² - 2bc * cos(α)`.
"""
function calculate_angle_facing_a(a::Real, b::Real, c::Real)::Float64
    if a < 0 || b < 0 || c < 0
        throw(ArgumentError("Edge lengths must be non-negative."))
    end

    if iszero(a)
        return 0.0
    end

    if iszero(b) || iszero(c)
        return π / 2.0
    end

    numerator = b^2 + c^2 - a^2
    denominator = 2 * b * c
    cos_alpha = clamp(numerator / denominator, -1.0, 1.0)
    return acos(cos_alpha)
end

export AbsFMS, PointsSpace, MPointsSpace, PermutMetric, SpherePMetric
export dist, dist_real, original, swap!, metric
export greedy_permutation_naive, greedy_permutation_vanity

end # module MetricSpace
