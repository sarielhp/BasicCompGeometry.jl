###############################################
### PntSeq type

"""
    AbsPntSeq{D, T}

Abstract supertype for all representations of a sequence of points (polygonal 
chain or polygon) in D dimensions with coordinate type T.
Subtypes must implement the `AbsFMS` interface (`dist` and `size`) and 
provide access to vertices.
"""
abstract type AbsPntSeq{D,T} <: AbsFMS end

"""
    PntSeq{D, T}

A wrapper around `Vector{Point{D, T}}` representing a sequence of points, 
a polygonal curve, or a classical polygon in D dimensions.
The points are stored sequentially in the `pnts` field.

`PntSeq` is a subtype of `AbsPntSeq{D, T}`, meaning it can be treated as a finite metric space
where its vertices are the points and the distance is Euclidean.
"""
struct PntSeq{D,T} <: AbsPntSeq{D,T}
    pnts::Vector{Point{D,T}}
end

"""
    MatPntSeq{D, T, M, V}

A representation of a sequence of points backed by an `AbstractMatrix{T}` 
of size `D x N`. It provides a zero-copy view of the matrix columns as 
`Point{D, T}` objects.
"""
struct MatPntSeq{D,T,M<:AbstractMatrix{T},V<:AbstractVector{Point{D,T}}} <: AbsPntSeq{D,T}
    mat::M
    pnts::V

    function MatPntSeq(mat::M) where {T,M<:AbstractMatrix{T}}
        D = size(mat, 1)
        # Use reinterpret and vec to get a view of columns as Points
        # For Array, vec() is zero-copy.
        pnts = reinterpret(Point{D,T}, vec(mat))
        return new{D,T,M,typeof(pnts)}(mat, pnts)
    end
end

"""
    dim(P::AbsPntSeq)

Return the dimension of the points in the point sequence.
"""
@inline dim(P::AbsPntSeq{D,T}) where {D,T} = D

"""
    dist(P::AbsPntSeq, i::Int, j::Int)

Distance between vertex `i` and vertex `j` of the point sequence.
Part of the `AbsFMS` interface.
"""
@inline function dist(P::AbsPntSeq, i::Int, j::Int)
    if i == j
        return 0.0
    end
    return dist(P.pnts[i], P.pnts[j])
end

# Constructors
"""
    PntSeq{D, T}()

Construct an empty point sequence for dimension `D` and type `T`.
"""
PntSeq{D,T}() where {D,T} = PntSeq{D,T}(Point{D,T}[])

# Note: Struct definition already provides: PntSeq(pnts::Vector{Point{D, T}})

# Accessors
"""
    Points(P::AbsPntSeq)

Return the underlying vector of points of the point sequence.
"""
@inline Points(P::AbsPntSeq) = P.pnts

"""
    cardin(P)

Return the number of vertices (cardinality) in the point sequence `P`.
"""
@inline cardin(P::AbsPntSeq) = length(P.pnts)

"""
    translate!(P, v)

Translate the point sequence `P` in-place by the vector `v`.
"""
function translate!(P::AbsPntSeq{D,T}, v::Point{D,T}) where {D,T}
    P.pnts .+= Ref(v)
    return P
end

"""
    scale!(P, s)

Scale the point sequence `P` in-place by factor `s` relative to the origin.
"""
function scale!(P::AbsPntSeq{D,T}, s::Real) where {D,T}
    P.pnts .*= s
    return P
end

"""
    geom_length(pnt_seq)

Calculate the total Euclidean length of the polygonal curve (sum of edge lengths).
"""
function geom_length(pnt_seq::AbsPntSeq{D,T}) where {D,T}
    len = 0.0
    pnts = pnt_seq.pnts
    for i = 1:(length(pnts)-1)
        len += dist(pnts[i], pnts[i+1])
    end
    return len
end

"""
    prefix_lengths(pnt_seq)

Return a vector containing the cumulative arc-length distances at each vertex of the point sequence.
The first element is always 0.0, and the last is the total geometric length.
"""
function prefix_lengths(pnt_seq::AbsPntSeq{D,T}) where {D,T}
    n = length(pnt_seq)
    v = Vector{Float64}(undef, n)
    n == 0 && return v

    v[1] = 0.0
    pnts = pnt_seq.pnts
    for i = 2:n
        v[i] = v[i-1] + dist(pnts[i-1], pnts[i])
    end
    return v
end

"""
    at(pnt_seq, t)

Return a point on the polygonal curve at the normalized parameter `t` in `[0, 1]`.
`t=0` is the first vertex, `t=1` is the last vertex.
The parameterization is based on arc-length (linear interpolation along edges).
"""
function at(pnt_seq::AbsPntSeq{D,T}, t::Real) where {D,T}
    n = length(pnt_seq)
    n == 0 && throw(ArgumentError("Empty point sequence"))
    n == 1 && return pnt_seq.pnts[1]

    plens = prefix_lengths(pnt_seq)
    total_len = plens[end]
    total_len == 0.0 && return pnt_seq.pnts[1]

    # Target absolute arc-length
    target = t * total_len
    if target <= 0
        return pnt_seq.pnts[1]
    elseif target >= total_len
        return pnt_seq.pnts[n]
    end

    # Find the segment containing the target distance
    i = searchsortedfirst(plens, target)
    i <= 1 && return pnt_seq.pnts[1]

    segment_len = plens[i] - plens[i-1]
    local_t = (target - plens[i-1]) / segment_len
    return convex_comb(pnt_seq.pnts[i-1], pnt_seq.pnts[i], local_t)
end

# Base Extensions for collection-like behavior
Base.length(P::AbsPntSeq) = length(P.pnts)
Base.getindex(P::AbsPntSeq, i::Int) = P.pnts[i]
Base.setindex!(P::AbsPntSeq, v, i::Int) = (P.pnts[i] = v)
Base.firstindex(P::AbsPntSeq) = 1
Base.lastindex(P::AbsPntSeq) = length(P.pnts)
Base.iterate(P::AbsPntSeq) = iterate(P.pnts)
Base.iterate(P::AbsPntSeq, state) = iterate(P.pnts, state)
Base.push!(P::PntSeq, p::Point) = push!(P.pnts, p)
Base.pop!(P::PntSeq) = pop!(P.pnts)
Base.empty!(P::PntSeq) = empty!(P.pnts)
Base.first(P::AbsPntSeq) = first(P.pnts)
Base.last(P::AbsPntSeq) = last(P.pnts)
# For AbsFMS interface, size should return the number of points as an Int
Base.size(P::AbsPntSeq) = length(P.pnts)
Base.size(P::AbsPntSeq, dim::Int) = dim == 1 ? length(P.pnts) : 1

function Base.show(io::IO, pnt_seq::PntSeq{D,T}) where {D,T}
    print(io, "PntSeq{$D,$T} with $(length(pnt_seq)) points")
end

function Base.show(io::IO, pnt_seq::MatPntSeq{D,T}) where {D,T}
    print(io, "MatPntSeq{$D,$T} with $(length(pnt_seq)) points")
end

"""
    Matrix(P::AbsPntSeq)

Convert the point sequence's vertices into a `D x N` matrix where `N` is the number of points.
"""
function Base.Matrix(P::AbsPntSeq{D,T}) where {D,T}
    m = Matrix{T}(undef, D, length(P))
    for (i, p) in enumerate(P.pnts)
        m[:, i] = p
    end
    return m
end

"""
    write_plt(filename, P)

Write the point sequence coordinates to a PLT-compatible delimited file (N rows, D columns).
"""
function write_plt(filename::String, P::AbsPntSeq)
    open(filename, "w") do io
        writedlm(io, Matrix(P)', ',')
    end
end

"""
    simplify(P, r)

Perform distance-based simplification of the point sequence `P`.
Vertices are skipped if they are within distance `r` of the last kept vertex.
Returns `(simplified_pnt_seq, kept_indices)`.
"""
function simplify(P::AbsPntSeq{D,T}, r::Real) where {D,T}
    n = length(P)
    n <= 2 && return P, collect(1:n)

    pout = PntSeq{D,T}()
    pindices = Int[]

    push!(pout, P[1])
    push!(pindices, 1)

    curr = P[1]
    for i = 2:(n-1)
        if dist(P[i], curr) > r
            push!(pout, P[i])
            push!(pindices, i)
            curr = P[i]
        end
    end

    push!(pout, P[n])
    push!(pindices, n)
    return pout, pindices
end

"""
    sample_uniformly(P, n)

Resample the point sequence `P` into `n` vertices that are evenly spaced by arc-length.
"""
function sample_uniformly(P::AbsPntSeq{D,T}, n::Int) where {D,T}
    n <= 1 && throw(ArgumentError("n must be > 1"))
    new_P = PntSeq{D,T}()
    for i = 1:n
        push!(new_P, at(P, (i-1)/(n-1)))
    end
    return new_P
end

"""
    rand_pnt_seq(D, T, n)

Generate a random point sequence in D dimensions with `n` vertices, each coordinate in `[0, 1)`.
"""
function rand_pnt_seq(D::Int, T::Type{<:Real}, n::Int)
    return PntSeq([rand_point(D, T) for _ = 1:n])
end

"""
    pnt_seq_random_sphere(D, T, n)

Generate a random point sequence in D dimensions with `n` vertices sampled uniformly from the unit sphere.
"""
function pnt_seq_random_sphere(D::Int, T::Type{<:Real}, n::Int)
    pnts = Vector{Point{D,T}}(undef, n)
    for i = 1:n
        p = randn(T, D)
        mag = norm(p)
        if mag > 0
            p ./= mag
        end
        pnts[i] = Point{D,T}(p)
    end
    return PntSeq(pnts)
end

"""
    PntSeq2I

Alias for a 2D point sequence with `Int64` coordinates.
"""
const PntSeq2I = PntSeq{2,Int64}

"""
    PntSeq2F

Alias for a 2D point sequence with `Float64` coordinates.
"""
const PntSeq2F = PntSeq{2,Float64}

"""
    PntSeq3F

Alias for a 3D point sequence with `Float64` coordinates.
"""
const PntSeq3F = PntSeq{3,Float64}

# Aliases for backward compatibility
"""
    AbsPolygon

Alias for `AbsPntSeq`.
"""
const AbsPolygon = AbsPntSeq

"""
    Polygon

Alias for `PntSeq`.
"""
const Polygon = PntSeq

"""
    MatPolygon

Alias for `MatPntSeq`.
"""
const MatPolygon = MatPntSeq

const Polygon2I = PntSeq2I
const Polygon2F = PntSeq2F
const Polygon3F = PntSeq3F

"""
    rand_polygon(D, T, n)

Alias for `rand_pnt_seq`.
"""
const rand_polygon = rand_pnt_seq

"""
    Polygon_random_sphere(D, T, n)

Alias for `pnt_seq_random_sphere`.
"""
const Polygon_random_sphere = pnt_seq_random_sphere

export PntSeq, AbsPntSeq, MatPntSeq, PntSeq2I, PntSeq2F, PntSeq3F
export AbsPolygon, Polygon, MatPolygon, Polygon2I, Polygon2F, Polygon3F
export Points, cardin, dim, geom_length, prefix_lengths
export translate!, scale!, simplify, sample_uniformly, at
export rand_pnt_seq, pnt_seq_random_sphere, write_plt
export rand_polygon, Polygon_random_sphere
