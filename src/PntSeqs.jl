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
    GPntSeq{D, T, P}

A sequence of points of type `P` in D dimensions with coordinate type T: a polygonal
curve or a classical polygon. The points are stored sequentially in the `pnts` field.
`P` must be a concrete `AbstractVector{T}` type; `PntSeq{D, T}` is the case
`P = Point{D, T}`. `P` may also be a virtual point type ([`AbsVirtualPoint`](@ref), such as
[`VPoint`](@ref)), points known only through distances: the sequence is then a
`GPntSeq{AnyDim, Float64, P}`, and functions that need coordinates do not accept it.

Choosing the point type: static points (`Point{D, T}`) have their loops unrolled and
compiled for each D, and creating one costs nothing. Points whose length is known
only at run time (`Vector{T}`) are compiled once for all dimensions and run plain
vectorized loops, but each new point is a heap allocation. Measured on distance
kernels (FrechetDist's `example_higher/point_types.jl`, AMD Ryzen 7 5700G, Julia
1.13): static points are faster up to D = 16, the two are close for D = 32 to 96,
and from D = 128 on `Vector` points compute distances 2 to 4.5 times faster. Static
code also takes seconds to compile at D = 1000, and minutes for a large algorithm at
D = 6000. Creating points stays faster with static points up to about D = 512.

`D` is either the dimension or [`AnyDim`](@ref), for a dimension known only at run
time (`dim(P)` is then the length of the first point). Points whose length is not part
of their type are checked on construction and by `push!`: each must have D
coordinates, or, for `AnyDim`, as many as the first point. Static points need no
check.

`GPntSeq` is a subtype of `AbsPntSeq{D, T}`, meaning it can be treated as a finite
metric space where its vertices are the points and the distance is Euclidean.
"""
struct GPntSeq{D,T,P<:Union{AbstractVector{T},AbsVirtualPoint}} <: AbsPntSeq{D,T}
    pnts::Vector{P}

    # The default constructor, plus the length check (a no-op for static points).
    function GPntSeq{D,T,P}(pnts) where {D,T,P<:Union{AbstractVector{T},AbsVirtualPoint}}
        v = convert(Vector{P}, pnts)
        check_point_lengths(D, T, v)
        return new{D,T,P}(v)
    end
end

# Every point of `pnts` has D coordinates (for AnyDim: as many as the first one).
# Static points carry their length in their type: nothing to check.
check_point_lengths(_, _, ::Vector{<:StaticVector}) = nothing
function check_point_lengths(D, _, pnts::Vector{<:AbstractVector})
    isempty(pnts) && return nothing
    n = D === AnyDim ? length(first(pnts)) : D
    for (j, p) in enumerate(pnts)
        length(p) == n || throw_point_length(j, length(p), n)
    end
    return nothing
end

# Virtual points have no coordinates: a sequence of them has D = AnyDim and Float64
# distances. The check is on types only.
function check_point_lengths(D, T, ::Vector{<:AbsVirtualPoint})
    D === AnyDim && T === Float64 && return nothing
    throw(
        ArgumentError(
            "a sequence of virtual points is a GPntSeq{AnyDim,Float64,P}, not {$D,$T}",
        ),
    )
end

@noinline throw_point_length(j, len, n) =
    throw(DimensionMismatch("point $j has $len coordinates, expected $n"))

"""
    PntSeq{D, T}

A sequence of `Point{D, T}` (static vectors), stored contiguously in the `pnts`
field: a polygonal curve or a classical polygon in D dimensions. Alias for
`GPntSeq{D, T, Point{D, T}}`.
"""
const PntSeq{D,T} = GPntSeq{D,T,Point{D,T}}

PntSeq(pnts::Vector{Point{D,T}}) where {D,T} = PntSeq{D,T}(pnts)

"""
    MatPntSeq(mat::AbstractMatrix{T}; copy::Bool = false)
    MatPntSeq{D, T, M, V}

A sequence of points stored in a `D × N` matrix: column `j` is point `j`, read as a
`Point{D, T}`. `D = size(mat, 1)`.

**By default the sequence shares the matrix's memory; nothing is copied.**

- Changes to the matrix are seen by the sequence: after `mat[1, 2] = 5.0`,
  `P[2][1] == 5.0`.
- Changes through the sequence are written into the matrix: `P[j] = p`,
  [`translate!`](@ref) and [`scale!`](@ref) modify `mat`.
- The same holds for a view or an adjoint: the sequence then shares the memory of
  the view's parent array.
- Points read from the sequence are copies (`P[j]` is a static vector built from
  column `j`), so they cannot be written through; changing a point later does not
  change the sequence.
- The number of points is fixed: `push!`, `pop!` and `empty!` are not defined.

**`copy = true` copies the matrix** into a new dense `Matrix{T}` first. The sequence
is then independent: later changes to `mat` are not seen, and changes through the
sequence do not reach `mat`. Use it when the matrix will be modified afterwards, or
when the sequence will be translated or scaled and the matrix must stay as it was.
`Matrix(P)` always returns a copy.

Since `D` is read from the matrix at run time, this constructor is not type-stable;
pass the result to a function (a function barrier) before using it in a hot loop.
"""
struct MatPntSeq{D,T,M<:AbstractMatrix{T},V<:AbstractVector{Point{D,T}}} <: AbsPntSeq{D,T}
    mat::M
    pnts::V

    function MatPntSeq(mat::AbstractMatrix{T}; copy::Bool = false) where {T}
        # `copy` names the option, so this calls the Matrix constructor, which
        # always allocates a new array.
        stored = copy ? Matrix{T}(mat) : mat
        D = size(stored, 1)
        # reinterpret views the column-major coordinates as consecutive points;
        # vec() of an Array is zero-copy, so the matrix memory is shared.
        #
        # Resizing is deliberately not offered. A vector that shares a matrix's
        # memory (vec(mat)) behaves inconsistently when resized: push! copies it to
        # new memory, but on Julia 1.13 a pop! followed by push! writes into the
        # matrix again. A copy-on-modify scheme would need its own bookkeeping.
        pnts = reinterpret(Point{D,T}, vec(stored))
        return new{D,T,typeof(stored),typeof(pnts)}(stored, pnts)
    end
end

"""
    AnyDim

Marker in the dimension slot of `AbsPntSeq{D, T}` for a sequence whose dimension is
known only at run time, such as [`MatVecPntSeq`](@ref). Code compiled for such a
sequence serves every dimension. Use [`dim`](@ref) to get the dimension; functions
that use `D` as a number (for example to build a `Point{D, T}`) do not accept these
sequences.
"""
struct AnyDim end

"""
    ColumnPoint{T, M}

Column `j` of a matrix, as a read-only point (`AbstractVector{T}`) that shares the
matrix's memory: reading it reads the matrix, nothing is copied, and writing through
it (`p[i] = x`) is an error. The points of a [`MatVecPntSeq`](@ref) are of this type.
Arithmetic on them (`p - q`, `2p`) returns a new `Vector{T}`; `copy(p)` and
`collect(p)` give a `Vector{T}` that no longer shares memory.

`ColumnPoint(mat, j)` checks that `mat` has a column `j` and 1-based indices.
"""
struct ColumnPoint{T,M<:AbstractMatrix{T}} <: AbstractVector{T}
    mat::M
    j::Int

    # Unchecked: for ColumnPoints, which checks j itself (or runs under @inbounds).
    ColumnPoint{T,M}(mat::M, j::Int, ::Val{:unchecked}) where {T,M<:AbstractMatrix{T}} =
        new{T,M}(mat, j)
end

function ColumnPoint(mat::AbstractMatrix{T}, j::Integer) where {T}
    Base.require_one_based_indexing(mat)
    checkbounds(mat, :, j)
    return ColumnPoint{T,typeof(mat)}(mat, Int(j), Val(:unchecked))
end

Base.size(p::ColumnPoint) = (size(p.mat, 1),)
Base.IndexStyle(::Type{<:ColumnPoint}) = IndexLinear()
Base.@propagate_inbounds function Base.getindex(p::ColumnPoint, i::Int)
    @boundscheck checkbounds(p, i)
    # In bounds: i was checked above (or by the caller under @inbounds), and j when the
    # point was created.
    return @inbounds p.mat[i, p.j]
end

# The columns of a matrix as a vector of `ColumnPoint`s, created on access.
struct ColumnPoints{T,M<:AbstractMatrix{T}} <: AbstractVector{ColumnPoint{T,M}}
    mat::M
end

Base.size(c::ColumnPoints) = (size(c.mat, 2),)
Base.IndexStyle(::Type{<:ColumnPoints}) = IndexLinear()
Base.@propagate_inbounds function Base.getindex(c::ColumnPoints{T,M}, j::Int) where {T,M}
    @boundscheck checkbounds(c, j)
    return ColumnPoint{T,M}(c.mat, j, Val(:unchecked))
end

"""
    MatVecPntSeq(mat::AbstractMatrix{T}; copy::Bool = false)

A sequence of points stored in a `D × N` matrix, like [`MatPntSeq`](@ref), whose
points are read-only views of the columns ([`ColumnPoint`](@ref)) instead of static
vectors. The dimension is `size(mat, 1)`, read at run time: the type is
`AbsPntSeq{AnyDim, T}` (see [`AnyDim`](@ref)), so the constructor is type-stable and
code compiled for one dimension serves all of them. This is the representation for
high dimensions: loops over the coordinates are not unrolled, and reading a point
copies nothing. In low dimensions `MatPntSeq` is faster (see the guide's "Static or
`Vector` points"). Functions written for any point type accept this sequence (among
them `dist`, `dist_sq` and `convex_comb` on its points, and `geom_length`, `prefix_lengths` and `split_edges` on the sequence); those that need `Point{D, T}` do not.

Memory is shared with the matrix exactly as for `MatPntSeq`: changes to the matrix are
seen by the sequence, and `P[j] = p`, [`translate!`](@ref) and [`scale!`](@ref) write
into the matrix. Points cannot be written through, and the number of points is fixed.
`copy = true` first copies the matrix into a new dense `Matrix{T}`, making the
sequence independent of it.
"""
struct MatVecPntSeq{T,M<:AbstractMatrix{T}} <: AbsPntSeq{AnyDim,T}
    mat::M
    pnts::ColumnPoints{T,M}

    function MatVecPntSeq(mat::AbstractMatrix{T}; copy::Bool = false) where {T}
        Base.require_one_based_indexing(mat)
        # `copy` names the option, so this calls the Matrix constructor, which
        # always allocates a new array. Resizing is not offered, as for MatPntSeq.
        stored = copy ? Matrix{T}(mat) : mat
        return new{T,typeof(stored)}(stored, ColumnPoints(stored))
    end
end

@inline dim(P::MatVecPntSeq) = size(P.mat, 1)

@inline function dist(P::MatVecPntSeq{T}, i::Int, j::Int) where {T}
    i == j && return zero(float(T))
    A = P.mat
    checkbounds(A, :, i)
    checkbounds(A, :, j)
    s = zero(float(T))
    # In bounds: columns i and j were checked above.
    @inbounds for k in axes(A, 1)
        s += abs2(A[k, i] - A[k, j])
    end
    return sqrt(s)
end

function translate_columns!(P::MatVecPntSeq, v::AbstractVector)
    length(v) == dim(P) || throw(
        DimensionMismatch(
            "translation of length $(length(v)) for points of dimension $(dim(P))",
        ),
    )
    P.mat .+= v
    return P
end

translate!(P::MatVecPntSeq, v::AbstractVector) = translate_columns!(P, v)
# Resolves the ambiguity with translate!(::AbsPntSeq{D,T}, ::Point{D,T}) at D = AnyDim.
translate!(P::MatVecPntSeq{T}, v::Point{AnyDim,T}) where {T} = translate_columns!(P, v)

function scale!(P::MatVecPntSeq, s::Real)
    P.mat .*= s
    return P
end

Base.setindex!(P::MatVecPntSeq, v::AbstractVector, i::Int) = (P.mat[:, i] = v; v)
Base.eltype(::Type{<:MatVecPntSeq{T,M}}) where {T,M} = ColumnPoint{T,M}
Base.eltype(::MatVecPntSeq{T,M}) where {T,M} = ColumnPoint{T,M}
Base.Matrix(P::MatVecPntSeq{T}) where {T} = Matrix{T}(P.mat)

function Base.show(io::IO, P::MatVecPntSeq{T}) where {T}
    print(io, "MatVecPntSeq{$T} with $(length(P)) points in dimension $(dim(P))")
end

"""
    dim(P::AbsPntSeq)

Return the dimension of the points in the point sequence.
"""
@inline dim(P::AbsPntSeq{D,T}) where {D,T} = D
# The dimension of an AnyDim sequence of separate points: that of its first point.
function dim(P::GPntSeq{AnyDim})
    isempty(P.pnts) && throw(ArgumentError("dim of an empty sequence with AnyDim"))
    return length(first(P.pnts))
end
dim(::GPntSeq{AnyDim,T,<:AbsVirtualPoint}) where {T} =
    throw(ArgumentError("virtual points have no coordinates, so no dimension"))

"""
    dist(P::AbsPntSeq, i::Int, j::Int)

Distance between vertex `i` and vertex `j` of the point sequence.
Part of the `AbsFMS` interface.
"""
@inline function dist(P::AbsPntSeq{D,T}, i::Int, j::Int) where {D,T}
    if i == j
        return zero(float(T))
    end
    return dist(P.pnts[i], P.pnts[j])
end

# Constructors
"""
    PntSeq{D, T}()

Construct an empty point sequence for dimension `D` and type `T`.
"""
PntSeq{D,T}() where {D,T} = PntSeq{D,T}(Point{D,T}[])

"""
    GPntSeq{D, T, P}()

Construct an empty point sequence with points of type `P`.
"""
GPntSeq{D,T,P}() where {D,T,P<:Union{AbstractVector{T},AbsVirtualPoint}} =
    GPntSeq{D,T,P}(P[])


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
    for i in 1:(length(pnts)-1)
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
    for i in 2:n
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
    OutT = float(T)
    n == 1 && return Point{D,OutT}(pnt_seq.pnts[1])

    plens = prefix_lengths(pnt_seq)
    total_len = plens[end]
    total_len == 0.0 && return Point{D,OutT}(pnt_seq.pnts[1])

    # Target absolute arc-length
    target = t * total_len
    if target <= 0
        return Point{D,OutT}(pnt_seq.pnts[1])
    elseif target >= total_len
        return Point{D,OutT}(pnt_seq.pnts[n])
    end

    # Find the segment containing the target distance
    i = searchsortedfirst(plens, target)
    i <= 1 && return Point{D,OutT}(pnt_seq.pnts[1])

    segment_len = plens[i] - plens[i-1]
    local_t = OutT((target - plens[i-1]) / segment_len)
    return convex_comb(
        Point{D,OutT}(pnt_seq.pnts[i-1]),
        Point{D,OutT}(pnt_seq.pnts[i]),
        local_t,
    )
end

# Base Extensions for collection-like behavior
Base.length(P::AbsPntSeq) = length(P.pnts)
Base.getindex(P::AbsPntSeq, i::Int) = P.pnts[i]
Base.setindex!(P::AbsPntSeq, v, i::Int) = (P.pnts[i] = v)
Base.firstindex(P::AbsPntSeq) = 1
Base.lastindex(P::AbsPntSeq) = length(P.pnts)
Base.iterate(P::AbsPntSeq) = iterate(P.pnts)
Base.iterate(P::AbsPntSeq, state) = iterate(P.pnts, state)
Base.eltype(::Type{<:AbsPntSeq{D,T}}) where {D,T} = Point{D,T}
Base.eltype(::AbsPntSeq{D,T}) where {D,T} = Point{D,T}
Base.eltype(::Type{<:GPntSeq{D,T,P}}) where {D,T,P} = P
Base.eltype(::GPntSeq{D,T,P}) where {D,T,P} = P
Base.push!(P::GPntSeq, p::AbstractVector) = push!(P.pnts, p)
# Points whose length is not in their type: check it against D (or the first point).
function Base.push!(P::GPntSeq{D,T,V}, p::AbstractVector) where {D,T,V<:Vector}
    n = D === AnyDim ? (isempty(P.pnts) ? length(p) : length(first(P.pnts))) : D
    length(p) == n || throw_point_length(length(P.pnts) + 1, length(p), n)
    return push!(P.pnts, p)
end
Base.push!(P::GPntSeq{D,T,V}, p::V) where {D,T,V<:AbsVirtualPoint} = push!(P.pnts, p)
Base.pop!(P::GPntSeq) = pop!(P.pnts)
Base.empty!(P::GPntSeq) = empty!(P.pnts)
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
    for i in 2:(n-1)
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
    m = length(P)
    m == 0 && throw(ArgumentError("Empty point sequence"))
    OutT = float(T)
    out = Vector{Point{D,OutT}}(undef, n)
    if m == 1
        fill!(out, Point{D,OutT}(P[1]))
        return PntSeq(out)
    end
    plens = prefix_lengths(P)
    total_len = plens[end]
    if total_len == 0.0
        fill!(out, Point{D,OutT}(P[1]))
        return PntSeq(out)
    end
    out[1] = Point{D,OutT}(P[1])
    out[n] = Point{D,OutT}(P[m])
    seg = 2
    for i in 2:(n-1)
        target = ((i - 1) / (n - 1)) * total_len
        while seg < m && plens[seg] < target
            seg += 1
        end
        denom = plens[seg] - plens[seg-1]
        local_t = denom == 0.0 ? zero(OutT) : OutT((target - plens[seg-1]) / denom)
        out[i] = convex_comb(Point{D,OutT}(P[seg-1]), Point{D,OutT}(P[seg]), local_t)
    end
    return PntSeq(out)
end

"""
    push_distinct!(P, p)

Append `p` to the `PntSeq` `P` unless it equals the last vertex. Returns `true` if
`p` was appended.
"""
function push_distinct!(P::PntSeq{D,T}, p::Point{D,T}) where {D,T}
    if isempty(P.pnts) || dist_sq(last(P.pnts), p) > 0
        push!(P.pnts, p)
        return true
    end
    return false
end

# Any point type; push! checks the length of points of run-time length.
function push_distinct!(P::GPntSeq, p::AbstractVector)
    if isempty(P.pnts) || dist_sq(last(P.pnts), p) > 0
        push!(P, p)
        return true
    end
    return false
end

"""
    split_edges(P)

Return the curve through the vertices of `P` with the midpoint of every edge inserted.
Consecutive repeated points are dropped. Requires floating-point coordinates.
"""
function split_edges(P::AbsPntSeq{D,T}) where {D,T<:AbstractFloat}
    n = length(P)
    Q = PntSeq{D,T}()
    n == 0 && return Q
    for i in 1:(n-1)
        push_distinct!(Q, P[i])
        push_distinct!(Q, convex_comb(P[i], P[i+1], 0.5))
    end
    push_distinct!(Q, P[n])
    return Q
end

split_edges(::AbsPntSeq) =
    throw(ArgumentError("split_edges requires floating-point coordinates"))

# Points of run-time length: the result has `Vector` points (for a `MatVecPntSeq`, a
# new `GPntSeq{AnyDim}`; the matrix is not resized).
split_edges(P::GPntSeq{D,T,V}) where {D,T<:AbstractFloat,V<:Vector{T}} =
    split_edges_into!(GPntSeq{D,T,V}(), P)
split_edges(P::MatVecPntSeq{T}) where {T<:AbstractFloat} =
    split_edges_into!(GPntSeq{AnyDim,T,Vector{T}}(), P)

function split_edges_into!(Q::GPntSeq, P::AbsPntSeq)
    n = length(P)
    n == 0 && return Q
    for i in 1:(n-1)
        push_distinct!(Q, P[i])
        push_distinct!(Q, convex_comb(P[i], P[i+1], 0.5))
    end
    push_distinct!(Q, P[n])
    return Q
end

"""
    reverse(P::PntSeq)

Return the point sequence with its vertices in reverse order.
"""
Base.reverse(P::PntSeq{D,T}) where {D,T} = PntSeq{D,T}(reverse(P.pnts))

"""
    rand_pnt_seq(D, T, n)

Generate a random point sequence in D dimensions with `n` vertices, each coordinate in `[0, 1)`.
"""
function rand_pnt_seq(D::Int, T::Type{<:Real}, n::Int)
    return PntSeq([rand_point(D, T) for _ in 1:n])
end

"""
    pnt_seq_random_sphere(D, T, n)

Generate a random point sequence in D dimensions with `n` vertices sampled uniformly from the unit sphere.
"""
function pnt_seq_random_sphere(D::Int, T::Type{<:Real}, n::Int)
    pnts = Vector{Point{D,T}}(undef, n)
    for i in 1:n
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

"""
    Polygon2I

Alias for `PntSeq2I` (the pre-0.5 name of point sequences).
"""
const Polygon2I = PntSeq2I
"""
    Polygon2F

Alias for `PntSeq2F` (the pre-0.5 name of point sequences).
"""
const Polygon2F = PntSeq2F
"""
    Polygon3F

Alias for `PntSeq3F` (the pre-0.5 name of point sequences).
"""
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

export PntSeq,
    GPntSeq,
    AbsPntSeq,
    AbsVirtualPoint,
    MatPntSeq,
    MatVecPntSeq,
    AnyDim,
    PntSeq2I,
    PntSeq2F,
    PntSeq3F
export Polygon, Polygon2I, Polygon2F, Polygon3F
export dim, geom_length
export translate!, scale!, simplify, sample_uniformly
export rand_pnt_seq
