###############################################
### Hilbert spaces known through squared distances, and virtual points

# A space's input vertices 1:n are known only through the squared distances an oracle
# computes for two of them. Every such value is computed at most once and cached in
# square blocks over the unordered pairs, allocated when first touched. A `VPoint` is
# a point (1 - t) x_a + t x_b on the segment of two input vertices; geometry on such
# points uses only the cached squared distances.

const HS_BLOCK = 256                      # a cache block covers HS_BLOCK^2 pairs
const HS_BLOCK_BYTES = HS_BLOCK^2 * sizeof(Float64)
const HS_EMPTY = Matrix{Float64}(undef, 0, 0)  # the slot of a block not allocated

# A slot of the block directory. The block is published with a release store after it
# has been filled with NaN, and read with an acquire load, so no thread sees an
# uninitialized block.
mutable struct HSSlot
    @atomic blk::Matrix{Float64}
end

# The cache: blocks of the upper triangle of the n x n matrix of squared distances.
# Counters are updated only on misses, never on hits.
mutable struct HSCache
    const n::Int
    const nb::Int                         # blocks per side
    const slots::Vector{HSSlot}           # block (bi, bj), bi <= bj, at slot_index
    const lock::ReentrantLock             # guards block allocation and empty!
    const max_blocks::Int
    @atomic n_blocks::Int
    @atomic oracle_calls::Int
    @atomic uncached::Int                 # misses not stored (memory budget reached)
end

function HSCache(n::Int, max_bytes::Integer)
    nb = cld(n, HS_BLOCK)
    slots = [HSSlot(HS_EMPTY) for _ in 1:(nb*(nb+1)÷2)]
    max_blocks = Int(min(max_bytes ÷ HS_BLOCK_BYTES, typemax(Int)))
    return HSCache(n, nb, slots, ReentrantLock(), max_blocks, 0, 0, 0)
end

@inline hs_slot_index(bi::Int, bj::Int) = bj * (bj - 1) ÷ 2 + bi

"""
    HilbertSpace(items, dist2; max_bytes = typemax(Int))
    HilbertSpace(items; kernel = k, max_bytes = typemax(Int))
    HilbertSpace(points::AbstractVector{<:AbstractVector{<:Real}}; max_bytes = typemax(Int))
    HilbertSpace(pts::AbstractVector{<:VPoint}; max_bytes = typemax(Int))

A Hilbert space known only through the squared distances of its input vertices
`1:length(items)`: the points of a dataset, possibly in an infinite-dimensional space,
such as the feature space of a kernel. Curves through its vertices ([`vcurve`](@ref))
are sequences of virtual points ([`VPoint`](@ref)); geometry on them uses only squared
distances, each computed at most once and cached.

The squared distance of two items `x`, `y` is

- `dist2(x, y)`, given directly (the most accurate choice when it can be computed
  without cancellation: for a Gaussian kernel, `-2expm1(-|x - y|^2 / s^2)`);
- `k(x, x) + k(y, y) - 2k(x, y)` for a positive-definite `kernel` `k` (clamped at 0;
  the diagonal is computed once);
- the Euclidean squared distance, for items that are coordinate vectors;
- for virtual points of another space, their squared distance there: the points become
  the input vertices of the new space (re-basing).

The distance function must be pure and thread-safe: threads share the space and may
compute an entry concurrently. It must return finite non-negative values, and must come
from a Hilbert space (a positive-definite kernel): otherwise the geometry is meaningless.

The cache holds the computed squared distances in blocks of 256 x 256 pairs (512 KB),
allocated when first touched. With `max_bytes`, no blocks are allocated beyond that
size; further values are computed without being stored (`cache_stats` counts them).
`empty!(S)` frees the cache. `deepcopy` returns the space itself, so morphings copied
with `deepcopy` keep their points in the same space.

A `HilbertSpace` is a finite metric space (`AbsFMS`): `size(S)` and `dist(S, i, j)`.
"""
mutable struct HilbertSpace{O} <: AbsFMS
    const n::Int
    const oracle::O
    const cache::HSCache
end

function hilbert_space(oracle, n::Integer, max_bytes::Integer)
    0 <= n <= typemax(Int32) || throw(ArgumentError("a space has 0 to 2^31 - 1 vertices"))
    max_bytes >= 0 || throw(ArgumentError("max_bytes must be non-negative"))
    return HilbertSpace(Int(n), oracle, HSCache(Int(n), max_bytes))
end

# Oracles: callable on two vertex indices, returning a Float64.
struct HSDist2{V,F}
    items::V
    f::F
end
(o::HSDist2)(i::Int, j::Int) = Float64(o.f(o.items[i], o.items[j]))

struct HSKernel{V,F}
    items::V
    k::F
    diag::Vector{Float64}
end
(o::HSKernel)(i::Int, j::Int) =
    max(0.0, o.diag[i] + o.diag[j] - 2 * Float64(o.k(o.items[i], o.items[j])))

struct HSEuclid{V}
    items::V
end
function (o::HSEuclid)(i::Int, j::Int)
    x, y = o.items[i], o.items[j]
    length(x) == length(y) || throw(DimensionMismatch("points $i and $j"))
    s = 0.0
    for k in eachindex(x, y)
        s += Float64(x[k] - y[k])^2
    end
    return s
end

struct HSRebased{P}
    pts::Vector{P}
end
(o::HSRebased)(i::Int, j::Int) = dist_sq(o.pts[i], o.pts[j])

HilbertSpace(items::AbstractVector, dist2; max_bytes::Integer = typemax(Int)) =
    hilbert_space(HSDist2(items, dist2), length(items), max_bytes)

function HilbertSpace(
    items::AbstractVector;
    kernel = nothing,
    max_bytes::Integer = typemax(Int),
)
    if kernel === nothing
        eltype(items) <: AbstractVector{<:Real} || throw(
            ArgumentError(
                "give a kernel or a dist2 function for items that are not vectors",
            ),
        )
        return hilbert_space(HSEuclid(items), length(items), max_bytes)
    end
    diag = [Float64(kernel(x, x)) for x in items]
    return hilbert_space(HSKernel(items, kernel, diag), length(items), max_bytes)
end

"""
    VPoint(S::HilbertSpace, a, b, t)

The virtual point `(1 - t) x_a + t x_b` of the space `S`, on the segment of its input
vertices `a` and `b` (any two, not only consecutive ones), with `t` in [0, 1]; 24 bytes,
holding a reference to `S`. [`vpoint`](@ref)`(S, i)` is the input vertex `i`.

The representation is canonical: `a < b`, and an input vertex has `a == b` and `t == 0`.
So `==` and `hash` agree with geometry for points of the same pair of vertices (two
input vertices at distance zero are still different points).

`dist`, `dist_sq` and `convex_comb` take virtual points of one space; so do
[`local_frame`](@ref BasicCompGeometry.local_frame), `BasicCompGeometry.line_params`, `BasicCompGeometry.inner_diff`
and `BasicCompGeometry.vscale`. A convex combination of two points is representable
only if they lie on one segment of two input vertices; otherwise it throws an
`ArgumentError` (re-base the points into a new space with
[`HilbertSpace`](@ref)`(pts)` first).
"""
struct VPoint{S<:HilbertSpace} <: AbsVirtualPoint
    space::S
    a::Int32
    b::Int32
    t::Float64
    function VPoint(S::HilbertSpace, a::Integer, b::Integer, t::Real)
        (1 <= a <= S.n && 1 <= b <= S.n) || throw(BoundsError(S, (a, b)))
        0 <= t <= 1 || throw(ArgumentError("VPoint parameter $t is not in [0, 1]"))
        if a > b
            a, b, t = b, a, 1 - t
        end
        t == 1 && return new{typeof(S)}(S, b, b, 0.0)
        (t == 0 || a == b) && return new{typeof(S)}(S, a, a, 0.0)
        return new{typeof(S)}(S, a, b, Float64(t))
    end
end

"""
    vpoint(S::HilbertSpace, i)

The input vertex `i` of `S`, as a [`VPoint`](@ref).
"""
vpoint(S::HilbertSpace, i::Integer) = VPoint(S, i, i, 0.0)

"""
    vcurve(S::HilbertSpace, idx)

The curve through the input vertices `idx` of `S`, in order: a
`GPntSeq{AnyDim, Float64, VPoint}`.
"""
vcurve(S::HilbertSpace, idx) =
    GPntSeq{AnyDim,Float64,VPoint{typeof(S)}}([vpoint(S, i) for i in idx])

HilbertSpace(pts::AbstractVector{<:VPoint}; max_bytes::Integer = typemax(Int)) =
    hilbert_space(HSRebased(collect(pts)), length(pts), max_bytes)

Base.deepcopy_internal(S::HilbertSpace, ::IdDict) = S
Base.show(io::IO, S::HilbertSpace) = print(io, "HilbertSpace with ", S.n, " vertices")
Base.show(io::IO, x::VPoint) =
    x.a == x.b ? print(io, "VPoint(", x.a, ")") :
    print(io, "VPoint(", x.a, ", ", x.b, ", ", x.t, ")")

Base.:(==)(x::VPoint, y::VPoint) =
    x.space === y.space && x.a == y.a && x.b == y.b && x.t == y.t
Base.hash(x::VPoint, h::UInt) = hash((objectid(x.space), x.a, x.b, x.t), h)
# The constructor rejects NaN and infinite parameters.
Base.isnan(::VPoint) = false
Base.isfinite(::VPoint) = true

is_input(x::VPoint) = x.a == x.b

"""
    same_space(x::VPoint, y::VPoint)

Throws an `ArgumentError` unless `x` and `y` are points of the same space.
"""
@inline function same_space(x::VPoint, y::VPoint)
    x.space === y.space && return nothing
    throw(ArgumentError("virtual points of different spaces"))
end

"""
    is_vertex_curve(P)

Whether every vertex of the virtual curve `P` is an input vertex of its space.
"""
is_vertex_curve(P::GPntSeq{AnyDim,Float64,<:VPoint}) = all(is_input, P.pnts)

#########################################################################
# The cache

"""
    BasicCompGeometry.vdist2(S::HilbertSpace, i, j)

The squared distance of the input vertices `i` and `j`, from the cache, computed by the
oracle on the first request.
"""
@inline function vdist2(S::HilbertSpace, i::Integer, j::Integer)
    i == j && return 0.0
    i, j = minmax(Int(i), Int(j))
    c = S.cache
    bi, bj = (i - 1) ÷ HS_BLOCK + 1, (j - 1) ÷ HS_BLOCK + 1
    slot = @inbounds c.slots[hs_slot_index(bi, bj)]
    blk = @atomic :acquire slot.blk
    isempty(blk) && return vdist2_new_block(S, slot, i, j)
    k = ((j - 1) % HS_BLOCK) * HS_BLOCK + (i - 1) % HS_BLOCK + 1
    v = GC.@preserve blk unsafe_load(pointer(blk), k, :monotonic)
    isnan(v) || return v
    v = oracle_dist2(S, i, j)
    GC.@preserve blk unsafe_store!(pointer(blk), v, k, :monotonic)
    return v
end

function oracle_dist2(S::HilbertSpace, i::Int, j::Int)
    @atomic S.cache.oracle_calls += 1
    v = S.oracle(i, j)::Float64
    (isfinite(v) && v >= 0) ||
        throw(ArgumentError("squared distance of vertices $i and $j is $v"))
    return v
end

# The first request in a block: allocate it (unless the budget is reached), then look up.
@noinline function vdist2_new_block(S::HilbertSpace, slot::HSSlot, i::Int, j::Int)
    c = S.cache
    lock(c.lock) do
        isempty(@atomic :acquire slot.blk) || return nothing
        (@atomic c.n_blocks) < c.max_blocks || return nothing
        @atomic :release slot.blk = fill(NaN, HS_BLOCK, HS_BLOCK)
        @atomic c.n_blocks += 1
        return nothing
    end
    isempty(@atomic :acquire slot.blk) || return vdist2(S, i, j)
    @atomic c.uncached += 1
    return oracle_dist2(S, i, j)
end

"""
    empty!(S::HilbertSpace)

Frees the cached squared distances of `S`; they are computed again when needed.
"""
function Base.empty!(S::HilbertSpace)
    c = S.cache
    lock(c.lock) do
        foreach(s -> @atomic(:release, s.blk = HS_EMPTY), c.slots)
        @atomic c.n_blocks = 0
    end
    return S
end

"""
    cache_stats(S::HilbertSpace)

`(oracle_calls, uncached, blocks, bytes)`: the calls of the distance function so far,
those whose value was not stored (memory budget), and the allocated cache blocks and
their size.
"""
function cache_stats(S::HilbertSpace)
    c = S.cache
    nb = @atomic c.n_blocks
    return (
        oracle_calls = (@atomic c.oracle_calls),
        uncached = (@atomic c.uncached),
        blocks = nb,
        bytes = nb * HS_BLOCK_BYTES,
    )
end

Base.size(S::HilbertSpace) = S.n
dist(S::HilbertSpace, i::Integer, j::Integer) = sqrt(vdist2(S, i, j))

#########################################################################
# Geometry on virtual points

"""
    dist_sq(x::VPoint, y::VPoint)

The squared distance of two virtual points of one space, from the squared distances of
the (at most four) input vertices involved: `|Σ λ_i x_i|^2 = -Σ_{i<j} λ_i λ_j D_ij^2` for
weights summing to zero (the identity needs no distinct vertices), clamped at 0. Points of one pair of vertices give
`(t - s)^2 D_ab^2`. Accuracy: for points at distance δ, with vertices at distance up to
D, the absolute error is about ε D^2 / δ (inherent to computing from squared distances).
"""
function dist_sq(x::VPoint, y::VPoint)
    same_space(x, y)
    S = x.space
    if x.a == y.a && x.b == y.b
        return x.a == x.b ? 0.0 : (x.t - y.t)^2 * vdist2(S, x.a, x.b)
    end
    # x - y = Σ w_k v_k with Σ w_k = 0, so |x - y|^2 = -Σ_{k<l} w_k w_l D²(v_k, v_l).
    # Repeated vertices need no merging: their D² is 0.
    v = (x.a, x.b, y.a, y.b)
    w = (1 - x.t, x.t, -(1 - y.t), -y.t)
    s = 0.0
    for k in 1:4, l in (k+1):4
        s -= w[k] * w[l] * vdist2(S, v[k], v[l])
    end
    return max(s, 0.0)
end

dist(x::VPoint, y::VPoint) = sqrt(dist_sq(x, y))

# The position of x on the segment of the input vertices (u, v), u < v, or NaN.
@inline function param_on(x::VPoint, u::Int32, v::Int32)
    x.a == u && x.b == v && return x.t
    is_input(x) || return NaN
    return x.a == u ? 0.0 : x.a == v ? 1.0 : NaN
end

"""
    convex_comb(x::VPoint, y::VPoint, s)

The point `(1 - s) x + s y`. Both points must lie on one segment of two input vertices
(two input vertices always do); otherwise an `ArgumentError` is thrown.
"""
function convex_comb(x::VPoint, y::VPoint, s::Real)
    same_space(x, y)
    if is_input(x) && is_input(y)
        return VPoint(x.space, x.a, y.a, s)
    end
    u, v = is_input(x) ? (y.a, y.b) : (x.a, x.b)
    tx, ty = param_on(x, u, v), param_on(y, u, v)
    (isnan(tx) || isnan(ty)) && throw(
        ArgumentError(
            "$x and $y do not lie on one segment of two input vertices; " *
            "re-base the points with HilbertSpace(points)",
        ),
    )
    return VPoint(x.space, u, v, (1 - s) * tx + s * ty)
end

"""
    BasicCompGeometry.inner_diff(x, y, z, w)

The inner product `⟨x - y, z - w⟩` of virtual points of one space, by polarization:
`(|x - w|^2 + |y - z|^2 - |x - z|^2 - |y - w|^2) / 2`.
"""
inner_diff(x::VPoint, y::VPoint, z::VPoint, w::VPoint) =
    (dist_sq(x, w) + dist_sq(y, z) - dist_sq(x, z) - dist_sq(y, w)) / 2

"""
    BasicCompGeometry.VDiff

`x - y` for virtual points: a formal difference, accepted by `LinearAlgebra.dot`
(through [`BasicCompGeometry.inner_diff`](@ref)), so that code written for coordinate
points that forms differences only to take inner products works on virtual points.
"""
struct VDiff{S}
    x::VPoint{S}
    y::VPoint{S}
end

Base.:-(x::VPoint, y::VPoint) = (same_space(x, y); VDiff(x, y))
LinearAlgebra.dot(d::VDiff, e::VDiff) = inner_diff(d.x, d.y, e.x, e.y)

# Twice the area of a triangle with sides a, b, c (Kahan's formula, accurate for
# needle-like triangles).
function kahan_2area(a::Float64, b::Float64, c::Float64)
    a, b = max(a, b), min(a, b)
    b, c = max(b, c), min(b, c)
    a, b = max(a, b), min(a, b)
    x = (a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c))
    return x <= 0 ? 0.0 : sqrt(x) / 2
end

"""
    BasicCompGeometry.line_params(q, a, b) -> (t0, h2)

For virtual points `q`, `a != b` of one space: the parameter `t0` of the point of the
line through `a` and `b` closest to `q` (`a + t0 (b - a)`), and the squared distance
`h2` from `q` to that line. For three input vertices, from the law of cosines and
Kahan's triangle area; otherwise through [`local_frame`](@ref BasicCompGeometry.local_frame).
"""
function line_params(q::VPoint, a::VPoint, b::VPoint)
    same_space(q, a)
    same_space(q, b)
    if is_input(q) && is_input(a) && is_input(b)
        S = q.space
        L2 = vdist2(S, a.a, b.a)
        L2 > 0 || throw(ArgumentError("line_params of a degenerate segment"))
        da, db = vdist2(S, q.a, a.a), vdist2(S, q.a, b.a)
        h = kahan_2area(sqrt(da), sqrt(db), sqrt(L2)) / sqrt(L2)
        return ((da + L2 - db) / (2 * L2), h * h)
    end
    Q, A, B = local_frame(q, a, b)
    d = B - A
    L2 = dot(d, d)
    L2 > 0 || throw(ArgumentError("line_params of a degenerate segment"))
    t0 = dot(Q - A, d) / L2
    return (t0, max(0.0, dist_sq(Q, A + t0 * d)))
end

"""
    BasicCompGeometry.vscale(points::VPoint...)

The largest distance between the input vertices involved in `points`: the scale for
thresholds on values computed from them (errors are relative to it).
"""
function vscale(points::Vararg{VPoint,N}) where {N}
    S = first(points).space
    vs = sorted_vertices(points)
    m = count(!=(0), vs)
    s = 0.0
    for i in 1:m, j in (i+1):m
        s = max(s, vdist2(S, vs[i], vs[j]))
    end
    return sqrt(s)
end

# The distinct input vertices of the points, sorted, padded with 0 (at most 8).
@inline function sorted_vertices(points::NTuple{N,VPoint}) where {N}
    N <= 4 || throw(ArgumentError("at most 4 virtual points"))
    raw = ntuple(
        i -> i <= 2N ? (isodd(i) ? points[(i+1)÷2].a : points[i÷2].b) : typemax(Int32),
        8,
    )
    srt = sort(SVector(raw))
    # The k-th distinct value of srt (0 past the last one).
    function kth(k)
        c = 0
        for i in 1:8
            srt[i] == typemax(Int32) && return Int32(0)
            (i == 1 || srt[i] != srt[i-1]) && (c += 1) == k && return srt[i]
        end
        return Int32(0)
    end
    return SVector(ntuple(kth, 8))
end

"""
    local_frame(points::VPoint...) -> NTuple{N, Point{3, Float64}}

The points of one space placed isometrically in ℝ³: the distinct input vertices
involved (at most four) are sorted by index, the first is placed at the origin, and the
others from the Gram matrix of their differences (polarization of the cached squared
distances) by a pivoted Cholesky factorization, with negative pivots (rounding)
clamped to 0. Each point is then `(1 - t) v_a + t v_b`. Euclidean code can run on the
result: distances, closest points and inner products in the frame are those of the
space, up to rounding (relative to [`BasicCompGeometry.vscale`](@ref)).
"""
function local_frame(points::Vararg{VPoint,N}) where {N}
    S = first(points).space
    foreach(x -> same_space(x, first(points)), points)
    vs = sorted_vertices(points)
    m = count(!=(0), vs)
    m <= 4 || throw(ArgumentError("local_frame takes points of at most 4 input vertices"))
    coords = frame_coords(S, vs, m)
    pos(v) = coords[vertex_slot(vs, v)]
    return map(x -> (1 - x.t) * pos(x.a) + x.t * pos(x.b), points)
end

# The position of the vertex v in vs (present by construction).
@inline function vertex_slot(vs::SVector{8,Int32}, v::Int32)
    for i in 1:8
        vs[i] == v && return i
    end
    return 1
end

# Coordinates of the m <= 4 vertices vs[1:m] in ℝ³: vs[1] at the origin, the rest from
# the pivoted Cholesky factorization of G[i, j] = <v_i - v_1, v_j - v_1>.
function frame_coords(S::HilbertSpace, vs::SVector{8,Int32}, m::Int)
    k = m - 1
    function g(l)
        i, j = (l - 1) % 3 + 1, (l - 1) ÷ 3 + 1
        (i > k || j > k) && return 0.0
        d1i, d1j = vdist2(S, vs[1], vs[i+1]), vdist2(S, vs[1], vs[j+1])
        return i == j ? d1i : (d1i + d1j - vdist2(S, vs[i+1], vs[j+1])) / 2
    end
    L = pivoted_cholesky3(SMatrix{3,3,Float64,9}(ntuple(g, 9)), k)
    row(i) = Point{3,Float64}(L[i, 1], L[i, 2], L[i, 3])
    return ntuple(i -> (i == 1 || i > m) ? zero(Point{3,Float64}) : row(i - 1), 8)
end

# L with G = L L^T for the leading k x k block of a positive-semidefinite G (k <= 3),
# pivoting on the largest remaining diagonal; a non-positive pivot ends the
# factorization (the rest of the block is rank-deficient up to rounding).
function pivoted_cholesky3(G::SMatrix{3,3,Float64,9}, k::Int)
    A, L = G, zero(SMatrix{3,3,Float64,9})
    done = (false, false, false)
    for col in 1:k
        p, best = 0, 0.0
        for i in 1:k
            (done[i] || A[i, i] <= best) || ((p, best) = (i, A[i, i]))
        end
        p == 0 && break
        done = Base.setindex(done, true, p)
        r = sqrt(best)
        for i in 1:k
            (done[i] && i != p) ||
                (L = Base.setindex(L, i == p ? r : A[i, p] / r, i + 3(col - 1)))
        end
        for i in 1:k, j in 1:k
            (done[i] || done[j]) ||
                (A = Base.setindex(A, A[i, j] - L[i, col] * L[j, col], i + 3(j - 1)))
        end
    end
    return L
end

#########################################################################
# Curves of virtual points

push_distinct!(P::GPntSeq{AnyDim,Float64,V}, p::V) where {V<:VPoint} =
    (isempty(P.pnts) || dist_sq(last(P.pnts), p) > 0) ? (push!(P.pnts, p); true) : false

split_edges(P::GPntSeq{AnyDim,Float64,V}) where {V<:VPoint} =
    split_edges_into!(GPntSeq{AnyDim,Float64,V}(), P)

export HilbertSpace, VPoint, vpoint, vcurve
