"""
    PntSeqHierarchy

Hierarchical simplification of a point sequence (polygonal curve).

A *palette* records, for a curve `P[1..n]`, the error of a recursive subdivision
of `1:n`: a range `i:j` is split at an interior vertex `k` chosen by a
[`SplitRule`](@ref), and `err[k]` is the largest error of the two sub-ranges
(measured by a user supplied function `err(P, i:j)`, the error of replacing
`P[i:j]` by the segment `P[i]P[j]`). All vertices above an error threshold can then be
extracted in time proportional to their number.

A [`SimplificationHierarchy`](@ref) stores simplifications of the curve at
decreasing widths, built from a [`HierarchyMeasure`](@ref): the bundle of
functions that define the error being controlled. The default measure uses the
Hausdorff distance ([`hausdorff_measure`](@ref)); other measures (for example a
Frechet distance) only have to provide their own functions.

Names are not exported; use qualified calls such as
`BasicCompGeometry.PntSeqHierarchy.palette`.
"""
module PntSeqHierarchy

import ..BasicCompGeometry:
    AbsPntSeq,
    PntSeq,
    cardin,
    dist_point_segment,
    geom_length,
    hausdorff_dist_subseg,
    hausdorff_simplify

export SplitRule,
    MidpointSplit,
    FarthestSplit,
    Palette,
    palette,
    palette_extract,
    palette_level,
    HierarchyMeasure,
    hausdorff_measure,
    SimplificationHierarchy,
    hierarchy_init,
    hierarchy_push!,
    hierarchy_find,
    hierarchy_approx,
    hierarchy_build!,
    hierarchy_print

#######################################################################
# Split rules
#######################################################################

"""
    SplitRule

Supertype of the rules choosing the vertex at which a range `i:j` of a curve is
subdivided. `split_vertex(rule, P, i, j)` returns a vertex `k`; a range with no
valid interior vertex (`!(i < k < j)`) is a leaf.
"""
abstract type SplitRule end

"""
    MidpointSplit()

Split at the middle index `(i + j) >> 1`. The resulting tree is implicit, so a
palette needs only its error vector.
"""
struct MidpointSplit <: SplitRule end

"""
    FarthestSplit()

Split at the interior vertex farthest from the segment `P[i]P[j]`
(Douglas-Peucker). Ties go to the first such vertex.
"""
struct FarthestSplit <: SplitRule end

@inline split_vertex(::MidpointSplit, P, i::Int, j::Int) = (i + j) >> 1

function split_vertex(::FarthestSplit, P, i::Int, j::Int)
    best = 0
    best_d = -Inf
    for k = (i+1):(j-1)
        d = dist_point_segment(P[k], P[i], P[j])
        if d > best_d
            best_d = d
            best = k
        end
    end
    return best
end

#######################################################################
# Palette
#######################################################################

"""
    Palette{T, R}

Subdivision errors of a curve, indexed by vertex: `pal[k]` is the error recorded
where vertex `k` splits a range (zero for vertices that never split one).
It behaves as an `AbstractVector{T}`.

For [`MidpointSplit`](@ref) the subdivision tree is implicit and a plain vector of
errors can be wrapped in `Palette(errs)` at no cost. Other rules store the
children of every split vertex.
"""
struct Palette{T,R<:SplitRule} <: AbstractVector{T}
    err::Vector{T}
    kids::Vector{NTuple{2,Int}}   # (left, right) child split of vertex k; empty if implicit
    root::Int                     # split vertex of the whole curve; unused if implicit
    rule::R
end

Palette(err::Vector{T}) where {T} =
    Palette{T,MidpointSplit}(err, NTuple{2,Int}[], 0, MidpointSplit())

Base.size(pal::Palette) = size(pal.err)
Base.getindex(pal::Palette, i::Int) = pal.err[i]
Base.IndexStyle(::Type{<:Palette}) = IndexLinear()

@inline function _mid_or_leaf(i::Int, j::Int)
    m = (i + j) >> 1
    return (m == i || m == j) ? 0 : m
end

@inline _root(pal::Palette{T,MidpointSplit}, n::Int) where {T} = _mid_or_leaf(1, n)
@inline _root(pal::Palette, n::Int) = pal.root

@inline _kids(::Palette{T,MidpointSplit}, i::Int, k::Int, j::Int) where {T} =
    (_mid_or_leaf(i, k), _mid_or_leaf(k, j))
@inline _kids(pal::Palette, i::Int, k::Int, j::Int) = pal.kids[k]

function _fill!(errs, kids, P, err, rule, i::Int, j::Int)
    i > j && return (0.0, 0)
    k = split_vertex(rule, P, i, j)
    (k <= i || k >= j) && return (0.0, 0)

    v_l, k_l = _fill!(errs, kids, P, err, rule, i, k)
    v_r, k_r = _fill!(errs, kids, P, err, rule, k, j)
    v = err(P, i:j)

    errs[k] = max(v_l, v_r)
    isempty(kids) || (kids[k] = (k_l, k_r))
    return (v, k)
end

"""
    palette(P, err; rule = MidpointSplit())

Build the [`Palette`](@ref) of the curve `P`. `err(P, i:j)` is the error of
replacing `P[i:j]` by the segment `P[i]P[j]`, for example
`hausdorff_dist_subseg`. The coordinates must be floating point. With
[`FarthestSplit`](@ref) the construction recurses once per split and
costs O(n) per split, so a badly skewed subdivision (for example a spiral
with a very large number of vertices) takes quadratic time and deep recursion;
[`MidpointSplit`](@ref) has logarithmic depth.
"""
function palette(
    P::AbsPntSeq{D,T},
    err;
    rule::SplitRule = MidpointSplit(),
) where {D,T<:AbstractFloat}
    n = length(P)
    errs = zeros(T, n)
    kids = rule isa MidpointSplit ? NTuple{2,Int}[] : fill((0, 0), n)
    n < 2 && return Palette(errs, kids, 0, rule)

    _, root = _fill!(errs, kids, P, err, rule, 1, n)
    return Palette(errs, kids, root, rule)
end

function _p_extract!(pout, pout_is, P, pal, i::Int, j::Int, k::Int, w)
    k == 0 && return 0.0
    if pal[k] < w
        push!(pout, P[k])
        push!(pout_is, k)
        return pal[k]
    end
    k_l, k_r = _kids(pal, i, k, j)
    w_l = _p_extract!(pout, pout_is, P, pal, i, k, k_l, w)
    push!(pout, P[k])
    push!(pout_is, k)
    w_r = _p_extract!(pout, pout_is, P, pal, k, j, k_r, w)
    return max(w_l, w_r)
end

"""
    palette_extract(P, pal, w)

Extract an approximation of `P` from its palette (running time proportional to
the output size). Returns `(Q, indices, w_out)`: `Q` is the simplified curve,
`indices` are the vertices of `P` it keeps, and `w_out` bounds the error of
the subranges that were not subdivided further, which is below `w` whenever
`w` exceeds the largest error in the palette. `pal` may be a `Palette` or a
plain error vector, which is taken to use [`MidpointSplit`](@ref). A curve of at
most four vertices is returned whole, with `w_out = 0`, whatever `w` is.
"""
function palette_extract(P::AbsPntSeq{D,T}, pal::Palette, w::Real) where {D,T}
    len = length(pal)
    @assert len == cardin(P)

    pout = PntSeq{D,T}()
    pout_is = Vector{Int}()
    if cardin(P) <= 4
        for i = 1:len
            push!(pout, P[i])
            push!(pout_is, i)
        end
        return pout, pout_is, zero(T)
    end

    push!(pout, P[1])
    push!(pout_is, 1)

    w = _p_extract!(pout, pout_is, P, pal, 1, len, _root(pal, len), w)

    push!(pout, P[len])
    push!(pout_is, len)

    return pout, pout_is, w
end

palette_extract(P::AbsPntSeq{D,T}, plt::Vector{T}, w::Real) where {D,T} =
    palette_extract(P, Palette(plt), w)

function _p_extract_level!(pout, pout_is, P, pal, i::Int, j::Int, k::Int, level::Int)
    (k == 0 || level <= 0) && return nothing
    k_l, k_r = _kids(pal, i, k, j)
    _p_extract_level!(pout, pout_is, P, pal, i, k, k_l, level - 1)
    push!(pout, P[k])
    push!(pout_is, k)
    _p_extract_level!(pout, pout_is, P, pal, k, j, k_r, level - 1)
    return nothing
end

"""
    palette_level(P, pal, level)

Extract the approximation of `P` made of the first `level` levels of the
subdivision tree (the first and last vertices are always kept). Returns
`(Q, indices)`. For a one-vertex curve the vertex appears twice.
"""
function palette_level(P::AbsPntSeq{D,T}, pal::Palette, level::Int) where {D,T}
    len = cardin(P)
    @assert len > 0

    pout = PntSeq{D,T}()
    pout_is = Vector{Int}()

    push!(pout, P[1])
    push!(pout_is, 1)

    _p_extract_level!(pout, pout_is, P, pal, 1, len, _root(pal, len), level)

    push!(pout, P[len])
    push!(pout_is, len)
    return pout, pout_is
end

palette_level(P::AbsPntSeq{D,T}, plt::Vector{T}, level::Int) where {D,T} =
    palette_level(P, Palette(plt), level)

#######################################################################
# Measure
#######################################################################

"""
    HierarchyMeasure(err, simplify, certify, rule)

The functions defining the error controlled by a hierarchy:

- `err(P, i:j)`: error of replacing `P[i:j]` by the segment `P[i]P[j]`.
- `simplify(P, w)`: a simplification of `P` with error at most about `w`, as
  `(Q, indices)` where `Q[t] == P[indices[t]]`.
- `certify(P, Q, indices)`: a computed bound on the error between `P` and its
  simplification `Q`, with `Q[t] == P[indices[t]]`.
- `rule`: the [`SplitRule`](@ref) used for palettes.
"""
struct HierarchyMeasure{E,S,C,R<:SplitRule}
    err::E
    simplify::S
    certify::C
    rule::R
end

function _hausdorff_certify(P, R, indices)
    m = 0.0
    for t = 1:(length(indices)-1)
        m = max(m, hausdorff_dist_subseg(P, indices[t]:indices[t+1]))
    end
    return m
end

"""
    hausdorff_measure()

The default [`HierarchyMeasure`](@ref): subcurve error is the Hausdorff
distance to its chord, simplification is `hausdorff_simplify`, and palettes
split at the vertex farthest from the chord.
"""
hausdorff_measure() = HierarchyMeasure(
    hausdorff_dist_subseg,
    hausdorff_simplify,
    _hausdorff_certify,
    FarthestSplit(),
)

#######################################################################
# Hierarchy
#######################################################################

"""
    SimplificationHierarchy{D, T, M}

Simplifications of the curve `P` stored from the widest (fewest vertices) to
the narrowest: `polys[i]` has error `widths[i]`, and `widths` is
non-increasing. `plt` is the palette of `P` and `len` its length. Access is
guarded by `lock`.
"""
mutable struct SimplificationHierarchy{D,T,M<:HierarchyMeasure}
    len::T
    P::PntSeq{D,T}
    plt::Palette{T}
    polys::Vector{PntSeq{D,T}}
    widths::Vector{Float64}
    lock::ReentrantLock
    measure::M
end

"""
    hierarchy_init(P, measure = hausdorff_measure(); lock = ReentrantLock())

Create the hierarchy of `P`: compute its palette and store the two-vertex
simplification as the initial level.
"""
function hierarchy_init(
    P::PntSeq{D,T},
    measure::HierarchyMeasure = hausdorff_measure();
    lock::ReentrantLock = ReentrantLock(),
) where {D,T<:AbstractFloat}
    h = SimplificationHierarchy(
        convert(T, geom_length(P)),
        P,
        palette(P, measure.err; rule = measure.rule),
        PntSeq{D,T}[],
        Float64[],
        lock,
        measure,
    )
    hierarchy_push!(h, PntSeq{D,T}([P[1], P[end]]), Float64(measure.err(P, 1:length(P))))
    return h
end

"""
    hierarchy_push!(h, Q, w)

Add the simplification `Q` of error `w`. A previous last level that is no
narrower and no smaller than `Q` is dropped, and the levels stay sorted by
decreasing width.
"""
function hierarchy_push!(
    h::SimplificationHierarchy{D,T},
    Q::PntSeq{D,T},
    w::Float64,
) where {D,T}
    @assert cardin(Q) > 1

    lock(h.lock) do
        if length(h.widths) > 0
            if last(h.widths) >= w && cardin(last(h.polys)) >= cardin(Q)
                pop!(h.polys)
                pop!(h.widths)
            end
        end
        Ps = h.polys
        Ws = h.widths
        push!(Ps, Q)
        push!(Ws, w)

        # sort in decreasing order...
        i = length(Ps)
        while i > 1 && Ws[i-1] < Ws[i]
            Ws[i-1], Ws[i] = Ws[i], Ws[i-1]
            Ps[i-1], Ps[i] = Ps[i], Ps[i-1]
            i -= 1
        end
    end
    return h
end

"""
    hierarchy_find(widths, w, resolution)

Index of the first level whose width is at most `w` and at least
`w / resolution`, or `-1` if the first level of width at most `w` is narrower
than that.
"""
function hierarchy_find(widths::Vector{Float64}, w, resolution)
    for i in eachindex(widths)
        if widths[i] <= w
            return w <= resolution * widths[i] ? i : -1
        end
    end
    return -1
end

function _push_target!(
    h::SimplificationHierarchy{D,T},
    w::Float64,
    P::PntSeq,
    lmt::Int,
    w_extra::Real,
) where {D,T}
    R, R_indices = h.measure.simplify(P, convert(T, w))
    cardin(R) > lmt && return -1.0

    new_w = h.measure.certify(P, R, R_indices) + w_extra
    hierarchy_push!(h, R, new_w)
    return new_w
end

"""
    hierarchy_approx(h, w, n_max = 25, resolution = 4.0)

An approximation of the curve with error at most about `w`, returned as
`(Q, width)`. If a stored level has width within a factor `resolution` below
`w`, the finest stored level with fewer than `n_max` vertices is returned
(`resolution` is used only to look levels up).
Otherwise a level is computed, stored in `h`, and returned: the new level is
built from a coarser level, or from the palette, of width about `w / 4`. This
factor is fixed.

The stored levels are read under `h.lock`, so concurrent calls and pushes are safe.
"""
function hierarchy_approx(
    h::SimplificationHierarchy{D,T},
    w::Float64,
    n_max::Int = 25,
    resolution::Float64 = 4.0,
) where {D,T}
    # The levels are read under the lock, but the work below is done outside it.
    found = lock(h.lock) do
        i = hierarchy_find(h.widths, w, resolution)
        i > 0 || return nothing
        len = length(h.polys)
        while i < len && cardin(h.polys[i+1]) < n_max
            i += 1
        end
        return (h.polys[i], h.widths[i])
    end
    found === nothing || return found

    new_w = w / 4.0
    coarse = lock(h.lock) do
        i = hierarchy_find(h.widths, new_w, resolution)
        i > 0 ? (h.polys[i], h.widths[i]) : nothing
    end
    if coarse === nothing
        Z, _, wZ = palette_extract(h.P, h.plt, new_w)
    else
        Z, wZ = coarse
    end

    X, X_indices = h.measure.simplify(Z, convert(T, max(w - wZ, 0.0)))
    w_out = h.measure.certify(Z, X, X_indices) + wZ

    hierarchy_push!(h, X, w_out)
    return X, w_out
end

"""
    hierarchy_build!(h, ratio = 4.0, lmt = 500)

Fill the hierarchy with levels whose widths decrease by roughly a factor
`ratio`, until a level has about as many vertices as the curve (at most `lmt`), or
until a round adds nothing new.
"""
function hierarchy_build!(
    h::SimplificationHierarchy{D,T},
    ratio::Float64 = 4.0,
    lmt::Int = 500,
) where {D,T}
    P = h.P
    lmt = min(round(Int, cardin(P) * 0.95), lmt)

    w = h.widths[1]
    for _ = 1:200
        cardin(last(h.polys)) >= lmt && break
        w_a = w / ratio

        Z, _, wZ = palette_extract(P, h.plt, w_a / 4.0)

        w_prev = w
        w = _push_target!(h, max(w_a - wZ, 0.0), Z, lmt, wZ)
        w < -0.5 && break
        # The next round would repeat this one exactly.
        w == w_prev && break
    end
    return h
end

"""
    hierarchy_print(h)

Print the number of vertices and width of every level.
"""
function hierarchy_print(h::SimplificationHierarchy)
    println("--------------------------------------------")
    println("|P|: ", cardin(h.P))
    for i in eachindex(h.widths)
        println("w[", i, "] :", h.widths[i], "  #: ", cardin(h.polys[i]))
    end
    println("----------------------------------------------\n")
end

end # module PntSeqHierarchy
