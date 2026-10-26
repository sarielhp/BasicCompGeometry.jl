# Shared static cyclic-chain geometry on exact planar points.
#
# Occurrence convention: vertex occurrence i (0 <= i <= n) is P[mod(i,n)+1];
# occurrences 0 and n share coordinates but remain distinct identities.
# Original edge i (1 <= i <= n) joins occurrence i-1 to occurrence i.

const _ExactScalar = Union{Integer,Rational}

_cx(p::Tuple) = p[1]
_cy(p::Tuple) = p[2]
_cx(p::AbstractVector) = p[1]
_cy(p::AbstractVector) = p[2]

# Exact promotion: fixed-width integers/rationals are widened to BigInt-based
# scalars before any arithmetic, so predicates and constructions cannot overflow.
_exact(x::Integer) = BigInt(x)
_exact(x::Rational) = Rational{BigInt}(x)
_exact(x) = throw(ArgumentError("coordinates must be exact (Integer or Rational), got $(typeof(x))"))

function _check_exact_point(p)
    (p isa Tuple || p isa AbstractVector) && length(p) == 2 ||
        throw(ArgumentError("points must be planar tuples or 2-vectors"))
    for c in (_cx(p), _cy(p))
        c isa _ExactScalar ||
            throw(ArgumentError("coordinates must be exact (Integer or Rational), got $(typeof(c))"))
        (c isa Rational && iszero(denominator(c))) &&
            throw(ArgumentError("coordinates must be finite"))
    end
    return nothing
end

# Build a point of the same family as `like` from exact coordinates.
_mkpoint(::Tuple, x, y) = (x, y)
_mkpoint(::AbstractVector, x, y) = SVector{2}(x, y)

_lex_less(a, b) = (_cx(a) < _cx(b)) || (_cx(a) == _cx(b) && _cy(a) < _cy(b))
_lex_greater(a, b) = _lex_less(b, a)
_yx_less(a, b) = (_cy(a) < _cy(b)) || (_cy(a) == _cy(b) && _cx(a) < _cx(b))
_yx_greater(a, b) = _yx_less(b, a)

"""
    ChainIndex

Static index over a cyclic point chain, built by [`build_chain_index`](@ref).
"""
struct ChainIndex{V,A,B,C,D}
    points::V
    xmin::A
    xmax::B
    ymin::C
    ymax::D
end

"""
    ChainWindow

O(1) view `[lo, hi]` of occurrences of a [`ChainIndex`](@ref); it holds only the
index and its two endpoints. Use [`nedges`](@ref) for its number of edges.
"""
struct ChainWindow{I<:ChainIndex}
    index::I
    lo::Int
    hi::Int
end

"""
    build_chain_index(points)

Index the nonempty vector `points` of exact planar points (tuples or 2-vectors
with `Integer`/`Rational` coordinates; repeated points allowed, no closing
duplicate). Occurrence `i`, `0 <= i <= n`, is `points[mod(i,n)+1]`; occurrences
`0` and `n` are distinct identities. Preprocessing is O(n) time and words (four
linear [`build_rmq`](@ref) indexes); the points are borrowed read-only and must
not be mutated while the index or its windows are in use.
Throws `ArgumentError` for empty or inexact input.
"""
function build_chain_index(points::AbstractVector)
    Base.require_one_based_indexing(points)
    isempty(points) && throw(ArgumentError("a chain index needs at least one point"))
    foreach(_check_exact_point, points)
    return ChainIndex(
        points,
        build_rmq(points; less = _lex_less),
        build_rmq(points; less = _lex_greater),
        build_rmq(points; less = _yx_less),
        build_rmq(points; less = _yx_greater),
    )
end

Base.length(index::ChainIndex) = length(index.points)

"""
    chain_window(index, lo, hi)

Window of occurrences `lo:hi`, `0 <= lo <= hi <= n`, with `hi-lo` edges. O(1).
"""
function chain_window(index::ChainIndex, lo::Integer, hi::Integer)
    0 <= lo <= hi <= length(index) || throw(BoundsError(index.points, lo:hi))
    return ChainWindow(index, Int(lo), Int(hi))
end

"""
    nedges(w)

Number of edges `hi-lo` of a window (it has `nedges(w)+1` vertex occurrences).
"""
nedges(w::ChainWindow) = w.hi - w.lo
Base.length(w::ChainWindow) = nedges(w)

"""
    vertex(index, i)
    vertex(w, i)

Point at occurrence `i` (`0 <= i <= n`, or `lo <= i <= hi` for a window).
"""
function vertex(index::ChainIndex, i::Integer)
    0 <= i <= length(index) || throw(BoundsError(index.points, i))
    return index.points[mod(i, length(index))+1]
end

function vertex(w::ChainWindow, i::Integer)
    w.lo <= i <= w.hi || throw(BoundsError(w, i))
    return vertex(w.index, i)
end

"""
    edge(w, i; direction=1)

Endpoints of original edge `i`, `lo < i <= hi`, as `(occurrence i-1, occurrence i)`;
`direction=-1` reverses them.
"""
function edge(w::ChainWindow, i::Integer; direction::Integer = 1)
    w.lo < i <= w.hi || throw(BoundsError(w, i))
    direction in (-1, 1) || throw(ArgumentError("edge direction must be -1 or 1"))
    a, b = vertex(w, i - 1), vertex(w, i)
    return direction == 1 ? (a, b) : (b, a)
end

# Occurrence in `w` of the extreme under `table` (rightmost on equivalent keys).
# The closing occurrence n is point 1 again but wins ties by being later, so a
# window reaching n is never turned into an index interval starting at 0.
function _extreme_occurrence(w::ChainWindow, table, less)
    n = length(w.index)
    w.lo == n && return n
    best = range_argmin(table, w.lo + 1, min(w.hi, n - 1) + 1) - 1
    if w.hi == n && !less(vertex(w, best), vertex(w, n))
        return n
    end
    return best
end

"""
    leftmost(w)

`(occurrence, point)` minimizing `(x,y)` lexicographically over the window;
equal points choose the greatest occurrence. O(1).
"""
function leftmost(w::ChainWindow)
    occurrence = _extreme_occurrence(w, w.index.xmin, _lex_less)
    return (; occurrence, point = vertex(w, occurrence))
end

"""
    bounds(w)

Exact `(xmin, xmax, ymin, ymax)` over all vertices (hence edges) of the window.
O(1) comparisons.
"""
function bounds(w::ChainWindow)
    i = w.index
    return (;
        xmin = _cx(vertex(w, _extreme_occurrence(w, i.xmin, _lex_less))),
        xmax = _cx(vertex(w, _extreme_occurrence(w, i.xmax, _lex_greater))),
        ymin = _cy(vertex(w, _extreme_occurrence(w, i.ymin, _yx_less))),
        ymax = _cy(vertex(w, _extreme_occurrence(w, i.ymax, _yx_greater))),
    )
end

"""
    point_at(w, t, j)

Exact point at parameter `t` (Integer/Rational) with supplied edge position `j`:
either `lo <= j < hi` and `j <= t <= j+1` (interpolating occurrences `j`, `j+1`
with fraction `t-j`), or `j == t == hi` (occurrence `hi`). O(1) exact operations.
Constructed coordinates are `Rational{BigInt}`; fixed-width inputs are widened
first, so nothing overflows. Throws `BoundsError`/`ArgumentError` otherwise.
"""
function point_at(w::ChainWindow, t::_ExactScalar, j::Integer)
    (t isa Rational && iszero(denominator(t))) && throw(ArgumentError("t must be finite"))
    if j == t == w.hi
        p = vertex(w, w.hi)
        return _mkpoint(p, Rational{BigInt}(_exact(_cx(p))), Rational{BigInt}(_exact(_cy(p))))
    end
    w.lo <= j < w.hi || throw(BoundsError(w, j))
    j <= t <= j + 1 || throw(ArgumentError("parameter t=$t is not on edge position j=$j"))
    a, b = vertex(w, j), vertex(w, j + 1)
    f = _exact(t) - j
    ax, ay, bx, by = _exact(_cx(a)), _exact(_cy(a)), _exact(_cx(b)), _exact(_cy(b))
    return _mkpoint(a, Rational{BigInt}(ax + f * (bx - ax)), Rational{BigInt}(ay + f * (by - ay)))
end

"""
    point_at(w, t)

Convenience form that finds `j` itself by `floor` of `t`; this costs
O(log(hi-lo+2)) comparisons for arbitrary rationals (not O(1)). Prefer
`point_at(w, t, j)` when the edge position is known.
"""
function point_at(w::ChainWindow, t::_ExactScalar)
    (t isa Rational && iszero(denominator(t))) && throw(ArgumentError("t must be finite"))
    w.lo <= t <= w.hi || throw(BoundsError(w, t))
    j = t == w.hi ? w.hi : Int(floor(t))
    return point_at(w, t, j)
end
