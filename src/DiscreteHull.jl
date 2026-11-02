###################################################################
### Discrete (lattice) convex hulls of planar convex bodies
#
# Implementation of
#   S. Har-Peled, "An output sensitive algorithm for discrete convex hulls",
#   Comput. Geom. Theory Appl. 10 (1998) 125-138.
#
# The discrete hull of a convex body C is ConvexHull(C ∩ Z²). The body is
# accessed only through a line-intersection (chord) oracle, `chord_range`.
# Each hull edge is found by a continued-fraction search (FindEdgeDirection)
# that performs O(log diam(C)) ray shootings, so the whole hull takes
# O(|C'| log diam(C)) oracle calls.

"""
    AbsLatticeBody

Abstract supertype of planar convex bodies whose discrete hull can be computed
by [`discrete_hull`](@ref). A subtype implements [`chord_range`](@ref).
"""
abstract type AbsLatticeBody end

"""
    chord_range(C::AbsLatticeBody, o::Point2I, d::Point2I)

Intersect `C` with the line `o + t*d` (`d ≠ 0`). Returns `nothing` when the
line misses `C`; otherwise returns `(ceil(t_lo), floor(t_hi))` as `Int`s, where
`[t_lo, t_hi]` is the closed parameter interval of the chord. The integers in
`lo:hi` are exactly the `t` with `o + t*d ∈ C`; the range is empty (`lo > hi`)
when the chord holds no integer parameter.

This is the only access to `C` used by [`discrete_hull`](@ref). Exact
implementations make the hull exact.
"""
function chord_range end

###################################################################
### Disk with rational center and rational squared radius (exact)

"""
    LatticeDisk(center, radius; den = 2^16)
    LatticeDisk(cx_num, cy_num, r2_num, den)

Closed disk with exact rational data: center `(cx_num, cy_num) / den` and
squared radius `r2_num / den^2`. In the first form, `Integer` and `Rational`
inputs are represented exactly; floating-point center and radius are rounded to
multiples of `1/den` first (the disk then is that rational disk).
[`chord_range`](@ref) is exact for this type.
"""
struct LatticeDisk <: AbsLatticeBody
    cx::Int128
    cy::Int128
    r2::Int128
    den::Int128
    function LatticeDisk(cx::Integer, cy::Integer, r2::Integer, den::Integer)
        den > 0 || throw(ArgumentError("den must be positive"))
        r2 >= 0 || throw(ArgumentError("squared radius must be non-negative"))
        return new(Int128(cx), Int128(cy), Int128(r2), Int128(den))
    end
end

_exact_rational(x::Union{Integer,Rational}, _) = Rational{BigInt}(x)
_exact_rational(x::AbstractFloat, den) = Rational{BigInt}(round(BigInt, BigFloat(x) * den), den)

function LatticeDisk(center, radius::Real; den::Integer = 2^16)
    cx = _exact_rational(center[1], den)
    cy = _exact_rational(center[2], den)
    r = _exact_rational(radius, den)
    r >= 0 || throw(ArgumentError("radius must be non-negative"))
    common = lcm(denominator(cx), denominator(cy), denominator(r))
    return LatticeDisk(
        Int128(numerator(cx * common)),
        Int128(numerator(cy * common)),
        Int128(numerator((r * common)^2)),
        Int128(common),
    )
end

"""
    radius(D::LatticeDisk)

Radius of `D` as a `Float64`.
"""
radius(D::LatticeDisk) = sqrt(Float64(D.r2)) / Float64(D.den)

"""
    center(D::LatticeDisk)

Center of `D` as a `Point2F`.
"""
center(D::LatticeDisk) = Point2F(Float64(D.cx) / D.den, Float64(D.cy) / D.den)

Base.in(p::Point{2,<:Integer}, D::LatticeDisk) =
    (D.den * p[1] - D.cx)^2 + (D.den * p[2] - D.cy)^2 <= D.r2

# Sign of b² - a*c with a floating-point filter and a BigInt fallback.
function _disc_sign(a::Int128, b::Int128, c::Int128)
    af, bf, cf = Float64(a), Float64(b), Float64(c)
    df = bf * bf - af * cf
    abs(df) > 1e-12 * (bf * bf + abs(af * cf)) && return sign(df)
    return sign(big(b)^2 - big(a) * big(c))
end

function chord_range(D::LatticeDisk, o::Point2I, d::Point2I)
    W1, W2 = D.den * o[1] - D.cx, D.den * o[2] - D.cy
    V1, V2 = D.den * d[1], D.den * d[2]
    a = V1 * V1 + V2 * V2
    b = V1 * W1 + V2 * W2
    c = W1 * W1 + W2 * W2 - D.r2
    _disc_sign(a, b, c) < 0 && return nothing
    # q(t) = a t² + 2 b t + c; P(t) = W + t V; q'(t)/2 = V⋅P(t).
    inside(t) = (W1 + t * V1)^2 + (W2 + t * V2)^2 <= D.r2
    vdotp(t) = V1 * (W1 + t * V1) + V2 * (W2 + t * V2)
    ge_lo(t) = inside(t) || vdotp(t) >= 0     # t ≥ t_lo
    le_hi(t) = inside(t) || vdotp(t) <= 0     # t ≤ t_hi
    sq = sqrt(max(Float64(b)^2 - Float64(a) * Float64(c), 0.0))
    lo = ceil(Int, (-Float64(b) - sq) / Float64(a))
    hi = floor(Int, (-Float64(b) + sq) / Float64(a))
    while ge_lo(lo - 1)
        lo -= 1
    end
    while !ge_lo(lo)
        lo += 1
    end
    while le_hi(hi + 1)
        hi += 1
    end
    while !le_hi(hi)
        hi -= 1
    end
    return (lo, hi)
end

# The lattice point closest to the center lies in D whenever D ∩ Z² ≠ ∅.
function _default_seed(D::LatticeDisk)
    p = Point2I(round(Int, D.cx // D.den), round(Int, D.cy // D.den))
    return p ∈ D ? p : nothing
end

###################################################################
### Convex polygon with rational vertices (exact)

"""
    ConvexPolygonBody(vertices; den = 2^16)

Convex polygon with exact rational vertices, given in counterclockwise order (at
least three, not all collinear). `Integer` and `Rational` coordinates are exact;
floating-point coordinates are rounded to multiples of `1/den`.
[`chord_range`](@ref) is exact. The lattice points of a thin polygon need not be
lattice-connected, so [`discrete_hull`](@ref) finds its first vertex by a row
scan over the polygon's y-extent.
"""
struct ConvexPolygonBody <: AbsLatticeBody
    vertices::Vector{Point{2,Int128}}   # den * (true vertices)
    den::Int128
    function ConvexPolygonBody(vs; den::Integer = 2^16)
        R = [(_exact_rational(v[1], den), _exact_rational(v[2], den)) for v in vs]
        common = reduce(lcm, (denominator(c) for v in R for c in v); init = BigInt(1))
        V = [Point{2,Int128}(Int128(numerator(v[1] * common)), Int128(numerator(v[2] * common)))
             for v in R]
        n = length(V)
        n >= 3 || throw(ArgumentError("need at least three vertices"))
        sum(_lat_cross(V[i], V[mod1(i + 1, n)]) for i = 1:n) > 0 ||
            throw(ArgumentError("vertices must be counterclockwise and not all collinear"))
        all(_lat_cross(V[mod1(i + 1, n)] - V[i], V[mod1(i + 2, n)] - V[mod1(i + 1, n)]) >= 0
            for i = 1:n) || throw(ArgumentError("polygon is not convex"))
        return new(V, Int128(common))
    end
end

"""
    vertices(P::ConvexPolygonBody)

Vertices of `P` as exact `Rational{Int128}` pairs.
"""
vertices(P::ConvexPolygonBody) = [(v[1] // P.den, v[2] // P.den) for v in P.vertices]

function chord_range(P::ConvexPolygonBody, o::Point2I, d::Point2I)
    V = P.vertices
    n = length(V)
    O = P.den * Point{2,Int128}(o)
    lo_num, lo_den = Int128(-1), Int128(0)   # -∞ as (-1)/0
    hi_num, hi_den = Int128(1), Int128(0)    # +∞ as 1/0
    for i = 1:n
        e = V[mod1(i + 1, n)] - V[i]
        α = P.den * _lat_cross(e, d)         # constraint: α t ≥ β
        β = -_lat_cross(e, O - V[i])
        if α > 0
            (lo_den == 0 || widemul(β, lo_den) > widemul(lo_num, α)) &&
                ((lo_num, lo_den) = (β, α))
        elseif α < 0
            (hi_den == 0 || widemul(-β, hi_den) < widemul(hi_num, -α)) &&
                ((hi_num, hi_den) = (-β, -α))
        elseif β > 0
            return nothing
        end
    end
    widemul(lo_num, hi_den) > widemul(hi_num, lo_den) && return nothing
    return (Int(cld(lo_num, lo_den)), Int(fld(hi_num, hi_den)))
end

function Base.in(p::Point{2,<:Integer}, P::ConvexPolygonBody)
    V, n = P.vertices, length(P.vertices)
    q = P.den * Point{2,Int128}(p)
    return all(_lat_cross(V[mod1(i + 1, n)] - V[i], q - V[i]) >= 0 for i = 1:n)
end

###################################################################
### Generic body given by a floating-point chord oracle

"""
    FloatChordBody(f)

Convex body given by a floating-point chord oracle: `f(o::Point2F, d::Point2F)`
returns `nothing` or `(t_lo, t_hi)` with `{t : o + t*d ∈ C} = [t_lo, t_hi]`.
The discrete hull is exact only up to the rounding of `f` near integer
parameters; use an exact body (e.g. [`LatticeDisk`](@ref)) when that matters.
"""
struct FloatChordBody{F} <: AbsLatticeBody
    f::F
end

function chord_range(C::FloatChordBody, o::Point2I, d::Point2I)
    r = C.f(Point2F(o), Point2F(d))
    r === nothing && return nothing
    return (ceil(Int, r[1]), floor(Int, r[2]))
end

###################################################################
### The algorithm

@inline _lat_cross(u, v) = Int128(u[1]) * v[2] - Int128(u[2]) * v[1]

# Does p + v lie in C?
function _offset_in(C, p::Point2I, v::Point2I)
    r = chord_range(C, p, v)
    return r !== nothing && r[1] <= 1 <= r[2]
end

# Closed quadrant j (angles [90j, 90j + 90]) as a unimodular basis (L, U).
const _QUADRANT_BASES = (
    (Point2I(1, 0), Point2I(0, 1)),
    (Point2I(0, 1), Point2I(-1, 0)),
    (Point2I(-1, 0), Point2I(0, -1)),
    (Point2I(0, -1), Point2I(1, 0)),
)

# Half-open quadrant index of a nonzero vector: angle in [90j, 90j + 90).
function _quadrant(e::Point2I)
    x, y = e[1], e[2]
    x > 0 && y >= 0 && return 0
    x <= 0 && y > 0 && return 1
    x < 0 && y <= 0 && return 2
    return 3
end

"""
    find_edge_direction(C, p, L, U)

FindEdgeDirection restricted to the closed cone spanned by the lattice basis
`(L, U)` (`det(L, U) = 1`, both in one closed quadrant). Returns the primitive
lattice vector `e` in this cone with `p + e ∈ C` of smallest angle from `L`, or
`nothing`. Requires `p ∈ C`.

The search mirrors the Euclidean algorithm on the unknown `e`: it alternates a
ray shooting from `L` in direction `U` (column step, even convergents) and from
`U` in direction `L` (row step, odd convergents). Invariants: `(L, U)` is a
unimodular basis, `p + L ∉ C`, no feasible vector lies at angle `≤ L`, and every
feasible vector at angle `≥ U` is no better than `best`. Each pair of steps at
least doubles the basis vectors, so O(log diam(C)) chord queries suffice.
"""
function find_edge_direction(C, p::Point2I, L::Point2I, U::Point2I)
    _offset_in(C, p, L) && return L
    best = _offset_in(C, p, U) ? U : nothing
    for _ = 1:256
        # Column step: candidates p + L + tU, t ≥ 0. Feasible vectors αL + βU
        # (α ≥ 1) have β/α in [t_lo, t_hi]; c = ⌈t_lo⌉ ≥ 1 since p + L ∉ C.
        r = chord_range(C, p + L, U)
        (r === nothing || r[2] < 0) && return best
        c = max(r[1], 1)
        newU = L + c * U
        c <= r[2] && (best = newU)
        L, U = L + (c - 1) * U, newU
        # Row step: candidates p + sL + U, s ≥ 0. Feasible vectors αL + βU
        # (β ≥ 1) have α/β in [s_lo, s_hi]; a = ⌊s_hi⌋.
        r = chord_range(C, p + U, L)
        (r === nothing || r[2] < 0) && return best
        a = r[2]
        newU = a * L + U
        max(r[1], 0) <= a && (best = newU)
        L, U = newU + L, newU
    end
    error("find_edge_direction did not terminate; inconsistent chord oracle?")
end

# Next counterclockwise hull vertex after the vertex p, searching the quadrants
# q, q+1, q+2 in order (q = quadrant of the previous edge direction).
function _next_hull_vertex(C, p::Point2I, q::Int)
    for k = 0:2
        j = mod(q + k, 4)
        L, U = _QUADRANT_BASES[j+1]
        e = find_edge_direction(C, p, L, U)
        e === nothing && continue
        return p + chord_range(C, p, e)[2] * e, e
    end
    return nothing
end

# Clockwise-most of the candidate directions from an arbitrary hull vertex p.
function _first_hull_edge(C, p::Point2I)
    cands = Point2I[]
    for (L, U) in _QUADRANT_BASES
        e = find_edge_direction(C, p, L, U)
        e === nothing || push!(cands, e)
    end
    isempty(cands) && return nothing
    i = findfirst(d -> all(_lat_cross(d, c) >= 0 for c in cands), cands)
    i === nothing && throw(ArgumentError("start point is not a vertex of the discrete hull"))
    e = cands[i]
    return p + chord_range(C, p, e)[2] * e, e
end

_row_has_lattice_point(C, y::Int) =
    (r = chord_range(C, Point2I(0, y), Point2I(1, 0)); r !== nothing && r[1] <= r[2])

"""
    discrete_hull_start(C, seed::Point2I)

Lowest-leftmost lattice point of `C` by unbounded and binary search over rows
(Lemma 3.10), using O(log diam(C)) chord queries. Requires `seed ∈ C` and that
`C ∩ Z²` be lattice-connected (true for disks).
"""
function discrete_hull_start(C, seed::Point2I)
    _row_has_lattice_point(C, seed[2]) || throw(ArgumentError("seed is not in C"))
    hi, step = seed[2], 1
    while _row_has_lattice_point(C, seed[2] - step)
        hi = seed[2] - step
        step *= 2
    end
    lo = seed[2] - step                       # row lo is empty, row hi is not
    while hi - lo > 1
        mid = (lo + hi) >> 1
        _row_has_lattice_point(C, mid) ? (hi = mid) : (lo = mid)
    end
    return Point2I(chord_range(C, Point2I(0, hi), Point2I(1, 0))[1], hi)
end

# Row scan over the y-extent, for bodies without a connectivity guarantee.
function _scan_start(P::ConvexPolygonBody)
    ys = [v[2] for v in P.vertices]
    for y = fld(minimum(ys), P.den):cld(maximum(ys), P.den)
        _row_has_lattice_point(P, Int(y)) &&
            return Point2I(chord_range(P, Point2I(0, y), Point2I(1, 0))[1], y)
    end
    return nothing
end

"""
    discrete_hull(C; seed = nothing, start = nothing) -> Vector{Point2I}

Vertices of the discrete hull `ConvexHull(C ∩ Z²)` of the planar convex body
`C`, in counterclockwise order starting at the lowest-leftmost vertex (no
collinear vertices). Returns an empty vector when `C ∩ Z² = ∅`.

Implements the output-sensitive algorithm of Har-Peled (1998): O(|C'| log
diam(C)) calls to [`chord_range`](@ref), where `|C'|` is the output size. For a
disk of radius `r` this is O(r^(2/3) log r).

The first vertex is found from `start` (any vertex of the discrete hull), or
else from `seed` (any lattice point of `C`; requires `C ∩ Z²` to be
lattice-connected) by [`discrete_hull_start`](@ref). [`LatticeDisk`](@ref) and
[`ConvexPolygonBody`](@ref) supply defaults.
"""
function discrete_hull(C::AbsLatticeBody; seed = nothing, start = nothing)
    if start === nothing
        if seed === nothing
            start = _default_start(C)
            start === nothing && return Point2I[]
        else
            start = discrete_hull_start(C, Point2I(seed[1], seed[2]))
        end
    end
    p0 = Point2I(start[1], start[2])
    r = chord_range(C, p0, Point2I(1, 0))
    (r !== nothing && r[1] <= 0 <= r[2]) || throw(ArgumentError("start is not in C"))
    hull = [p0]
    step = _first_hull_edge(C, p0)
    step === nothing && return hull
    p, e = step
    while p != p0
        push!(hull, p)
        length(hull) > 10^9 && error("discrete_hull did not close")
        p, e = _next_hull_vertex(C, p, _quadrant(e))
    end
    k = argmin([(v[2], v[1]) for v in hull])
    return circshift(hull, 1 - k)
end

function _default_start(D::LatticeDisk)
    s = _default_seed(D)
    return s === nothing ? nothing : discrete_hull_start(D, s)
end
_default_start(P::ConvexPolygonBody) = _scan_start(P)
_default_start(::AbsLatticeBody) =
    throw(ArgumentError("pass `seed` (a lattice point of C) or `start` (a hull vertex)"))

###################################################################
### Baseline: hull of the column extremes (the paper's "IG" algorithm)

"""
    discrete_hull_columns(C, xs) -> Vector{Point2I}

Baseline O(|xs|) discrete hull: takes the lowest and highest lattice points of
`C` on every vertical line `x ∈ xs` and runs an exact monotone-chain hull. The
result has the same format as [`discrete_hull`](@ref). `xs` must cover the
x-extent of `C`. For a `LatticeDisk` the range defaults to its x-extent.
"""
function discrete_hull_columns(C::AbsLatticeBody, xs::AbstractUnitRange{<:Integer})
    pts = Point2I[]
    for x in xs
        r = chord_range(C, Point2I(x, 0), Point2I(0, 1))
        (r === nothing || r[1] > r[2]) && continue
        push!(pts, Point2I(x, r[1]))
        r[2] > r[1] && push!(pts, Point2I(x, r[2]))
    end
    return lattice_convex_hull(pts)
end

function discrete_hull_columns(D::LatticeDisk)
    c, r = center(D), radius(D)
    return discrete_hull_columns(D, floor(Int, c[1] - r)-1:ceil(Int, c[1] + r)+1)
end

"""
    lattice_convex_hull(pts) -> Vector{Point2I}

Exact convex hull of integer points (Andrew's monotone chain, integer
orientation tests), counterclockwise from the lowest-leftmost point, without
collinear vertices.
"""
function lattice_convex_hull(pts::AbstractVector)
    P = unique!(sort!([Point2I(p[1], p[2]) for p in pts]; by = p -> (p[1], p[2])))
    length(P) <= 2 && return sort!(P; by = p -> (p[2], p[1]))
    H = Point2I[]
    for pass = 1:2
        start = length(H)
        for p in (pass == 1 ? P : Iterators.reverse(P))
            while length(H) >= start + 2 && _lat_cross(H[end] - H[end-1], p - H[end-1]) <= 0
                pop!(H)
            end
            push!(H, p)
        end
        pop!(H)
    end
    k = argmin([(v[2], v[1]) for v in H])
    return circshift(H, 1 - k)
end
