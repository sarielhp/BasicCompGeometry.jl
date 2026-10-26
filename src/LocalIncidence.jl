# Exact vertical-ray contacts and angular-sector predicates.

"""Abstract supertype of [`NoHit`](@ref), [`PointHit`](@ref), [`OverlapHit`](@ref)."""
abstract type RayHit end

"""No intersection of the closed ray with the closed segment."""
struct NoHit <: RayHit end

"""Singleton intersection `point` at exact parameter `t >= 0`."""
struct PointHit{P,S} <: RayHit
    point::P
    t::S
end

"""Nontrivial collinear intersection: parameters `0 <= t_lo < t_hi`, both included."""
struct OverlapHit{S} <: RayHit
    t_lo::S
    t_hi::S
end

Base.:(==)(::NoHit, ::NoHit) = true
Base.:(==)(a::PointHit, b::PointHit) = a.t == b.t && all(a.point .== b.point)
Base.:(==)(a::OverlapHit, b::OverlapHit) = a.t_lo == b.t_lo && a.t_hi == b.t_hi

function _check_direction(direction)
    direction in (-1, 1) || throw(ArgumentError("direction must be -1 or 1"))
    return Int(direction)
end

"""
    vertical_ray_segment(p, direction, a, b)

Intersect the closed ray `p + (0,s)t`, `t >= 0`, `s = direction in {-1,1}`, with
the closed segment `ab` (vertical and zero-length segments allowed). Returns
[`NoHit`](@ref), [`PointHit`](@ref) or [`OverlapHit`](@ref) (never a degenerate
overlap). Points are tuples or 2-vectors with `Integer`/`Rational` coordinates;
all arithmetic is exact (widened to `BigInt`), with O(1) operations. Returned
`t` values and constructed coordinates are `Rational{BigInt}`.
"""
function vertical_ray_segment(p, direction::Integer, a, b)
    s = _check_direction(direction)
    px, py = _exact(_cx(p)), _exact(_cy(p))
    ax, ay, bx, by = _exact(_cx(a)), _exact(_cy(a)), _exact(_cx(b)), _exact(_cy(b))
    if ax == bx
        return _vertical_overlap(p, s, px, py, ax, ay, by)
    end
    if ax > bx
        ax, ay, bx, by = bx, by, ax, ay
    end
    (px < ax || px > bx) && return NoHit()
    y = ay + (px - ax) * (by - ay) // (bx - ax)
    t = s * (y - py)
    t < 0 && return NoHit()
    return PointHit(_mkpoint(p, Rational{BigInt}(px), Rational{BigInt}(y)), Rational{BigInt}(t))
end

function _vertical_overlap(p, s, px, py, ax, ay, by)
    px == ax || return NoHit()
    ylo, yhi = minmax(ay, by)
    tlo, thi = s == 1 ? (ylo - py, yhi - py) : (py - yhi, py - ylo)
    tlo = max(tlo, zero(tlo))
    thi < 0 && return NoHit()
    if tlo == thi
        point = _mkpoint(p, Rational{BigInt}(px), Rational{BigInt}(py + s * tlo))
        return PointHit(point, Rational{BigInt}(tlo))
    end
    return OverlapHit(Rational{BigInt}(tlo), Rational{BigInt}(thi))
end

"""
    first_positive_vertical_contact(segments, p, direction)

Scan records `(a, b, carrier)` of **nonvertical** segments (a vertical record
throws `ArgumentError`) and return `(point, t, carrier, record_index)` for the
minimum `t > 0` contact, or `nothing`. Ties in `t` go to the smallest record
index; `carrier` is returned untouched. Endpoint contacts count, an origin-only
contact (`t == 0`) does not. O(k) exact operations, O(1) extra space.
"""
function first_positive_vertical_contact(segments, p, direction::Integer)
    _check_direction(direction)
    best = nothing
    for (i, rec) in enumerate(segments)
        a, b, carrier = rec
        _cx(a) == _cx(b) &&
            throw(ArgumentError("record $i is a vertical segment; only nonvertical records are supported"))
        hit = vertical_ray_segment(p, direction, a, b)
        hit isa PointHit && hit.t > 0 || continue
        if best === nothing || hit.t < best.t
            best = (; point = hit.point, t = hit.t, carrier, record_index = i)
        end
    end
    return best
end

_cross(a, b) = _exact(_cx(a)) * _exact(_cy(b)) - _exact(_cy(a)) * _exact(_cx(b))
_dot(a, b) = _exact(_cx(a)) * _exact(_cx(b)) + _exact(_cy(a)) * _exact(_cy(b))

function _check_nonzero(v, name)
    (iszero(_cx(v)) && iszero(_cy(v))) && throw(ArgumentError("vector $name must be nonzero"))
    return nothing
end

function _check_coincident(coincident)
    coincident in (:full, :empty) ||
        throw(ArgumentError("coincident must be :full or :empty, got $(repr(coincident))"))
    return nothing
end

# Shared decision from signs: c = cross(v,u), dvu = dot(v,u), L = cross(v,d),
# R = cross(d,u), dvd = dot(v,d) (only used for a same-direction full turn).
function _sector_decision(c, dvu, L, R, dvd, coincident)
    c > 0 && return L > 0 && R > 0
    c < 0 && return L > 0 || R > 0
    dvu < 0 && return L > 0
    coincident == :empty && return false
    return L != 0 || dvd < 0
end

"""
    open_sector_contains(v, u, d; coincident=:full)

Whether direction `d` lies in the open counterclockwise sector from ray `v` to
ray `u` (all nonzero vectors at a common base point; reflex sectors included,
endpoints excluded). If `u` and `v` have the same direction, `coincident=:full`
(slit tip) means a full turn excluding the common ray and `:empty` means an empty
sector. Opposite rays give the open ccw half-plane. Exact, O(1) operations.
"""
function open_sector_contains(v, u, d; coincident::Symbol = :full)
    _check_coincident(coincident)
    foreach(x -> _check_nonzero(x...), ((v, "v"), (u, "u"), (d, "d")))
    return _sector_decision(
        sign(_cross(v, u)), sign(_dot(v, u)), sign(_cross(v, d)), sign(_cross(d, u)),
        sign(_dot(v, d)), coincident,
    )
end

# Sign of A + epsilon*B for all sufficiently small epsilon > 0.
_lexsign(A, B) = iszero(A) ? sign(B) : sign(A)

"""
    incoming_sector_contains(v, u, direction; coincident=:full)

[`open_sector_contains`](@ref) for the symbolic incoming direction
`d_eps = (-eps, -direction)`, `direction in {-1,1}`, as `eps -> 0+` (the incoming
direction tilted toward smaller x). Signs of degree-one expressions are taken
lexicographically (`A + eps*B` has the sign of `A`, else of `B`), so no numeric
epsilon is used. With nonzero `u`, `v` the symbolic direction never lies on a
boundary ray, so a full-turn sector contains it. Exact, O(1) operations.
"""
function incoming_sector_contains(v, u, direction::Integer; coincident::Symbol = :full)
    s = _check_direction(direction)
    _check_coincident(coincident)
    _check_nonzero(v, "v")
    _check_nonzero(u, "u")
    L = _lexsign(-s * _exact(_cx(v)), _exact(_cy(v)))
    R = _lexsign(s * _exact(_cx(u)), -_exact(_cy(u)))
    # dot(v, d_eps) only matters for a same-direction full turn, where L != 0 already.
    return _sector_decision(sign(_cross(v, u)), sign(_dot(v, u)), L, R, 0, coincident)
end
