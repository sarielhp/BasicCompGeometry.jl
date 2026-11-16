# Ear clipping for small simple polygons (3 <= n <= 64): the fast path of
# `triangulate_ear_cut`. Adapted from the small-polygon kernel of Triage
# (https://github.com/sarielhp/triage, MIT license). Reflex vertices live in one `UInt64`
# mask, ears are tested lazily while walking the ring (not all up front), and
# point-in-triangle tests start with a bounding-box rejection (`_point_in_ccw_triangle`).

const _SMALL_EAR_CUT_MAX = 64

@inline _vertex_bit(i::Int) = UInt64(1) << (i - 1)

# One pass over the ring: the reflex mask (collinear counts as reflex, so it blocks
# conservatively), twice the signed area (the same sum, in the same order, as `twice_area`),
# and a check that no two consecutive vertices are equal.
function _small_classify(kernel::AbstractKernel, p, n::Int)
    origin = point(0, 0)
    area2 = orient(kernel, origin, p[1], p[2])
    reflex = UInt64(0)
    for i in 1:n
        prev, next = i == 1 ? n : i - 1, i == n ? 1 : i + 1
        p[i] == p[next] && _throw_duplicate_vertex(i, next)
        i > 1 && (area2 += orient(kernel, origin, p[i], p[next]))
        orient_sign(kernel, p[prev], p[i], p[next]) <= 0 && (reflex |= _vertex_bit(i))
    end
    return reflex, area2
end

# Fan from vertex 1 of a strictly convex polygon.
function _small_fan(n::Int)
    diagonals = Tuple{Int,Int}[(1, i) for i in 3:(n-1)]
    triangles = NTuple{3,Int}[(1, i, i + 1) for i in 2:(n-1)]
    return PolygonTriangulation(diagonals, triangles)
end

# Whether a reflex vertex other than the ear's corners lies in the closed ear triangle.
function _small_blocked(kernel::AbstractKernel, p, ear::NTuple{3,Int}, reflex::UInt64)
    a, b, c = ear
    bits = reflex & ~(_vertex_bit(a) | _vertex_bit(b) | _vertex_bit(c))
    pa, pb, pc = p[a], p[b], p[c]
    while bits != 0
        i = trailing_zeros(bits) + 1
        bits &= bits - 1
        _point_in_ccw_triangle(kernel, p[i], pa, pb, pc) && return true
    end
    return false
end

@inline _set_reflex_bit(reflex::UInt64, i::Int, convex::Bool) =
    convex ? reflex & ~_vertex_bit(i) : reflex | _vertex_bit(i)

function _small_clip(kernel::AbstractKernel, p, n::Int, reflex::UInt64)
    prev = Int[i == 1 ? n : i - 1 for i in 1:n]
    next = Int[i == n ? 1 : i + 1 for i in 1:n]
    diagonals = Tuple{Int,Int}[]
    sizehint!(diagonals, n - 3)
    triangles = NTuple{3,Int}[]
    sizehint!(triangles, n - 2)
    current, remaining, misses = 1, n, 0
    while remaining > 3
        a, c = prev[current], next[current]
        if reflex & _vertex_bit(current) == 0 &&
           !_small_blocked(kernel, p, (a, current, c), reflex)
            push!(triangles, (a, current, c))
            _is_boundary_edge(a, c, n) || push!(diagonals, minmax(a, c))
            next[a], prev[c] = c, a
            # Recompute both neighbors' bits in both directions: on non-simple or numerically
            # degenerate input a convex vertex can turn reflex, and a stale convex bit would
            # let a clockwise ear through.
            pa, pc = p[a], p[c]
            reflex = _set_reflex_bit(reflex, a, orient_sign(kernel, p[prev[a]], pa, pc) > 0)
            reflex = _set_reflex_bit(reflex, c, orient_sign(kernel, pa, pc, p[next[c]]) > 0)
            remaining -= 1
            misses = 0
        else
            misses += 1
            misses < remaining || _throw_no_ear(remaining)
        end
        current = c
    end
    b = next[current]
    c = next[b]
    orient_sign(kernel, p[current], p[b], p[c]) >= 0 || _throw_bad_last_triangle()
    push!(triangles, (current, b, c))
    return PolygonTriangulation(diagonals, triangles)
end

# Fast path of `triangulate_ear_cut` for `3 <= n <= _SMALL_EAR_CUT_MAX`.
function _small_ear_cut(kernel::AbstractKernel, vertices)
    n = length(vertices)
    reflex, area2 = _small_classify(kernel, vertices, n)
    area2 > 0 || _throw_not_ccw(area2)
    n == 3 && return PolygonTriangulation(Tuple{Int,Int}[], NTuple{3,Int}[(1, 2, 3)])
    reflex == 0 && return _small_fan(n)
    return _small_clip(kernel, vertices, n, reflex)
end
