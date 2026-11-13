"""
    PolygonTriangulation

Result of triangulating a simple counterclockwise polygon.
Fields:
- `diagonals`: `k - 3` distinct internal diagonals as index pairs `(i, j)` with `i < j`.
- `triangles`: `k - 2` counterclockwise triangles as index triples `(i, j, m)`.
"""
struct PolygonTriangulation
    diagonals::Vector{Tuple{Int,Int}}
    triangles::Vector{NTuple{3,Int}}
end

Base.iterate(res::PolygonTriangulation) = (res.diagonals, 1)
Base.iterate(res::PolygonTriangulation, state) = state == 1 ? (res.triangles, 2) : nothing
Base.length(::PolygonTriangulation) = 2

"""
    MonotoneTriangulation

Compatibility alias for [`PolygonTriangulation`](@ref).
"""
const MonotoneTriangulation = PolygonTriangulation

@inline function _is_boundary_edge(i::Int, j::Int, k::Int)
    d = abs(i - j)
    return d == 1 || d == k - 1
end

"""
    monotone_chain_indices(p; kernel=FloatingKernel())

Return `(upper, lower)` where `upper` and `lower` are vectors of 1-based vertex indices
tracing the upper and lower boundary chains of counterclockwise x-monotone polygon `p`
from the leftmost vertex to the rightmost vertex.
"""
function monotone_chain_indices(p; kernel::AbstractKernel = FloatingKernel())
    k = length(p)
    k >= 3 || throw(ArgumentError("polygon must have at least 3 vertices; got k = $k"))

    lo = 1; hi = 1
    for i in 2:k
        p[i].x < p[lo].x && (lo = i)
        p[i].x > p[hi].x && (hi = i)
    end

    c_minus = Int[lo]; curr = lo
    while curr != hi
        curr = mod1(curr - 1, k)
        push!(c_minus, curr)
    end

    c_plus = Int[lo]; curr = lo
    while curr != hi
        curr = mod1(curr + 1, k)
        push!(c_plus, curr)
    end

    # In a counterclockwise polygon, step -1 is the upper chain if orient_sign(lo, plus, minus) > 0
    if orient_sign(kernel, p[lo], p[c_plus[2]], p[c_minus[2]]) > 0
        return (c_minus, c_plus)
    else
        return (c_plus, c_minus)
    end
end

"""
    monotone_chains(p; kernel=FloatingKernel())

Return `(upper, lower)` where `upper` and `lower` are vectors of `Point`s representing
the upper and lower boundary chains of counterclockwise x-monotone polygon `p`
from the leftmost vertex to the rightmost vertex.
"""
function monotone_chains(p; kernel::AbstractKernel = FloatingKernel())
    upper_idx, lower_idx = monotone_chain_indices(p; kernel = kernel)
    return ([p[i] for i in upper_idx], [p[i] for i in lower_idx])
end

"""
    merge_chains(upper, lower)

Merge two point sequences `upper` and `lower` sorted by x-coordinate into a single sequence
in linear time `O(k)`, eliminating duplicate endpoints.
"""
function merge_chains(upper::AbstractVector{T}, lower::AbstractVector{T}) where {T<:Point}
    order = T[]
    sizehint!(order, length(upper) + length(lower))
    i = 1; j = 1
    while i <= length(upper) || j <= length(lower)
        if j > length(lower) || (i <= length(upper) && upper[i].x <= lower[j].x)
            pt = upper[i]; i += 1
        else
            pt = lower[j]; j += 1
        end
        (isempty(order) || last(order) != pt) && push!(order, pt)
    end
    return order
end

"""
    triangulate_monotone(vertices; kernel=FloatingKernel()) -> MonotoneTriangulation

Triangulate a simple counterclockwise x-monotone polygon with `k >= 3` vertices in linear time `O(k)`.
Returns a `MonotoneTriangulation` with:
- `result.diagonals`: `k - 3` distinct internal diagonals as index pairs `(i, j)` with `i < j`.
- `result.triangles`: `k - 2` counterclockwise triangles as index triples `(i, j, m)`.

# Preconditions
- `vertices` must form a simple counterclockwise x-monotone polygon with `k >= 3` vertices.
- Original x-coordinates are distinct, and no three supplied vertices are collinear.
"""
function triangulate_monotone(
    vertices;
    kernel::AbstractKernel = FloatingKernel(),
)
    k = length(vertices)
    k >= 3 || throw(ArgumentError("polygon must have at least 3 vertices; got k = $k"))

    if k == 3
        t = orient_sign(kernel, vertices[1], vertices[2], vertices[3]) > 0 ? (1, 2, 3) : (1, 3, 2)
        return MonotoneTriangulation(Tuple{Int,Int}[], NTuple{3,Int}[t])
    end

    upper, lower = monotone_chain_indices(vertices; kernel = kernel)

    is_upper = falses(k)
    for v in @view upper[2:end-1]
        is_upper[v] = true
    end

    # Merge chains by increasing x in linear time O(k)
    order = Int[]
    sizehint!(order, k)
    i = 1; j = 1
    while i <= length(upper) || j <= length(lower)
        if j > length(lower) || (i <= length(upper) && vertices[upper[i]].x <= vertices[lower[j]].x)
            v = upper[i]; i += 1
        else
            v = lower[j]; j += 1
        end
        (isempty(order) || last(order) != v) && push!(order, v)
    end

    diagonals = Tuple{Int,Int}[]
    sizehint!(diagonals, k - 3)
    triangles = NTuple{3,Int}[]
    sizehint!(triangles, k - 2)

    @inline function emit_diagonal!(a::Int, b::Int)
        if a != b && !_is_boundary_edge(a, b, k)
            push!(diagonals, a < b ? (a, b) : (b, a))
        end
    end

    @inline function make_ccw(i1::Int, i2::Int, i3::Int)
        if orient_sign(kernel, vertices[i1], vertices[i2], vertices[i3]) < 0
            return (i1, i3, i2)
        else
            return (i1, i2, i3)
        end
    end

    stack = Int[order[1], order[2]]
    for j in 3:k-1
        v = order[j]
        if is_upper[v] != is_upper[last(stack)]
            while length(stack) >= 2
                s_top = pop!(stack)
                s_next = last(stack)
                push!(triangles, make_ccw(v, s_top, s_next))
                emit_diagonal!(v, s_top)
            end
            pop!(stack)
            stack = Int[order[j-1], v]
        else
            q = pop!(stack)
            while !isempty(stack)
                a = last(stack)
                turn = orient(kernel, vertices[a], vertices[q], vertices[v])
                (is_upper[v] ? turn < 0 : turn > 0) || break
                push!(triangles, make_ccw(a, q, v))
                emit_diagonal!(a, v)
                q = pop!(stack)
            end
            push!(stack, q, v)
        end
    end

    v_last = last(order)
    while length(stack) >= 2
        s_top = pop!(stack)
        s_next = last(stack)
        push!(triangles, make_ccw(v_last, s_top, s_next))
        emit_diagonal!(v_last, s_top)
    end

    return MonotoneTriangulation(diagonals, triangles)
end

"""
    stack_diagonals(vertices; kernel=FloatingKernel())

Return the internal diagonals of a counterclockwise x-monotone polygon as index pairs `(i, j)` with `i < j`.
"""
@inline function stack_diagonals(vertices; kernel::AbstractKernel = FloatingKernel())
    return triangulate_monotone(vertices; kernel = kernel).diagonals
end

"""
    is_x_monotone(p) -> Bool

Return `true` if simple polygon `p` is strictly x-monotone (i.e. any vertical line
intersects `p` in at most a single connected component with strictly monotonic x-coordinates).
"""
function is_x_monotone(p)
    n = length(p)
    n >= 3 || return false

    lo = 1; hi = 1
    for i in 2:n
        p[i].x < p[lo].x && (lo = i)
        p[i].x > p[hi].x && (hi = i)
    end
    lo == hi && return false

    curr = lo
    while curr != hi
        nxt = mod1(curr + 1, n)
        p[nxt].x > p[curr].x || return false
        curr = nxt
    end

    curr = lo
    while curr != hi
        nxt = mod1(curr - 1, n)
        p[nxt].x > p[curr].x || return false
        curr = nxt
    end

    return true
end

"""
    is_y_monotone(p) -> Bool

Return `true` if simple polygon `p` is strictly y-monotone (i.e. any horizontal line
intersects `p` in at most a single connected component with strictly monotonic y-coordinates).
"""
function is_y_monotone(p)
    n = length(p)
    n >= 3 || return false

    lo = 1; hi = 1
    for i in 2:n
        p[i].y < p[lo].y && (lo = i)
        p[i].y > p[hi].y && (hi = i)
    end
    lo == hi && return false

    curr = lo
    while curr != hi
        nxt = mod1(curr + 1, n)
        p[nxt].y > p[curr].y || return false
        curr = nxt
    end

    curr = lo
    while curr != hi
        nxt = mod1(curr - 1, n)
        p[nxt].y > p[curr].y || return false
        curr = nxt
    end

    return true
end

@inline function _point_in_ccw_triangle(kernel::AbstractKernel, p::Point, a::Point, b::Point, c::Point)
    return orient_sign(kernel, a, b, p) >= 0 &&
           orient_sign(kernel, b, c, p) >= 0 &&
           orient_sign(kernel, c, a, p) >= 0
end

abstract type _AbstractReflexIndex end

mutable struct _PlainReflexIndex <: _AbstractReflexIndex
    is_reflex::BitVector
    r_prev::Vector{Int}
    r_next::Vector{Int}
    r_head::Int
end

function _build_plain_reflex_index(is_reflex::BitVector)
    n = length(is_reflex)
    r_prev = zeros(Int, n)
    r_next = zeros(Int, n)
    r_head = 0
    @inbounds for i in n:-1:1
        if is_reflex[i]
            r_next[i] = r_head
            r_head != 0 && (r_prev[r_head] = i)
            r_head = i
        end
    end
    return _PlainReflexIndex(is_reflex, r_prev, r_next, r_head)
end

@inline _is_reflex_vertex(idx::_AbstractReflexIndex, i::Int)::Bool = idx.is_reflex[i]

@inline function _add_reflex!(idx::_PlainReflexIndex, i::Int)
    idx.is_reflex[i] && return nothing
    idx.is_reflex[i] = true
    idx.r_prev[i] = 0
    idx.r_next[i] = idx.r_head
    idx.r_head != 0 && (idx.r_prev[idx.r_head] = i)
    idx.r_head = i
    return nothing
end

@inline function _remove_reflex!(idx::_PlainReflexIndex, i::Int)
    idx.is_reflex[i] || return nothing
    idx.is_reflex[i] = false
    rp = idx.r_prev[i]
    rn = idx.r_next[i]
    if rp != 0
        idx.r_next[rp] = rn
    else
        idx.r_head = rn
    end
    rn != 0 && (idx.r_prev[rn] = rp)
    idx.r_prev[i] = 0
    idx.r_next[i] = 0
    return nothing
end

function _ear_blocked(
    idx::_PlainReflexIndex,
    poly,
    kernel::AbstractKernel,
    pr::Int,
    ear::Int,
    nx::Int,
    a::Point,
    b::Point,
    c::Point,
)::Bool
    curr = idx.r_head
    while curr != 0
        if curr != pr && curr != ear && curr != nx
            _point_in_ccw_triangle(kernel, poly[curr], a, b, c) && return true
        end
        curr = idx.r_next[curr]
    end
    return false
end

mutable struct _EarClipState{K<:AbstractKernel,P,R<:_AbstractReflexIndex}
    kernel::K
    poly::P
    prev_v::Vector{Int}
    next_v::Vector{Int}
    is_active::BitVector
    rindex::R
end

@inline function _init_ear_clip_ring(n::Int)
    prev_v = Vector{Int}(undef, n)
    next_v = Vector{Int}(undef, n)
    @inbounds for i in 1:n
        prev_v[i] = i == 1 ? n : i - 1
        next_v[i] = i == n ? 1 : i + 1
    end
    return prev_v, next_v, trues(n)
end

@inline function _update_reflex!(st::_EarClipState, i::Int)
    if orient_sign(st.kernel, st.poly[st.prev_v[i]], st.poly[i], st.poly[st.next_v[i]]) <= 0
        _add_reflex!(st.rindex, i)
    else
        _remove_reflex!(st.rindex, i)
    end
    return nothing
end

@inline function _is_ear_vertex(st::_EarClipState, i::Int)::Bool
    _is_reflex_vertex(st.rindex, i) && return false
    pr = st.prev_v[i]
    nx = st.next_v[i]
    poly = st.poly
    return !_ear_blocked(st.rindex, poly, st.kernel, pr, i, nx, poly[pr], poly[i], poly[nx])
end

function _find_fallback_ear(st::_EarClipState, hint::Int)::Int
    start_v = st.is_active[hint] ? hint : findfirst(st.is_active)
    start_v === nothing && return 0
    curr = start_v
    while true
        _is_ear_vertex(st, curr) && return curr
        curr = st.next_v[curr]
        curr == start_v && return 0
    end
end

# Validate a simple-polygon input: `n >= 3`, positive (counterclockwise) signed area and,
# when `check_duplicates`, pairwise distinct vertices (O(n log n) sort). Returns `n`.
function _validate_ccw_polygon(kernel::AbstractKernel, vertices; check_duplicates::Bool = true)
    n = length(vertices)
    n >= 3 || throw(ArgumentError("polygon must have at least 3 vertices; got n = $n"))

    area2 = twice_area(vertices; kernel = kernel)
    area2 > 0 || throw(
        ArgumentError(
            "polygon must be oriented counterclockwise with positive area; got twice_area = $area2",
        ),
    )
    check_duplicates && _check_distinct_vertices(kernel, vertices)
    return n
end

# Duplicate detection compares raw coordinates (exact equality), independent of the kernel,
# so it adds no predicate calls.
function _check_distinct_vertices(kernel::AbstractKernel, vertices)
    n = length(vertices)
    order = sortperm([(vertices[i][1], vertices[i][2]) for i in 1:n]; alg = MergeSort)
    for i in 1:n-1
        a = vertices[order[i]]
        b = vertices[order[i + 1]]
        (a[1] == b[1] && a[2] == b[2]) && throw(
            ArgumentError("polygon contains duplicate vertex $(order[i + 1]) (equal to vertex $(order[i]))"),
        )
    end
    return nothing
end

@noinline function _throw_no_ear(active_count::Int)
    throw(
        ArgumentError(
            "ear clipping found no ear with $active_count vertices remaining; the polygon is not " *
            "simple, or is numerically degenerate for the kernel (try kernel = ExactKernel() or algorithm = :sweep)",
        ),
    )
end

# Ear queue. It is seeded with the polygon's initial ears (verified up front, in index
# order). When an ear is clipped, each of its two neighbors is (re)appended to the back of
# the queue as an unverified candidate, superseding any earlier entry, and is tested only
# when that entry is dequeued. A candidate whose neighbors change several times before it
# is dequeued is thus tested once, while the round-robin order (changed vertices go to the
# back) keeps clipped triangles local.
# `entry_ver[i]` is the `vert_ver[i]` of `i`'s live queue entry (0 = not queued);
# `verified[i]` is the `vert_ver[i]` at which `i` was last confirmed an ear.
function _ear_clip_run!(st::_EarClipState)
    n = length(st.poly)
    ear_queue = Int[]
    sizehint!(ear_queue, 3n)
    entry_ver = zeros(Int, n)
    vert_ver = ones(Int, n)
    verified = zeros(Int, n)
    for i in 1:n
        if _is_ear_vertex(st, i)
            push!(ear_queue, i)
            entry_ver[i] = 1
            verified[i] = 1
        end
    end

    diagonals = Tuple{Int,Int}[]
    sizehint!(diagonals, n - 3)
    triangles = NTuple{3,Int}[]
    sizehint!(triangles, n - 2)

    active_count = n
    active_hint = 1
    q_pos = 1
    pushed_ver = Int[]  # vert_ver at push time, parallel to ear_queue
    sizehint!(pushed_ver, 3n)
    append!(pushed_ver, ones(Int, length(ear_queue)))

    while active_count > 3
        if q_pos > length(ear_queue)
            fb = _find_fallback_ear(st, active_hint)
            fb == 0 && _throw_no_ear(active_count)
            push!(ear_queue, fb)
            push!(pushed_ver, vert_ver[fb])
            entry_ver[fb] = vert_ver[fb]
        end

        ear = ear_queue[q_pos]
        ver = pushed_ver[q_pos]
        q_pos += 1
        (entry_ver[ear] == ver && st.is_active[ear]) || continue  # superseded entry
        entry_ver[ear] = 0
        (verified[ear] == vert_ver[ear] || _is_ear_vertex(st, ear)) || continue

        pr = st.prev_v[ear]
        nx = st.next_v[ear]

        push!(triangles, (pr, ear, nx))
        d = abs(pr - nx)
        if d != 1 && d != n - 1
            push!(diagonals, pr < nx ? (pr, nx) : (nx, pr))
        end

        _remove_reflex!(st.rindex, ear)
        st.next_v[pr] = nx
        st.prev_v[nx] = pr
        st.is_active[ear] = false
        active_count -= 1
        active_hint = nx

        _update_reflex!(st, pr)
        _update_reflex!(st, nx)
        for nb in (pr, nx)
            vert_ver[nb] += 1
            push!(ear_queue, nb)
            push!(pushed_ver, vert_ver[nb])
            entry_ver[nb] = vert_ver[nb]
        end
    end

    if active_count == 3
        v1 = st.is_active[active_hint] ? active_hint : findfirst(st.is_active)
        if v1 !== nothing
            v2 = st.next_v[v1]
            v3 = st.next_v[v2]
            push!(triangles, (v1, v2, v3))
        end
    end

    return PolygonTriangulation(diagonals, triangles)
end

# Shared ear-clipping driver: classify reflex vertices (collinear counts as reflex, so it
# blocks conservatively), build the reflex index with `make_index(is_reflex, reflex_cnt)`,
# and run the ear queue.
function _ear_clip_with_index(make_index, poly, n::Int, kernel::AbstractKernel)
    if n == 3
        return PolygonTriangulation(Tuple{Int,Int}[], NTuple{3,Int}[(1, 2, 3)])
    end
    prev_v, next_v, is_active = _init_ear_clip_ring(n)
    is_reflex = falses(n)
    reflex_cnt = 0
    @inbounds for i in 1:n
        if orient_sign(kernel, poly[prev_v[i]], poly[i], poly[next_v[i]]) <= 0
            is_reflex[i] = true
            reflex_cnt += 1
        end
    end
    st = _EarClipState(kernel, poly, prev_v, next_v, is_active, make_index(is_reflex, reflex_cnt))
    return _ear_clip_run!(st)
end

_ear_clipping_from_valid(poly, n::Int; kernel::AbstractKernel = FloatingKernel()) =
    _ear_clip_with_index((is_reflex, _) -> _build_plain_reflex_index(is_reflex), poly, n, kernel)

"""
    triangulate_ear_cut(vertices; kernel=FloatingKernel()) -> PolygonTriangulation

Triangulate a simple counterclockwise polygon with `n >= 3` vertices using reflex-vertex
ear clipping in `O(n r)` time, where `r` is the number of reflex vertices (`O(n)` fast path
on strictly x-monotone or y-monotone polygons).

# Preconditions
- `vertices` must form a simple counterclockwise polygon with `n >= 3` vertices.
"""
function triangulate_ear_cut(
    vertices;
    kernel::AbstractKernel = FloatingKernel(),
)
    n = _validate_ccw_polygon(kernel, vertices)

    if n == 3
        return PolygonTriangulation(Tuple{Int,Int}[], NTuple{3,Int}[(1, 2, 3)])
    end

    if is_x_monotone(vertices)
        return triangulate_monotone(vertices; kernel = kernel)
    elseif is_y_monotone(vertices)
        rotated = [point(p.y, -p.x) for p in vertices]
        res = triangulate_monotone(rotated; kernel = kernel)
        return PolygonTriangulation(res.diagonals, res.triangles)
    else
        return _ear_clipping_from_valid(vertices, n; kernel = kernel)
    end
end

"""
    triangulate_polygon(vertices; kernel=FloatingKernel(), algorithm=:sweep, rng=Random.default_rng()) -> PolygonTriangulation

Triangulate a simple counterclockwise polygon with `n >= 3` vertices.
Returns a `PolygonTriangulation` with:
- `result.diagonals`: `n - 3` distinct internal diagonals as index pairs `(i, j)` with `i < j`.
- `result.triangles`: `n - 2` counterclockwise triangles as index triples `(i, j, m)`.

Supported `algorithm` options:
- `:sweep` (default): Plane-sweep monotone partitioning and triangulation running in worst-case `O(n log n)` time.
- `:seidel`: Raimund Seidel's randomized incremental algorithm running in expected `O(n log* n)` time.
- `:ear_cut_z_order`: Ear clipping indexing reflex vertices along a Morton Z-order curve (`O(n r)` worst case; queries prune only by the ear's bounding box, so long thin ears and far outlier vertices degrade it toward quadratic in practice).
- `:ear_cut_kd_tree`: Ear clipping indexing active reflex vertices in a `Kd_tree` (`BBT.FlatTree2D`) spatial tree (`O(n r)` worst case; measured roughly `n^1.2` on random star polygons).
- `:ear_cut_u_kd_tree`: Ear clipping indexing active reflex vertices in an unshrunk-cell `Kd_tree` (`BBT.FlatUTree2D`) spatial tree with uncapped upward cell-containment climbing (`O(n r)` worst case; subquadratic in practice).
- `:ear_cut` (also accepted as `:ear_clipping`, its pre-0.5.15 name): Ear clipping with reflex-vertex acceleration running in `O(n r)` time (`O(n)` on x-monotone or y-monotone polygons).

# Preconditions
- `vertices` must form a simple counterclockwise polygon with `n >= 3` vertices.
"""
function triangulate_polygon(
    vertices;
    kernel::AbstractKernel = FloatingKernel(),
    algorithm::Symbol = :sweep,
    rng::AbstractRNG = Random.default_rng(),
)
    if algorithm == :sweep
        return triangulate_sweep(vertices; kernel = kernel)
    elseif algorithm == :seidel
        return triangulate_seidel(vertices; kernel = kernel, rng = rng)
    elseif algorithm == :ear_cut_z_order
        return triangulate_ear_cut_z_order(vertices; kernel = kernel)
    elseif algorithm == :ear_cut_kd_tree
        return triangulate_ear_cut_kd_tree(vertices; kernel = kernel)
    elseif algorithm == :ear_cut_u_kd_tree
        return triangulate_ear_cut_u_kd_tree(vertices; kernel = kernel)
    elseif algorithm == :ear_cut || algorithm == :ear_clipping
        return triangulate_ear_cut(vertices; kernel = kernel)
    else
        throw(ArgumentError("unknown triangulation algorithm: $algorithm (supported: :sweep, :seidel, :ear_cut_z_order, :ear_cut_kd_tree, :ear_cut_u_kd_tree, :ear_cut)"))
    end
end

"""
    polygon_diagonals(vertices; kernel=FloatingKernel(), algorithm=:sweep, rng=Random.default_rng())

Return the internal diagonals of a counterclockwise simple polygon as index pairs `(i, j)` with `i < j`.
"""
@inline function polygon_diagonals(
    vertices;
    kernel::AbstractKernel = FloatingKernel(),
    algorithm::Symbol = :sweep,
    rng::AbstractRNG = Random.default_rng(),
)
    return triangulate_polygon(vertices; kernel = kernel, algorithm = algorithm, rng = rng).diagonals
end

export PolygonTriangulation, MonotoneTriangulation
export triangulate_polygon, triangulate_monotone, triangulate_ear_cut
export is_x_monotone, is_y_monotone
export monotone_chains, monotone_chain_indices, merge_chains
export stack_diagonals, polygon_diagonals
