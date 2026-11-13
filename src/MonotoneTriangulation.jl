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

mutable struct _EarClipState{K<:AbstractKernel,P}
    kernel::K
    poly::P
    prev_v::Vector{Int}
    next_v::Vector{Int}
    is_active::BitVector
    is_reflex::BitVector
    r_prev::Vector{Int}
    r_next::Vector{Int}
    r_head::Int
end

@inline function _add_reflex!(st::_EarClipState, i::Int)
    st.is_reflex[i] && return nothing
    st.is_reflex[i] = true
    st.r_prev[i] = 0
    st.r_next[i] = st.r_head
    st.r_head != 0 && (st.r_prev[st.r_head] = i)
    st.r_head = i
    return nothing
end

@inline function _remove_reflex!(st::_EarClipState, i::Int)
    st.is_reflex[i] || return nothing
    st.is_reflex[i] = false
    rp = st.r_prev[i]
    rn = st.r_next[i]
    if rp != 0
        st.r_next[rp] = rn
    else
        st.r_head = rn
    end
    rn != 0 && (st.r_prev[rn] = rp)
    st.r_prev[i] = 0
    st.r_next[i] = 0
    return nothing
end

@inline function _update_reflex!(st::_EarClipState, i::Int)
    if orient_sign(st.kernel, st.poly[st.prev_v[i]], st.poly[i], st.poly[st.next_v[i]]) <= 0
        _add_reflex!(st, i)
    else
        _remove_reflex!(st, i)
    end
    return nothing
end

function _is_ear_vertex(st::_EarClipState, i::Int)
    pr = st.prev_v[i]
    nx = st.next_v[i]
    orient_sign(st.kernel, st.poly[pr], st.poly[i], st.poly[nx]) > 0 || return false

    poly = st.poly
    kernel = st.kernel
    a = poly[pr]
    b = poly[i]
    c = poly[nx]
    curr = st.r_head
    while curr != 0
        if curr != pr && curr != i && curr != nx
            _point_in_ccw_triangle(kernel, poly[curr], a, b, c) && return false
        end
        curr = st.r_next[curr]
    end
    return true
end

function _ear_clipping_triangulate(
    poly;
    kernel::AbstractKernel = FloatingKernel(),
)
    n = length(poly)
    n >= 3 || throw(ArgumentError("polygon must have at least 3 vertices; got n = $n"))

    area2 = twice_area(poly; kernel = kernel)
    area2 > 0 || throw(ArgumentError("polygon must be oriented counterclockwise with positive area; got twice_area = $area2"))

    if n == 3
        return PolygonTriangulation(Tuple{Int,Int}[], NTuple{3,Int}[(1, 2, 3)])
    end

    st = _EarClipState(
        kernel,
        poly,
        [mod1(i - 1, n) for i in 1:n],
        [mod1(i + 1, n) for i in 1:n],
        trues(n),
        falses(n),
        zeros(Int, n),
        zeros(Int, n),
        0,
    )

    for i in 1:n
        _update_reflex!(st, i)
    end

    ear_queue = Int[]
    sizehint!(ear_queue, 3n)
    in_queue = falses(n)
    for i in 1:n
        if !st.is_reflex[i] && _is_ear_vertex(st, i)
            push!(ear_queue, i)
            in_queue[i] = true
        end
    end

    diagonals = Tuple{Int,Int}[]
    sizehint!(diagonals, n - 3)
    triangles = NTuple{3,Int}[]
    sizehint!(triangles, n - 2)

    @inline function emit_diagonal!(u::Int, v::Int)
        d = abs(u - v)
        if d != 1 && d != n - 1
            push!(diagonals, u < v ? (u, v) : (v, u))
        end
    end

    active_count = n
    q_pos = 1

    while active_count > 3
        if q_pos > length(ear_queue)
            curr = findfirst(st.is_active)
            curr === nothing && break
            start_v = curr
            found = false
            while true
                if _is_ear_vertex(st, curr)
                    push!(ear_queue, curr)
                    in_queue[curr] = true
                    found = true
                    break
                end
                curr = st.next_v[curr]
                curr == start_v && break
            end
            found || break
        end

        ear = ear_queue[q_pos]
        q_pos += 1
        if !in_queue[ear] || !st.is_active[ear] || !_is_ear_vertex(st, ear)
            in_queue[ear] = false
            continue
        end
        in_queue[ear] = false

        pr = st.prev_v[ear]
        nx = st.next_v[ear]

        push!(triangles, (pr, ear, nx))
        emit_diagonal!(pr, nx)

        _remove_reflex!(st, ear)
        st.next_v[pr] = nx
        st.prev_v[nx] = pr
        st.is_active[ear] = false
        active_count -= 1

        _update_reflex!(st, pr)
        _update_reflex!(st, nx)

        for nb in (pr, nx)
            if _is_ear_vertex(st, nb)
                if !in_queue[nb]
                    push!(ear_queue, nb)
                    in_queue[nb] = true
                end
            else
                in_queue[nb] = false
            end
        end
    end

    if active_count == 3
        v1 = findfirst(st.is_active)
        if v1 !== nothing
            v2 = st.next_v[v1]
            v3 = st.next_v[v2]
            push!(triangles, (v1, v2, v3))
        end
    end

    return PolygonTriangulation(diagonals, triangles)
end

"""
    triangulate_polygon(vertices; kernel=FloatingKernel(), algorithm=:ear_clipping, rng=Random.default_rng()) -> PolygonTriangulation

Triangulate a simple counterclockwise polygon with `n >= 3` vertices.
Returns a `PolygonTriangulation` with:
- `result.diagonals`: `n - 3` distinct internal diagonals as index pairs `(i, j)` with `i < j`.
- `result.triangles`: `n - 2` counterclockwise triangles as index triples `(i, j, m)`.

Supported `algorithm` options:
- `:ear_clipping` (default): Fast ear clipping with reflex-vertex acceleration running in `O(n r)` time.
- `:sweep`: Plane-sweep monotone partitioning and triangulation running in worst-case `O(n log n)` time.
- `:seidel`: Raimund Seidel's randomized incremental algorithm running in expected `O(n log* n)` time.

For x-monotone or y-monotone polygons, triangulates in linear time `O(n)`.

# Preconditions
- `vertices` must form a simple counterclockwise polygon with `n >= 3` vertices.
"""
function triangulate_polygon(
    vertices;
    kernel::AbstractKernel = FloatingKernel(),
    algorithm::Symbol = :ear_clipping,
    rng::AbstractRNG = Random.default_rng(),
)
    if algorithm == :seidel
        return triangulate_seidel(vertices; kernel = kernel, rng = rng)
    elseif algorithm == :sweep
        return triangulate_sweep(vertices; kernel = kernel)
    elseif algorithm != :ear_clipping
        throw(ArgumentError("unknown triangulation algorithm: $algorithm (supported: :ear_clipping, :sweep, :seidel)"))
    end

    n = length(vertices)
    n >= 3 || throw(ArgumentError("polygon must have at least 3 vertices; got n = $n"))

    area2 = twice_area(vertices; kernel = kernel)
    area2 > 0 || throw(ArgumentError("polygon must be oriented counterclockwise with positive area; got twice_area = $area2"))

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
        return _ear_clipping_triangulate(vertices; kernel = kernel)
    end
end

"""
    polygon_diagonals(vertices; kernel=FloatingKernel(), algorithm=:ear_clipping, rng=Random.default_rng())

Return the internal diagonals of a counterclockwise simple polygon as index pairs `(i, j)` with `i < j`.
"""
@inline function polygon_diagonals(
    vertices;
    kernel::AbstractKernel = FloatingKernel(),
    algorithm::Symbol = :ear_clipping,
    rng::AbstractRNG = Random.default_rng(),
)
    return triangulate_polygon(vertices; kernel = kernel, algorithm = algorithm, rng = rng).diagonals
end

export PolygonTriangulation, MonotoneTriangulation
export triangulate_polygon, triangulate_monotone
export is_x_monotone, is_y_monotone
export monotone_chains, monotone_chain_indices, merge_chains
export stack_diagonals, polygon_diagonals
