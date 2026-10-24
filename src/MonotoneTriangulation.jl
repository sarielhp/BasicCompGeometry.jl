"""
    MonotoneTriangulation

Result of triangulating a simple counterclockwise x-monotone polygon.
Fields:
- `diagonals`: `k - 3` distinct internal diagonals as index pairs `(i, j)` with `i < j`.
- `triangles`: `k - 2` counterclockwise triangles as index triples `(i, j, m)`.
"""
struct MonotoneTriangulation
    diagonals::Vector{Tuple{Int,Int}}
    triangles::Vector{NTuple{3,Int}}
end

Base.iterate(res::MonotoneTriangulation) = (res.diagonals, 1)
Base.iterate(res::MonotoneTriangulation, state) = state == 1 ? (res.triangles, 2) : nothing
Base.length(::MonotoneTriangulation) = 2

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

export MonotoneTriangulation
export triangulate_monotone, monotone_chains, monotone_chain_indices, merge_chains, stack_diagonals
