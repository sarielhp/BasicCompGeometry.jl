module BBT

using ..BasicCompGeometry
using ..BasicCompGeometry.VirtArray
using LinearAlgebra
using DataStructures

"""
    Node{D, T, S, V}

A node in the BoundingBoxTree.
"""
mutable struct Node{D,T,S,V}
    bb::BBox{D,T}
    r::UnitRange{Int}
    left::Union{Nothing,Node{D,T,S,V}}
    right::Union{Nothing,Node{D,T,S,V}}
    f_leaf::Bool
    diam::S
    id::Int
    min_index::Union{Nothing,Int}
    max_index::Union{Nothing,Int}
    split_dim::Union{Nothing,Int}
    split_val::Union{Nothing,T}
end

"""
    Tree{D, T, S, V}

BoundingBoxTree structure for efficient spatial queries on a set of points.
"""
mutable struct Tree{D,T,S,V}
    PS::VArray{Point{D,T},V}
    root::Union{Nothing,Node{D,T,S,V}}
    id_counter::Int
end

function get_min_max_orig_index(
    tree::Tree{D,T,S,V},
    node::Node{D,T,S,V},
    min_max_func,
    field_name::Symbol,
)::Int where {D,T,S,V}
    val = getfield(node, field_name)
    if !isnothing(val)
        return val
    end

    if node.f_leaf
        idx = orig_index(tree.PS, first(node.r))
        setfield!(node, field_name, idx)
        return idx
    end

    if !isnothing(node.left) && !isnothing(node.right)
        v_l = get_min_max_orig_index(tree, node.left, min_max_func, field_name)
        v_r = get_min_max_orig_index(tree, node.right, min_max_func, field_name)
        val = min_max_func((v_l, v_r))
    else
        val = min_max_func(orig_index(tree.PS, i) for i in node.r)
    end

    setfield!(node, field_name, val)
    return val
end

function get_min_orig_index(tree::Tree{D,T,S,V}, node::Node{D,T,S,V})::Int where {D,T,S,V}
    return get_min_max_orig_index(tree, node, minimum, :min_index)
end

function get_max_orig_index(tree::Tree{D,T,S,V}, node::Node{D,T,S,V})::Int where {D,T,S,V}
    return get_min_max_orig_index(tree, node, maximum, :max_index)
end

function node_init(tree::Tree{D,T,S,V}, range::UnitRange{Int}) where {D,T,S,V}
    bb = BBox{D,T}()
    tree.id_counter += 1
    bound!(bb, tree.PS[range])
    d = diam(bb)
    return Node{D,T,typeof(d),V}(
        bb,
        range,
        nothing,
        nothing,
        false,
        d,
        tree.id_counter,
        nothing,
        nothing,
        nothing,
        nothing,
    )
end

function Tree_init_inner(p::AbsPntSeq{D,T}) where {D,T}
    pnts = Points(p)
    varr = VArray(pnts)
    V = typeof(pnts)
    # Determine type of diam from an empty bbox or a bound one
    d_temp = diam(bound!(BBox{D,T}(), pnts))
    S = typeof(d_temp)
    tree = Tree{D,T,S,V}(varr, nothing, 0)
    tree.root = node_init(tree, 1:length(p))
    return tree
end

function is_identical(arr::VArray{T,V}, r::UnitRange{Int})::Bool where {T,V}
    length(r) <= 1 && return true
    s = first(r)
    for i = (s+1):last(r)
        if arr[s] != arr[i]
            return false
        end
    end
    return true
end

function pnt_varray_partition!(
    P::VArray{Point{D,T},V},
    r::UnitRange{Int},
    dim,
    cutoff,
) where {D,T,V}
    @assert 1 <= first(r) && last(r) <= length(P)
    top = last(r)
    i = first(r)
    @inbounds while i <= top
        if P[i][dim] > cutoff
            swap!(P, i, top)
            top -= 1
        else
            i += 1
        end
    end
    return first(r):top, (top+1):last(r)
end

function node_split(node::Node{D,T,S,V}, tree::Tree{D,T,S,V}) where {D,T,S,V}
    if !isnothing(node.left) && !isnothing(node.right)
        return node
    end

    if node.f_leaf
        return node
    end

    if is_identical(tree.PS, node.r)
        node.f_leaf = true
        return node
    end

    # Find dimension with maximum width
    w_max = -1.0
    dim_max = 1
    for i = 1:D
        w = width(node.bb, i)
        if w > w_max
            w_max = w
            dim_max = i
        end
    end

    m_min = node.bb.mini[dim_max]
    m_max = node.bb.maxi[dim_max]
    cutoff = (T <: Integer) ? T(fld(m_min + m_max, 2)) : (m_min + m_max) / 2

    r_l, r_r = pnt_varray_partition!(tree.PS, node.r, dim_max, cutoff)

    # If partition failed to split, mark as leaf
    if length(r_l) == 0 || length(r_r) == 0
        node.f_leaf = true
        return node
    end

    node.split_dim = dim_max
    node.split_val = cutoff
    node.left = node_init(tree, r_l)
    node.right = node_init(tree, r_r)

    return node
end

function node_expand(v::Node{D,T,S,V}, tree::Tree{D,T,S,V}) where {D,T,S,V}
    v.f_leaf && return
    node_split(v, tree)
end

function Tree_init(PS::AbsPntSeq{D,T}) where {D,T}
    tree = Tree_init_inner(PS)
    node_expand(tree.root, tree)
    return tree
end

function fully_expand(node::Node{D,T,S,V}, tree::Tree{D,T,S,V}) where {D,T,S,V}
    node.f_leaf && return

    if isnothing(node.left) || isnothing(node.right)
        node_expand(node, tree)
    end

    if !isnothing(node.left)
        fully_expand(node.left, tree)
    end
    if !isnothing(node.right)
        fully_expand(node.right, tree)
    end
end

function Tree_fully_expand(tree::Tree{D,T,S,V}) where {D,T,S,V}
    fully_expand(tree.root, tree)
end

function depth(node::Union{Nothing, Node{D,T,S,V}})::Int where {D,T,S,V}
    if isnothing(node)
        return 0
    end
    if node.f_leaf
        return 1
    end
    return 1 + max(depth(node.left), depth(node.right))
end

function Tree_refine_node(tree, node)
    node_split(node, tree)
end

"""
    exact_naive_scan(tree::Tree{D, T, S, V}, q::Point{D, T})

Compute the exact nearest neighbor of query point `q` in the tree's points
using a brute-force linear scan.
Returns `(min_dist, closest_point, original_index)`.
"""
function exact_naive_scan(tree::Tree{D,T,S,V}, q::Point{D}) where {D,T,S,V}
    min_dist = Inf
    best_p = nothing
    best_orig_idx = -1
    for i in 1:length(tree.PS)
        p = tree.PS[i]
        d = dist(p, q)
        if d < min_dist
            min_dist = d
            best_p = p
            best_orig_idx = orig_index(tree.PS, i)
        end
    end
    return min_dist, best_p, best_orig_idx
end

"""
    exact_naive_scan(P::AbsPntSeq{D, T}, q::Point{D, T})

Compute the exact nearest neighbor of query point `q` in point sequence `P`
using a brute-force linear scan.
Returns `(min_dist, closest_point, index)`.
"""
function exact_naive_scan(P::AbsPntSeq{D,T}, q::Point{D}) where {D,T}
    min_dist = Inf
    best_p = nothing
    best_idx = -1
    for (i, p) in enumerate(P)
        d = dist(p, q)
        if d < min_dist
            min_dist = d
            best_p = p
            best_idx = i
        end
    end
    return min_dist, best_p, best_idx
end

"""
    approx_nn(tree::Tree{D, T, S, V}, q::Point{D, T}, c::Real)

Compute an approximate nearest neighbor to query point `q` in the tree.
The parameter `c >= 1.0` controls the approximation quality.
Returns `(best_dist, closest_point, original_index)`.
"""
function approx_nn(tree::Tree{D,T,S,V}, q::Point{D}, c::Real) where {D,T,S,V}
    @assert(c >= 1.0)
    if isnothing(tree.root)
        error("Tree is empty")
    end

    # 1. Arbitrary point distance
    idx_first = first(tree.root.r)
    best_pnt = tree.PS[idx_first]
    best_dist = dist(q, best_pnt)
    best_orig_index = orig_index(tree.PS, idx_first)

    return approx_nn_with_start(tree, q, c, best_dist, best_pnt, best_orig_index)
end

"""
    hybrid_nn(tree::Tree{D, T, S, V}, q::Point{D, T}, c::Real)

A hybrid nearest neighbor search that first calls `silly_nn` to get an initial
candidate, and then uses that to start `approx_nn`.
"""
function hybrid_nn(tree::Tree{D,T,S,V}, q::Point{D}, c::Real) where {D,T,S,V}
    # 1. Start with silly_nn
    silly_d, silly_p, silly_idx = silly_nn(tree, q)

    # 2. Use it as starting point for approx_nn
    return approx_nn_with_start(tree, q, c, silly_d, silly_p, silly_idx)
end

"""
    approx_nn_with_start(tree, q, c, best_dist, best_pnt, best_orig_index)

Core implementation of approximate nearest neighbor search with a given initial candidate.
"""
function approx_nn_with_start(
    tree::Tree{D,T,S,V},
    q::Point{D},
    c::Real,
    best_dist,
    best_pnt,
    best_orig_index,
) where {D,T,S,V}
    best_dist_over_c = best_dist / c

    # 2. Distance to root bounding box
    d_root = dist(q, tree.root.bb)

    # Priority Queue for nodes to visit (min-priority on distance)
    # Use Float64 for distances to handle cases where tree is T=Int but query is Float
    pq = PriorityQueue{Node{D,T,S,V},Float64}()

    if d_root * c < best_dist
        push!(pq, tree.root => Float64(d_root))
    end

    while !isempty(pq)
        node, d_node_bb = first(pq)
        popfirst!(pq)

        # If it is not going to improve anything, no point looking...
        if d_node_bb >= best_dist
            continue
        end

        # "Whenever a node of the BBT is poped the distance of its point to the query point is computed"
        p_node = tree.PS[first(node.r)]
        d_p = dist(q, p_node)
        if d_p < best_dist
            best_dist = d_p
            best_dist_over_c = best_dist / c
            best_pnt = p_node
            best_orig_index = orig_index(tree.PS, first(node.r))
        end
        if d_node_bb >= best_dist_over_c
            continue
        end

        # If it's not a leaf, push children
        if !node.f_leaf
            # Ensure children are expanded
            if isnothing(node.left) || isnothing(node.right)
                node_expand(node, tree)
            end

            if !isnothing(node.left)
                d_l = dist(q, node.left.bb)
                if d_l * c < best_dist
                    pq[node.left] = d_l
                end
            end
            if !isnothing(node.right)
                d_r = dist(q, node.right.bb)
                if d_r * c < best_dist
                    pq[node.right] = d_r
                end
            end
        end
    end

    return best_dist, best_pnt, best_orig_index
end

"""
    silly_nn(tree::Tree{D, T, S, V}, q::Point{D, T})

A "silly" nearest neighbor search that simply descends the tree to a leaf
based on the split dimension and value at each node.
Returns `(dist, closest_point, original_index)`.
"""
function silly_nn(tree::Tree{D,T,S,V}, q::Point{D}) where {D,T,S,V}
    if isnothing(tree.root)
        error("Tree is empty")
    end

    curr = tree.root
    while !curr.f_leaf
        # Ensure children are expanded if they aren't already
        if isnothing(curr.left) || isnothing(curr.right)
            node_expand(curr, tree)
            if curr.f_leaf
                break
            end
        end

        t = curr.split_dim
        val = curr.split_val

        # According to partitioning logic: points <= cutoff go left
        if q[t] <= val
            curr = curr.left
        else
            curr = curr.right
        end
    end

    # At leaf, return the distance to the point stored there
    idx = first(curr.r)
    pnt = tree.PS[idx]
    d = dist(q, pnt)
    orig_idx = orig_index(tree.PS, idx)

    return d, pnt, orig_idx
end

"""
    Tree_draw(tree, filename)

A visualization utility for BoundingBoxTrees. 
Draws the tree levels into a multi-page PDF.
Requires `Cairo` and `Colors` to be loaded.
"""
function Tree_draw end

export Node, Tree
export Tree_draw, Tree_init, Tree_fully_expand
export Tree_refine_node, depth
export get_min_orig_index, get_max_orig_index, approx_nn, silly_nn, hybrid_nn, exact_naive_scan

end # module BBT
