module BBT

using ..BasicCompGeometry
using ..BasicCompGeometry.VirtArray
using LinearAlgebra

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
    )
end

function Tree_init_inner(p::AbsPolygon{D,T}) where {D,T}
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

    cutoff = middle(node.bb, dim_max)
    r_l, r_r = pnt_varray_partition!(tree.PS, node.r, dim_max, cutoff)

    # If partition failed to split, mark as leaf
    if length(r_l) == 0 || length(r_r) == 0
        node.f_leaf = true
        return node
    end

    node.left = node_init(tree, r_l)
    node.right = node_init(tree, r_r)

    return node
end

function node_expand(v::Node{D,T,S,V}, tree::Tree{D,T,S,V}) where {D,T,S,V}
    v.f_leaf && return
    node_split(v, tree)
end

function Tree_init(PS::AbsPolygon{D,T}) where {D,T}
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

# Visualization code commented out as Cairo and Colors are not dependencies
#=
using Cairo, Colors
function   bbox_draw( context, bb, color )
    ...
end
=#

export Node, Tree
export Tree_draw, Tree_init, Tree_fully_expand
export Tree_refine_node, depth
export get_min_orig_index, get_max_orig_index

end # module BBT
