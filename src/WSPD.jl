module WSPD

using ..BasicCompGeometry
using ..BasicCompGeometry.BBT
using ..BasicCompGeometry.VirtArray
using DataStructures

"""
    WSPair{D, T, S, V}

A pair of nodes in the BBT that may be well-separated.
`S` is the type for distance/diameter (usually Float64).
"""
struct WSPair{D,T,S,V}
    left::BBT.Node{D,T,S,V}
    right::BBT.Node{D,T,S,V}

    dist::S             # Distance between the boxes
    max_sides_diam::S   # Maximum diameter of the two boxes
    diam_ub::S          # Upper bound on distance between any two points in the boxes
end

function WSPair_init(left::BBT.Node{D,T,S,V}, right::BBT.Node{D,T,S,V}) where {D,T,S,V}
    d = dist(left.bb, right.bb)
    max_side = max(diam(left.bb), diam(right.bb))
    ub = max_dist(left.bb, right.bb)
    return WSPair{D,T,S,V}(left, right, d, max_side, ub)
end

const MAX_SEP = 1e18 # Representation of very large separation

function separation(pair::WSPair)
    if pair.dist <= 0.0
        return MAX_SEP
    end
    return pair.max_sides_diam / pair.dist
end

# Ordering for the priority queue (max-heap by diam_ub)
struct WSPDOrder <: Base.Order.Ordering end
Base.Order.lt(::WSPDOrder, a::WSPair, b::WSPair) = a.diam_ub > b.diam_ub

const PairInt = Tuple{Int,Int}

"""
    PD{D, T, S, V}

WSPD construction state.
"""
mutable struct PD{D,T,S,V}
    pnts::AbsPntSeq{D,T}
    finals::Vector{WSPair{D,T,S,V}}
    heap::BinaryHeap{WSPair{D,T,S,V},WSPDOrder}
    tree::BBT.Tree{D,T,S,V}
    sep::S  # Desired quality of separation (s-well-separated)
    hash_pairs::Dict{PairInt,Bool}
end

@inline function get_id(u::BBT.Node, v::BBT.Node)::PairInt
    return (min(u.id, v.id), max(u.id, v.id))
end

function push_pair!(W::PD{D,T,S,V}, u::BBT.Node{D,T,S,V}, v::BBT.Node{D,T,S,V}) where {D,T,S,V}
    id = get_id(u, v)
    if haskey(W.hash_pairs, id)
        return nothing
    end

    pair = WSPair_init(u, v)
    push!(W.heap, pair)
    W.hash_pairs[id] = true
    return pair
end

"""
    top_refine!(W)

Takes the top of the active pairs and refines it by splitting the larger node.
"""
function top_refine!(W::PD{D,T,S,V}) where {D,T,S,V}
    isempty(W.heap) && return nothing

    top = pop!(W.heap)
    BBT.Tree_refine_node(W.tree, top.left)
    BBT.Tree_refine_node(W.tree, top.right)

    l_diam = diam(top.left.bb)
    r_diam = diam(top.right.bb)

    # If both are single points (diam=0), check if they are the same point
    if l_diam == 0.0 && r_diam == 0.0
        if top.dist > 0.0
            push!(W.finals, top)
        end
        return nothing
    end

    # Refine the node with the larger diameter
    if l_diam >= r_diam
        push_pair!(W, top.left.left, top.right)
        push_pair!(W, top.left.right, top.right)
    else
        push_pair!(W, top.left, top.right.left)
        push_pair!(W, top.left, top.right.right)
    end
    return nothing
end

"""
    top_finalize!(W)

Moves the top of the heap to the finals list.
"""
function top_finalize!(W::PD)
    isempty(W.heap) && return nothing
    push!(W.finals, pop!(W.heap))
    return nothing
end

"""
    init(pnts, sep)

Initialize the WSPD construction state.
"""
function init(pnts::AbsPntSeq{D,T}, sep::S) where {D,T,S}
    tree = BBT.Tree_init(pnts)
    V = typeof(Points(pnts))

    PairT = WSPair{D,T,S,V}
    heap = BinaryHeap{PairT}(WSPDOrder(), PairT[])

    W = PD{D,T,S,V}(pnts, PairT[], heap, tree, sep, Dict{PairInt,Bool}())

    # Start with (root, root) pair
    push_pair!(W, tree.root, tree.root)

    return W
end

"""
    expand!(W)

Construct the whole WSPD.
"""
function expand!(W::PD)
    while !isempty(W.heap)
        top = first(W.heap)
        if separation(top) > W.sep
            top_refine!(W)
        else
            top_finalize!(W)
        end
    end
    return W.finals
end

"""
    expand_bichromatic!(W, index_cut)

Construct WSPD for pairs where one point is in 1:index_cut and the other is in index_cut+1:end.
"""
function expand_bichromatic!(W::PD, index_cut::Int)
    while !isempty(W.heap)
        top = first(W.heap)

        i_min = min_orig_index(W, top)
        i_max = max_orig_index(W, top)

        # If both nodes are entirely on one side of the cut, ignore this pair
        f_boring = (i_max <= index_cut) || (i_min > index_cut)

        if f_boring
            pop!(W.heap)
            continue
        end

        if separation(top) > W.sep
            top_refine!(W)
        else
            top_finalize!(W)
        end
    end
    return W.finals
end

"""
    get_top(W)

Return the top pair from the priority queue without removing it.
"""
function get_top(W::PD)
    return first(W.heap)
end

"""
    top_diam_ub(W)

Return the diameter upper bound of the top pair in the priority queue.
"""
function top_diam_ub(W::PD{D,T,S,V}) where {D,T,S,V}
    if isempty(W.heap)
        return zero(S)
    end
    return first(W.heap).diam_ub
end

function get_reps(W::PD, pair::WSPair)
    return W.tree.PS[first(pair.left.r)], W.tree.PS[first(pair.right.r)]
end

function reps_orig_indexes(W::PD, pair::WSPair)
    return orig_index(W.tree.PS, first(pair.left.r)),
    orig_index(W.tree.PS, first(pair.right.r))
end

function min_orig_index(W::PD, pair::WSPair)
    return min(
        BBT.get_min_orig_index(W.tree, pair.left),
        BBT.get_min_orig_index(W.tree, pair.right),
    )
end

function max_orig_index(W::PD, pair::WSPair)
    return max(
        BBT.get_max_orig_index(W.tree, pair.left),
        BBT.get_max_orig_index(W.tree, pair.right),
    )
end

function get_orig_ranges(W::PD, pair::WSPair)
    l_min = BBT.get_min_orig_index(W.tree, pair.left)
    l_max = BBT.get_max_orig_index(W.tree, pair.left)
    r_min = BBT.get_min_orig_index(W.tree, pair.right)
    r_max = BBT.get_max_orig_index(W.tree, pair.right)
    return l_min:l_max, r_min:r_max
end

export WSPair, PD
export init, expand!, expand_bichromatic!, top_refine!
export get_reps, reps_orig_indexes, get_orig_ranges
export min_orig_index, max_orig_index, get_top, top_diam_ub

end # module WSPD
