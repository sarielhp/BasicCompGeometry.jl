module BBT

using ..BasicCompGeometry
import ..BasicCompGeometry: Points, bound!
import ..BasicCompGeometry: width
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
    parent::Union{Nothing,Node{D,T,S,V}}
    f_leaf::Bool
    diam::S
    id::Int
    min_index::Union{Nothing,Int}
    max_index::Union{Nothing,Int}
    split_dim::Union{Nothing,Int}
    split_val::Union{Nothing,T}

    function Node{D,T,S,V}(
        bb::BBox{D,T},
        r::UnitRange{Int},
        left::Union{Nothing,Node{D,T,S,V}},
        right::Union{Nothing,Node{D,T,S,V}},
        f_leaf::Bool,
        diam::S,
        id::Int,
        min_index::Union{Nothing,Int},
        max_index::Union{Nothing,Int},
        split_dim::Union{Nothing,Int},
        split_val::Union{Nothing,T},
        parent::Union{Nothing,Node{D,T,S,V}} = nothing,
    ) where {D,T,S,V}
        return new{D,T,S,V}(
            bb,
            r,
            left,
            right,
            parent,
            f_leaf,
            diam,
            id,
            min_index,
            max_index,
            split_dim,
            split_val,
        )
    end
end

"""
    Tree{D, T, S, V}

BoundingBoxTree structure for efficient spatial queries on a set of points.
"""
mutable struct Tree{D,T,S,V}
    PS::VArray{Point{D,T},V}
    root::Union{Nothing,Node{D,T,S,V}}
    id_counter::Int
    leaf_size::Int
    compute_diam::Bool

    function Tree{D,T,S,V}(
        PS::VArray{Point{D,T},V},
        root::Union{Nothing,Node{D,T,S,V}},
        id_counter::Int,
        leaf_size::Int = 1,
        compute_diam::Bool = true,
    ) where {D,T,S,V}
        return new{D,T,S,V}(PS, root, id_counter, max(1, leaf_size), compute_diam)
    end
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

@inline function _bound_varray_range!(
    bb::BBox{D,T},
    PS::VArray{Point{D,T},V},
    r::UnitRange{Int},
) where {D,T,V}
    s = first(r)
    e = last(r)
    s > e && return bb
    vec = PS.vec
    loc = PS.loc
    @inbounds begin
        p0 = vec[loc[s]]
        if D == 2
            minx = p0.x
            maxx = p0.x
            miny = p0.y
            maxy = p0.y
            for idx = (s+1):e
                p = vec[loc[idx]]
                px = p.x
                py = p.y
                px < minx && (minx = px)
                px > maxx && (maxx = px)
                py < miny && (miny = py)
                py > maxy && (maxy = py)
            end
            bb.f_init = true
            bb.mini[1] = minx
            bb.mini[2] = miny
            bb.maxi[1] = maxx
            bb.maxi[2] = maxy
        else
            for d = 1:D
                bb.mini[d] = p0[d]
                bb.maxi[d] = p0[d]
            end
            for idx = (s+1):e
                p = vec[loc[idx]]
                for d = 1:D
                    v = p[d]
                    v < bb.mini[d] && (bb.mini[d] = v)
                    v > bb.maxi[d] && (bb.maxi[d] = v)
                end
            end
            bb.f_init = true
        end
    end
    return bb
end

function node_init(
    tree::Tree{D,T,S,V},
    range::UnitRange{Int},
    parent::Union{Nothing,Node{D,T,S,V}} = nothing,
) where {D,T,S,V}
    bb = BBox{D,T}()
    tree.id_counter += 1
    _bound_varray_range!(bb, tree.PS, range)
    d = tree.compute_diam ? S(diam(bb)) : zero(S)
    return Node{D,T,S,V}(
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
        parent,
    )
end

function Tree_init_inner(
    p::AbsPntSeq{D,T};
    leaf_size::Int = 1,
    compute_diam::Bool = true,
) where {D,T}
    pnts = Points(p)
    varr = VArray(pnts)
    V = typeof(pnts)
    S = isempty(pnts) ? Float64 : typeof(dist(pnts[1], pnts[1]))
    tree = Tree{D,T,S,V}(varr, nothing, 0, leaf_size, compute_diam)
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
    dim::Int,
    cutoff,
) where {D,T,V}
    @assert 1 <= first(r) && last(r) <= length(P)
    vec = P.vec
    loc = P.loc
    top = last(r)
    i = first(r)
    @inbounds while i <= top
        if vec[loc[i]][dim] > cutoff
            loc[i], loc[top] = loc[top], loc[i]
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

    node.f_leaf && return node

    if length(node.r) <= tree.leaf_size
        node.f_leaf = true
        return node
    end

    # Find dimension with maximum width; zero width means all points in r are equal.
    w_max = zero(T)
    dim_max = 1
    @inbounds for i = 1:D
        w = node.bb.maxi[i] - node.bb.mini[i]
        if w > w_max
            w_max = w
            dim_max = i
        end
    end
    if !(w_max > zero(T))
        node.f_leaf = true
        return node
    end

    @inbounds m_min = node.bb.mini[dim_max]
    @inbounds m_max = node.bb.maxi[dim_max]
    cutoff = (T <: Integer) ? T(fld(m_min + m_max, 2)) : (m_min + m_max) / 2

    r_l, r_r = pnt_varray_partition!(tree.PS, node.r, dim_max, cutoff)

    # If partition failed to split, mark as leaf
    if length(r_l) == 0 || length(r_r) == 0
        node.f_leaf = true
        return node
    end

    node.split_dim = dim_max
    node.split_val = cutoff
    node.left = node_init(tree, r_l, node)
    node.right = node_init(tree, r_r, node)

    return node
end

function node_expand(v::Node{D,T,S,V}, tree::Tree{D,T,S,V}) where {D,T,S,V}
    v.f_leaf && return
    node_split(v, tree)
end

function Tree_init(
    PS::AbsPntSeq{D,T};
    leaf_size::Int = 1,
    compute_diam::Bool = true,
) where {D,T}
    tree = Tree_init_inner(PS; leaf_size = leaf_size, compute_diam = compute_diam)
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

depth(::Nothing)::Int = 0

function depth(node::Node{D,T,S,V})::Int where {D,T,S,V}
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
    for i = 1:length(tree.PS)
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

# Push `child` onto `pq` if it could still contain a point improving `best_dist` by factor `c`.
function _push_child_if_useful!(pq, child, q, c, best_dist)
    isnothing(child) && return nothing
    d_child = dist(q, child.bb)
    if d_child * c < best_dist
        pq[child] = d_child
    end
    return nothing
end

# Make sure the children of `node` exist.
function _ensure_children!(node, tree)
    if isnothing(node.left) || isnothing(node.right)
        node_expand(node, tree)
    end
    return nothing
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

        if !node.f_leaf
            _ensure_children!(node, tree)
            _push_child_if_useful!(pq, node.left, q, c, best_dist)
            _push_child_if_useful!(pq, node.right, q, c, best_dist)
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

# ==============================================================================
# Flat array-backed 2D Kd-tree (cache-line packed for high-throughput queries)
# ==============================================================================

"""
    FlatNode2D{T}

A compact node in a flat array-backed 2D [`FlatTree2D`](@ref) (`Kd_tree.FlatTree2D`).

Fields are ordered so `FlatNode2D{Float64}` packs cleanly into 64 bytes (one cache
line). `left`, `right`, and `parent` are `Int32` node indices (`0` denotes none).
"""
struct FlatNode2D{T}
    x0::T
    x1::T
    y0::T
    y1::T
    split_val::T
    left::Int32
    right::Int32
    parent::Int32
    r_start::Int32
    r_end::Int32
    split_dim::Int8
    is_leaf::Bool
end

"""
    FlatTree2D{T}

A flat array-backed 2D Bounding Box Tree / Kd-tree storing points in Structure-of-Arrays
(`px`, `py`) layout alongside a permutation `perm` mapping leaf slots back to caller-supplied
integer identifiers, plus reverse lookup arrays (`leaf_of`, `slot_of`).

Parent navigation uses index-based `parent(ft, id)`, [`parent_node`](@ref) (`parent_node(ft, id)`),
and [`is_root`](@ref) (`is_root(ft, id)`), where `0` denotes no parent.
"""
struct FlatTree2D{T}
    nodes::Vector{FlatNode2D{T}}
    px::Vector{T}
    py::Vector{T}
    perm::Vector{Int32}
    leaf_of::Vector{Int32}
    slot_of::Vector{Int32}
    leaf_size::Int
    n_total::Int
end

Base.isempty(ft::FlatTree2D) = isempty(ft.nodes)
Base.length(ft::FlatTree2D) = length(ft.px)

"""
    FlatUNode2D{T}

A node in an unshrunk-cell flat 2D Kd-tree ([`FlatUTree2D`](@ref), `Kd_tree.FlatUTree2D`),
storing both the tight point bounding box `(x0, x1, y0, y1)` and the unshrunk split cell
`(ux0, ux1, uy0, uy1)` alongside `bound_mask::UInt8` (bits `0..3` mark constrained
`ux0, ux1, uy0, uy1` split walls).

`FlatUNode2D{Float64}` occupies 96 bytes.
"""
struct FlatUNode2D{T}
    x0::T
    x1::T
    y0::T
    y1::T
    split_val::T
    ux0::T
    ux1::T
    uy0::T
    uy1::T
    left::Int32
    right::Int32
    parent::Int32
    r_start::Int32
    r_end::Int32
    split_dim::Int8
    bound_mask::UInt8
    is_leaf::Bool
end

"""
    FlatUTree2D{T}

A flat array-backed 2D Bounding Box Tree / Kd-tree using [`FlatUNode2D`](@ref) nodes that
maintain both tight point bounding boxes and unshrunk split-cell rectangles `(ux0, ux1] x (uy0, uy1]`
for faster bottom-up subtree stopping.

Parent navigation uses index-based `parent(ft, id)`, [`parent_node`](@ref) (`parent_node(ft, id)`),
and [`is_root`](@ref) (`is_root(ft, id)`), where `0` denotes no parent.
"""
struct FlatUTree2D{T}
    nodes::Vector{FlatUNode2D{T}}
    px::Vector{T}
    py::Vector{T}
    perm::Vector{Int32}
    leaf_of::Vector{Int32}
    slot_of::Vector{Int32}
    leaf_size::Int
    n_total::Int
end

const _AnyFlatTree2D = Union{FlatTree2D,FlatUTree2D}
Base.isempty(ft::FlatUTree2D) = isempty(ft.nodes)
Base.length(ft::FlatUTree2D) = length(ft.px)

@inline function _flat_midpoint(a::T, b::T) where {T<:AbstractFloat}
    return T(0.5) * a + T(0.5) * b
end
@inline function _flat_midpoint(a::T, b::T) where {T<:Integer}
    both_odd = (mod(a, 2) == 1 && mod(b, 2) == 1) ? 1 : 0
    return T(fld(a, 2) + fld(b, 2) + both_odd)
end
@inline function _flat_midpoint(a::T, b::T) where {T}
    return (a + b) / 2
end

@inline function _sort5(a::T, b::T, c::T, d::T, e::T) where {T}
    if a > b
        a, b = b, a
    end
    if d > e
        d, e = e, d
    end
    if c > e
        c, e = e, c
    end
    if c > d
        c, d = d, c
    end
    if a > d
        a, d = d, a
    end
    if a > c
        a, c = c, a
    end
    if b > e
        b, e = e, b
    end
    if b > d
        b, d = d, b
    end
    if b > c
        b, c = c, b
    end
    return a, b, c, d, e
end

@inline function _pick_sample_split(
    arr::Vector{T},
    s::Int,
    e::Int,
    m_min::T,
    m_max::T,
) where {T}
    span = e - s
    @inbounds begin
        v1 = arr[s]
        v2 = arr[s + (span >> 2)]
        v3 = arr[(s + e) >> 1]
        v4 = arr[s + ((3 * span) >> 2)]
        v5 = arr[e]
    end
    s1, s2, s3, _, _ = _sort5(v1, v2, v3, v4, v5)
    # m_min <= s1 <= s2 <= s3 and m_min < m_max guarantee any candidate < m_max splits non-emptily
    if s3 < m_max
        return s3
    elseif s2 < m_max
        return s2
    elseif s1 < m_max
        return s1
    end
    return m_min
end

@inline function _partition_flat_xy!(
    px::Vector{T},
    py::Vector{T},
    perm::Vector{Int32},
    s::Int,
    e::Int,
    dim::Int,
    cutoff::T,
)::Int where {T}
    arr = dim == 1 ? px : py
    i = s
    j = e
    @inbounds while i <= j
        while i <= j && arr[i] <= cutoff
            i += 1
        end
        while i <= j && arr[j] > cutoff
            j -= 1
        end
        if i < j
            px[i], px[j] = px[j], px[i]
            py[i], py[j] = py[j], py[i]
            perm[i], perm[j] = perm[j], perm[i]
            i += 1
            j -= 1
        end
    end
    return j
end

@inline _root_ubounds(::Type{FlatNode2D{T}}) where {T} = (zero(T), zero(T), zero(T), zero(T), UInt8(0))
@inline function _root_ubounds(::Type{FlatUNode2D{T}}) where {T}
    if T <: AbstractFloat
        return (T(-Inf), T(Inf), T(-Inf), T(Inf), UInt8(0))
    end
    return (zero(T), zero(T), zero(T), zero(T), UInt8(0))
end

@inline function _make_flat_leaf(
    ::Type{FlatNode2D{T}},
    x0::T,
    x1::T,
    y0::T,
    y1::T,
    ux0::T,
    ux1::T,
    uy0::T,
    uy1::T,
    bmask::UInt8,
    parent_id::Int32,
    s::Int,
    e::Int,
) where {T}
    return FlatNode2D{T}(
        x0,
        x1,
        y0,
        y1,
        zero(T),
        Int32(0),
        Int32(0),
        parent_id,
        Int32(s),
        Int32(e),
        Int8(0),
        true,
    )
end

@inline function _make_flat_leaf(
    ::Type{FlatUNode2D{T}},
    x0::T,
    x1::T,
    y0::T,
    y1::T,
    ux0::T,
    ux1::T,
    uy0::T,
    uy1::T,
    bmask::UInt8,
    parent_id::Int32,
    s::Int,
    e::Int,
) where {T}
    return FlatUNode2D{T}(
        x0,
        x1,
        y0,
        y1,
        zero(T),
        ux0,
        ux1,
        uy0,
        uy1,
        Int32(0),
        Int32(0),
        parent_id,
        Int32(s),
        Int32(e),
        Int8(0),
        bmask,
        true,
    )
end

@inline function _make_flat_internal(
    ::Type{FlatNode2D{T}},
    x0::T,
    x1::T,
    y0::T,
    y1::T,
    cutoff::T,
    ux0::T,
    ux1::T,
    uy0::T,
    uy1::T,
    bmask::UInt8,
    left_id::Int,
    right_id::Int,
    parent_id::Int32,
    s::Int,
    e::Int,
    dim::Int,
) where {T}
    return FlatNode2D{T}(
        x0,
        x1,
        y0,
        y1,
        cutoff,
        Int32(left_id),
        Int32(right_id),
        parent_id,
        Int32(s),
        Int32(e),
        Int8(dim),
        false,
    )
end

@inline function _make_flat_internal(
    ::Type{FlatUNode2D{T}},
    x0::T,
    x1::T,
    y0::T,
    y1::T,
    cutoff::T,
    ux0::T,
    ux1::T,
    uy0::T,
    uy1::T,
    bmask::UInt8,
    left_id::Int,
    right_id::Int,
    parent_id::Int32,
    s::Int,
    e::Int,
    dim::Int,
) where {T}
    return FlatUNode2D{T}(
        x0,
        x1,
        y0,
        y1,
        cutoff,
        ux0,
        ux1,
        uy0,
        uy1,
        Int32(left_id),
        Int32(right_id),
        parent_id,
        Int32(s),
        Int32(e),
        Int8(dim),
        bmask,
        false,
    )
end

@inline function _split_flat_ubounds(
    ::Type{FlatNode2D{T}},
    dim::Int,
    cutoff::T,
    ux0::T,
    ux1::T,
    uy0::T,
    uy1::T,
    bmask::UInt8,
) where {T}
    z = zero(T)
    return (z, z, z, z, UInt8(0), z, z, z, z, UInt8(0))
end

@inline function _split_flat_ubounds(
    ::Type{FlatUNode2D{T}},
    dim::Int,
    cutoff::T,
    ux0::T,
    ux1::T,
    uy0::T,
    uy1::T,
    bmask::UInt8,
) where {T}
    if dim == 1
        return (ux0, cutoff, uy0, uy1, bmask | UInt8(0x02), cutoff, ux1, uy0, uy1, bmask | UInt8(0x01))
    end
    return (ux0, ux1, uy0, cutoff, bmask | UInt8(0x08), ux0, ux1, cutoff, uy1, bmask | UInt8(0x04))
end

@inline _wrap_flat_tree(::Type{FlatNode2D{T}}, nodes, px, py, perm, leaf_of, slot_of, leaf_size, n_total) where {T} =
    FlatTree2D{T}(nodes, px, py, perm, leaf_of, slot_of, leaf_size, n_total)

@inline _wrap_flat_tree(::Type{FlatUNode2D{T}}, nodes, px, py, perm, leaf_of, slot_of, leaf_size, n_total) where {T} =
    FlatUTree2D{T}(nodes, px, py, perm, leaf_of, slot_of, leaf_size, n_total)

function _build_flat_node!(
    nodes::Vector{N},
    px::Vector{T},
    py::Vector{T},
    perm::Vector{Int32},
    leaf_of::Vector{Int32},
    slot_of::Vector{Int32},
    id::Int,
    next_free::Int,
    parent_id::Int32,
    s::Int,
    e::Int,
    leaf_size::Int,
    ux0::T,
    ux1::T,
    uy0::T,
    uy1::T,
    bmask::UInt8,
)::Int where {T,N<:Union{FlatNode2D{T},FlatUNode2D{T}}}
    @inbounds begin
        x0 = px[s]
        x1 = x0
        y0 = py[s]
        y1 = y0
        for k in (s + 1):e
            xk = px[k]
            yk = py[k]
            if xk < x0
                x0 = xk
            elseif xk > x1
                x1 = xk
            end
            if yk < y0
                y0 = yk
            elseif yk > y1
                y1 = yk
            end
        end
    end

    n = e - s + 1
    dx = x1 - x0
    dy = y1 - y0
    if n <= leaf_size || (dx <= zero(T) && dy <= zero(T))
        @inbounds nodes[id] =
            _make_flat_leaf(N, x0, x1, y0, y1, ux0, ux1, uy0, uy1, bmask, parent_id, s, e)
        @inbounds for k in s:e
            orig = perm[k]
            leaf_of[orig] = Int32(id)
            slot_of[orig] = Int32(k)
        end
        return next_free
    end

    dim = dx >= dy ? 1 : 2
    m_min = dim == 1 ? x0 : y0
    m_max = dim == 1 ? x1 : y1
    cutoff = _flat_midpoint(m_min, m_max)

    mid = _partition_flat_xy!(px, py, perm, s, e, dim, cutoff)
    n_left = mid - s + 1
    n_right = e - mid
    if n_left == 0 || n_right == 0 || (n_left << 3) < n || (n_right << 3) < n
        arr = dim == 1 ? px : py
        cand = _pick_sample_split(arr, s, e, m_min, m_max)
        if cand != cutoff || n_left == 0 || n_right == 0
            mid = _partition_flat_xy!(px, py, perm, s, e, dim, cand)
            cutoff = cand
        end
    end

    if mid < s || mid >= e
        @inbounds nodes[id] =
            _make_flat_leaf(N, x0, x1, y0, y1, ux0, ux1, uy0, uy1, bmask, parent_id, s, e)
        @inbounds for k in s:e
            orig = perm[k]
            leaf_of[orig] = Int32(id)
            slot_of[orig] = Int32(k)
        end
        return next_free
    end

    left_id = next_free
    right_id = next_free + 1
    next_free_after_children = next_free + 2

    @inbounds nodes[id] = _make_flat_internal(
        N,
        x0,
        x1,
        y0,
        y1,
        cutoff,
        ux0,
        ux1,
        uy0,
        uy1,
        bmask,
        left_id,
        right_id,
        parent_id,
        s,
        e,
        dim,
    )

    lux0, lux1, luy0, luy1, lmask, rux0, rux1, ruy0, ruy1, rmask =
        _split_flat_ubounds(N, dim, cutoff, ux0, ux1, uy0, uy1, bmask)

    next_free = _build_flat_node!(
        nodes,
        px,
        py,
        perm,
        leaf_of,
        slot_of,
        left_id,
        next_free_after_children,
        Int32(id),
        s,
        mid,
        leaf_size,
        lux0,
        lux1,
        luy0,
        luy1,
        lmask,
    )
    return _build_flat_node!(
        nodes,
        px,
        py,
        perm,
        leaf_of,
        slot_of,
        right_id,
        next_free,
        Int32(id),
        mid + 1,
        e,
        leaf_size,
        rux0,
        rux1,
        ruy0,
        ruy1,
        rmask,
    )
end

function _build_flat_tree!(
    ::Type{N},
    px::Vector{T},
    py::Vector{T},
    perm::Vector{Int32},
    n_total::Int,
    leaf_size::Int,
) where {T,N<:Union{FlatNode2D{T},FlatUNode2D{T}}}
    leaf_size >= 1 || throw(ArgumentError("leaf_size must be at least 1, got $leaf_size"))
    m = length(px)
    length(py) == m || throw(DimensionMismatch("px and py must have the same length"))
    length(perm) == m || throw(DimensionMismatch("perm must match px length"))
    leaf_of = zeros(Int32, n_total)
    slot_of = zeros(Int32, n_total)
    if m == 0
        return _wrap_flat_tree(N, N[], px, py, perm, leaf_of, slot_of, leaf_size, n_total)
    end
    nodes = Vector{N}(undef, 2 * m - 1)
    ux0, ux1, uy0, uy1, bmask = _root_ubounds(N)
    used =
        _build_flat_node!(
            nodes,
            px,
            py,
            perm,
            leaf_of,
            slot_of,
            1,
            2,
            Int32(0),
            1,
            m,
            leaf_size,
            ux0,
            ux1,
            uy0,
            uy1,
            bmask,
        ) - 1
    resize!(nodes, used)
    return _wrap_flat_tree(N, nodes, px, py, perm, leaf_of, slot_of, leaf_size, n_total)
end

_build_flat_tree!(
    px::Vector{T},
    py::Vector{T},
    perm::Vector{Int32},
    n_total::Int,
    leaf_size::Int,
) where {T} = _build_flat_tree!(FlatNode2D{T}, px, py, perm, n_total, leaf_size)

function _extract_xy_perm(pts::Union{AbsPntSeq{2,T},AbstractVector{<:Point{2,T}}}) where {T}
    m = length(pts)
    px = Vector{T}(undef, m)
    py = Vector{T}(undef, m)
    perm = Vector{Int32}(undef, m)
    @inbounds for i in 1:m
        p = pts[i]
        px[i] = p.x
        py[i] = p.y
        perm[i] = Int32(i)
    end
    return px, py, perm, m
end

"""
    Tree_init_flat(P::AbsPntSeq{2,T}; leaf_size::Int = 16) -> FlatTree2D{T}
    Tree_init_flat(pts::AbstractVector{<:Point{2,T}}; leaf_size::Int = 16) -> FlatTree2D{T}
    Tree_init_flat(px::AbstractVector{T}, py::AbstractVector{T}, orig_ids::AbstractVector{<:Integer}, n_total::Integer; leaf_size::Int = 16) -> FlatTree2D{T}

Build and fully expand a flat array-backed 2D Kd-tree ([`FlatTree2D`](@ref)) over the
input 2D points using cache-line-packed [`FlatNode2D`](@ref) nodes, midpoint splits with
5-point sample median skew fallback, and leaf stopping threshold `leaf_size`.
"""
function Tree_init_flat(P::AbsPntSeq{2,T}; leaf_size::Int = 16) where {T}
    px, py, perm, m = _extract_xy_perm(P)
    return _build_flat_tree!(FlatNode2D{T}, px, py, perm, m, leaf_size)
end

function Tree_init_flat(pts::AbstractVector{<:Point{2,T}}; leaf_size::Int = 16) where {T}
    px, py, perm, m = _extract_xy_perm(pts)
    return _build_flat_tree!(FlatNode2D{T}, px, py, perm, m, leaf_size)
end

function Tree_init_flat(
    px::AbstractVector{T},
    py::AbstractVector{T},
    orig_ids::AbstractVector{<:Integer},
    n_total::Integer;
    leaf_size::Int = 16,
) where {T}
    px_c = Vector{T}(px)
    py_c = Vector{T}(py)
    perm_c = Int32[Int32(i) for i in orig_ids]
    return _build_flat_tree!(FlatNode2D{T}, px_c, py_c, perm_c, Int(n_total), leaf_size)
end

"""
    Tree_init_flat_u(P::AbsPntSeq{2,T}; leaf_size::Int = 16) -> FlatUTree2D{T}
    Tree_init_flat_u(pts::AbstractVector{<:Point{2,T}}; leaf_size::Int = 16) -> FlatUTree2D{T}
    Tree_init_flat_u(px::AbstractVector{T}, py::AbstractVector{T}, orig_ids::AbstractVector{<:Integer}, n_total::Integer; leaf_size::Int = 16) -> FlatUTree2D{T}

Build and fully expand a flat array-backed unshrunk-cell 2D Kd-tree ([`FlatUTree2D`](@ref))
over the input 2D points using [`FlatUNode2D`](@ref) nodes that maintain both tight point
bounding boxes and unshrunk split-cell rectangles `(ux0, ux1] x (uy0, uy1]`.
"""
function Tree_init_flat_u(P::AbsPntSeq{2,T}; leaf_size::Int = 16) where {T}
    px, py, perm, m = _extract_xy_perm(P)
    return _build_flat_tree!(FlatUNode2D{T}, px, py, perm, m, leaf_size)
end

function Tree_init_flat_u(pts::AbstractVector{<:Point{2,T}}; leaf_size::Int = 16) where {T}
    px, py, perm, m = _extract_xy_perm(pts)
    return _build_flat_tree!(FlatUNode2D{T}, px, py, perm, m, leaf_size)
end

function Tree_init_flat_u(
    px::AbstractVector{T},
    py::AbstractVector{T},
    orig_ids::AbstractVector{<:Integer},
    n_total::Integer;
    leaf_size::Int = 16,
) where {T}
    px_c = Vector{T}(px)
    py_c = Vector{T}(py)
    perm_c = Int32[Int32(i) for i in orig_ids]
    return _build_flat_tree!(FlatUNode2D{T}, px_c, py_c, perm_c, Int(n_total), leaf_size)
end

import Base: parent

"""
    parent(node::Node) -> Union{Nothing,Node}
    parent(ft::Union{FlatTree2D,FlatUTree2D}, id::Integer) -> Int

Return the parent node of `node` (or `nothing` if root), or the parent node index of
`id` in `ft` (or `0` if `id` is the root node).
"""
parent(node::Node) = node.parent
parent(ft::_AnyFlatTree2D, id::Integer)::Int = @inbounds Int(ft.nodes[id].parent)

"""
    parent_node(node::Node) -> Union{Nothing,Node}
    parent_node(ft::Union{FlatTree2D,FlatUTree2D}, id::Integer) -> Int

Return the parent of `node` / node `id`. Equivalent to `parent`.
"""
parent_node(node::Node) = node.parent
parent_node(ft::_AnyFlatTree2D, id::Integer)::Int = @inbounds Int(ft.nodes[id].parent)

"""
    is_root(node::Node) -> Bool
    is_root(ft::Union{FlatTree2D,FlatUTree2D}, id::Integer) -> Bool

Return `true` if `node` (or node `id` in `ft`) has no parent.
"""
is_root(node::Node) = isnothing(node.parent)
is_root(ft::_AnyFlatTree2D, id::Integer)::Bool = @inbounds ft.nodes[id].parent == Int32(0)

export Node, Tree, FlatNode2D, FlatTree2D, FlatUNode2D, FlatUTree2D
export parent_node, is_root
export Tree_draw, Tree_init, Tree_init_flat, Tree_init_flat_u, Tree_fully_expand
export Tree_refine_node, depth
export get_min_orig_index,
    get_max_orig_index, approx_nn, silly_nn, hybrid_nn, exact_naive_scan

end # module BBT
