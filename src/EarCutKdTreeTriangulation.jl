# Ear-clipping polygon triangulation using the flat array-backed unshrunk-cell 2D Kd-tree
# (`Kd_tree.FlatUTree2D`) over active reflex
# vertices, with swap-to-end live-leaf tracking (O(depth) per update, for the per-node
# active counts), parent climbing, and cache-line-packed halfplane pruning.

mutable struct _KdTreeReflexIndex{T,FT<:Kd_tree.FlatUTree2D} <: _AbstractReflexIndex
    is_reflex::Vector{Bool}
    tree::FT
    live_end::Vector{Int32}
    active_count::Vector{Int32}
    overflow::Vector{Int}
    stack::Vector{Int32}
    blocker_hint::Int32
end

function _build_kdtree_reflex_index(vertices, is_reflex::Vector{Bool}, reflex_cnt::Int)
    n = length(vertices)
    T = if isconcretetype(eltype(vertices))
        typeof(@inbounds vertices[1].x)
    else
        promote_type((typeof(p.x) for p in vertices)...)
    end
    px = Vector{T}(undef, reflex_cnt)
    py = Vector{T}(undef, reflex_cnt)
    orig_ids = Vector{Int32}(undef, reflex_cnt)
    k = 0
    @inbounds for i in 1:n
        if is_reflex[i]
            k += 1
            p = vertices[i]
            px[k] = p.x
            py[k] = p.y
            orig_ids[k] = Int32(i)
        end
    end

    tree = Kd_tree._build_flat_tree!(Kd_tree.FlatUNode2D{T}, px, py, orig_ids, n, 16)
    num_nodes = length(tree.nodes)
    live_end = Vector{Int32}(undef, num_nodes)
    active_count = Vector{Int32}(undef, num_nodes)
    nodes = tree.nodes
    @inbounds for id in 1:num_nodes
        nd = nodes[id]
        live_end[id] = nd.r_end
        active_count[id] = nd.r_end - nd.r_start + Int32(1)
    end
    stack = Int32[]
    sizehint!(stack, 64)

    return _KdTreeReflexIndex{T,typeof(tree)}(
        is_reflex,
        tree,
        live_end,
        active_count,
        Int[],
        stack,
        Int32(0),
    )
end

@inline function _swap_tree_slots!(ft::Kd_tree.FlatUTree2D{T}, p1::Int, p2::Int) where {T}
    px = ft.px
    py = ft.py
    perm = ft.perm
    slot_of = ft.slot_of
    @inbounds begin
        v1 = perm[p1]
        v2 = perm[p2]
        px[p1], px[p2] = px[p2], px[p1]
        py[p1], py[p2] = py[p2], py[p1]
        perm[p1] = v2
        perm[p2] = v1
        slot_of[v1] = Int32(p2)
        slot_of[v2] = Int32(p1)
    end
    return nothing
end

function _add_reflex!(idx::_KdTreeReflexIndex{T}, v::Int) where {T}
    @inbounds idx.is_reflex[v] && return nothing
    @inbounds idx.is_reflex[v] = true
    ft = idx.tree
    c = @inbounds Int(ft.leaf_of[v])
    if c != 0
        @inbounds begin
            le_next = Int(idx.live_end[c]) + 1
            pos = Int(ft.slot_of[v])
            if pos != le_next
                _swap_tree_slots!(ft, pos, le_next)
            end
            idx.live_end[c] = Int32(le_next)
            ac = idx.active_count
            curr = c
            while curr != 0
                ac[curr] += Int32(1)
                curr = Kd_tree.parent(ft, curr)
            end
        end
    else
        push!(idx.overflow, v)
    end
    return nothing
end

@inline function _remove_reflex!(idx::_KdTreeReflexIndex{T}, i::Int) where {T}
    @inbounds idx.is_reflex[i] || return nothing
    @inbounds idx.is_reflex[i] = false
    if idx.blocker_hint == Int32(i)
        idx.blocker_hint = Int32(0)
    end
    ft = idx.tree
    c = @inbounds Int(ft.leaf_of[i])
    if c != 0
        @inbounds begin
            le = Int(idx.live_end[c])
            pos = Int(ft.slot_of[i])
            if pos < le
                _swap_tree_slots!(ft, pos, le)
            end
            idx.live_end[c] = Int32(le - 1)
            ac = idx.active_count
            curr = c
            while curr != 0
                ac[curr] -= Int32(1)
                curr = Kd_tree.parent(ft, curr)
            end
        end
    elseif !isempty(idx.overflow)
        ov = idx.overflow
        @inbounds for j in 1:length(ov)
            if ov[j] == i
                ov[j] = ov[end]
                pop!(ov)
                break
            end
        end
    end
    return nothing
end

struct _FlatTriQueryF64
    ax::Float64
    ay::Float64
    dx1::Float64
    dy1::Float64
    adx1::Float64
    ady1::Float64
    bx::Float64
    by::Float64
    dx2::Float64
    dy2::Float64
    adx2::Float64
    ady2::Float64
    cx::Float64
    cy::Float64
    dx3::Float64
    dy3::Float64
    adx3::Float64
    ady3::Float64
    min_tx::Float64
    max_tx::Float64
    min_ty::Float64
    max_ty::Float64
end

@inline function _make_tri_query(
    kernel::FloatingKernel{Float64},
    a::Point{2,Float64},
    b::Point{2,Float64},
    c::Point{2,Float64},
    min_tx::Float64,
    max_tx::Float64,
    min_ty::Float64,
    max_ty::Float64,
)
    dx1 = b.x - a.x
    dy1 = b.y - a.y
    dx2 = c.x - b.x
    dy2 = c.y - b.y
    dx3 = a.x - c.x
    dy3 = a.y - c.y
    return _FlatTriQueryF64(
        a.x,
        a.y,
        dx1,
        dy1,
        abs(dx1),
        abs(dy1),
        b.x,
        b.y,
        dx2,
        dy2,
        abs(dx2),
        abs(dy2),
        c.x,
        c.y,
        dx3,
        dy3,
        abs(dx3),
        abs(dy3),
        min_tx,
        max_tx,
        min_ty,
        max_ty,
    )
end

struct _FlatTriQueryGeneric{K<:AbstractKernel,T}
    kernel::K
    a::Point{2,T}
    b::Point{2,T}
    c::Point{2,T}
    min_tx::T
    max_tx::T
    min_ty::T
    max_ty::T
end

@inline function _make_tri_query(
    kernel::AbstractKernel,
    a::Point{2,T},
    b::Point{2,T},
    c::Point{2,T},
    min_tx::T,
    max_tx::T,
    min_ty::T,
    max_ty::T,
) where {T}
    return _FlatTriQueryGeneric(kernel, a, b, c, min_tx, max_tx, min_ty, max_ty)
end

@inline function _box_disjoint_edge_f64(
    ax::Float64,
    ay::Float64,
    dx::Float64,
    dy::Float64,
    adx::Float64,
    ady::Float64,
    x0::Float64,
    x1::Float64,
    y0::Float64,
    y1::Float64,
)::Bool
    px = dy <= 0.0 ? x1 : x0
    py = dx >= 0.0 ? y1 : y0
    sx = max(abs(x0 - ax), abs(x1 - ax))
    sy = max(abs(y0 - ay), abs(y1 - ay))
    tol = 1e-11 * (adx * sy + ady * sx)
    return dx * (py - ay) - dy * (px - ax) < -tol
end

@inline function _box_hits(q::_FlatTriQueryF64, nd::Kd_tree.FlatUNode2D{Float64})::Bool
    x0 = nd.x0
    x1 = nd.x1
    (x1 < q.min_tx || x0 > q.max_tx) && return false
    y0 = nd.y0
    y1 = nd.y1
    (y1 < q.min_ty || y0 > q.max_ty) && return false
    _box_disjoint_edge_f64(q.ax, q.ay, q.dx1, q.dy1, q.adx1, q.ady1, x0, x1, y0, y1) &&
        return false
    _box_disjoint_edge_f64(q.bx, q.by, q.dx2, q.dy2, q.adx2, q.ady2, x0, x1, y0, y1) &&
        return false
    _box_disjoint_edge_f64(q.cx, q.cy, q.dx3, q.dy3, q.adx3, q.ady3, x0, x1, y0, y1) &&
        return false
    return true
end

@inline function _box_disjoint_edge_generic(
    kernel::AbstractKernel,
    a::Point{2,T},
    b::Point{2,T},
    x0::T,
    x1::T,
    y0::T,
    y1::T,
)::Bool where {T}
    px = (b.y <= a.y) ? x1 : x0
    py = (b.x >= a.x) ? y1 : y0
    return orient_sign(kernel, a, b, Point{2,T}(px, py)) < 0
end

@inline function _box_hits(
    q::_FlatTriQueryGeneric{K,T},
    nd::Kd_tree.FlatUNode2D,
)::Bool where {K,T}
    x0 = T(nd.x0)
    x1 = T(nd.x1)
    (x1 < q.min_tx || x0 > q.max_tx) && return false
    y0 = T(nd.y0)
    y1 = T(nd.y1)
    (y1 < q.min_ty || y0 > q.max_ty) && return false
    _box_disjoint_edge_generic(q.kernel, q.a, q.b, x0, x1, y0, y1) && return false
    _box_disjoint_edge_generic(q.kernel, q.b, q.c, x0, x1, y0, y1) && return false
    _box_disjoint_edge_generic(q.kernel, q.c, q.a, x0, x1, y0, y1) && return false
    return true
end

@inline function _pt_in_tri(q::_FlatTriQueryF64, px::Float64, py::Float64)::Bool
    return (q.dx1 * (py - q.ay) - q.dy1 * (px - q.ax) >= 0.0) &&
           (q.dx2 * (py - q.by) - q.dy2 * (px - q.bx) >= 0.0) &&
           (q.dx3 * (py - q.cy) - q.dy3 * (px - q.cx) >= 0.0)
end

@inline function _pt_in_tri(q::_FlatTriQueryGeneric{K,T}, px, py)::Bool where {K,T}
    return _point_in_ccw_triangle(q.kernel, Point{2,T}(T(px), T(py)), q.a, q.b, q.c)
end

@inline function _in_unshrunk_cell(
    nd::Kd_tree.FlatUNode2D{T},
    min_tx::T,
    max_tx::T,
    min_ty::T,
    max_ty::T,
)::Bool where {T<:AbstractFloat}
    return nd.ux0 < min_tx && max_tx <= nd.ux1 && nd.uy0 < min_ty && max_ty <= nd.uy1
end

@inline function _in_unshrunk_cell(
    nd::Kd_tree.FlatUNode2D{T},
    min_tx::T,
    max_tx::T,
    min_ty::T,
    max_ty::T,
)::Bool where {T}
    bm = nd.bound_mask
    return ((bm & 0x01) == 0 || nd.ux0 < min_tx) &&
           ((bm & 0x02) == 0 || max_tx <= nd.ux1) &&
           ((bm & 0x04) == 0 || nd.uy0 < min_ty) &&
           ((bm & 0x08) == 0 || max_ty <= nd.uy1)
end

@inline function _climb_to_subtree(
    ft::Kd_tree.FlatUTree2D{T},
    start::Int,
    min_tx::T,
    max_tx::T,
    min_ty::T,
    max_ty::T,
)::Int where {T}
    start == 0 && return 1
    nodes = ft.nodes
    curr = start
    @inbounds while curr > 1
        _in_unshrunk_cell(nodes[curr], min_tx, max_tx, min_ty, max_ty) && return curr
        curr = Kd_tree.parent(ft, curr)
    end
    return 1
end

# Branch-free (SIMD-friendly) pre-scan of a leaf's live points for the Float64 query:
# returns `false` when no live point lies in the closed triangle other than its endpoints
# `a` (= `pr`) and `c` (= `nx`), which are recognized by their coordinates (input vertices
# are pairwise distinct; the ear tip `b` is convex and thus never live). Only on `true`
# does the caller run the exact scalar scan that identifies the blocker.
@inline function _leaf_may_block(
    q::_FlatTriQueryF64,
    px::Vector{Float64},
    py::Vector{Float64},
    rs::Int,
    le::Int,
)::Bool
    hit = false
    @inbounds @simd for k in rs:le
        x = px[k]
        y = py[k]
        inside =
            (q.dx1 * (y - q.ay) - q.dy1 * (x - q.ax) >= 0.0) &
            (q.dx2 * (y - q.by) - q.dy2 * (x - q.bx) >= 0.0) &
            (q.dx3 * (y - q.cy) - q.dy3 * (x - q.cx) >= 0.0)
        is_end = ((x == q.ax) & (y == q.ay)) | ((x == q.cx) & (y == q.cy))
        hit |= inside & !is_end
    end
    return hit
end

@inline _leaf_may_block(q, px, py, rs::Int, le::Int)::Bool = true

function _dfs_ear_blocked!(
    idx::_KdTreeReflexIndex{T},
    q,
    pr::Int,
    ear::Int,
    nx::Int,
    b_x,
    b_y,
    u::Int,
)::Bool where {T}
    ft = idx.tree
    nodes = ft.nodes
    ac = idx.active_count
    live_end = idx.live_end
    px_arr = ft.px
    py_arr = ft.py
    perm = ft.perm
    stack = idx.stack
    empty!(stack)
    push!(stack, Int32(u))

    min_tx = q.min_tx
    max_tx = q.max_tx
    min_ty = q.min_ty
    max_ty = q.max_ty

    @inbounds while !isempty(stack)
        nid = Int(pop!(stack))
        nd = nodes[nid]
        if nd.is_leaf
            rs = Int(nd.r_start)
            le = Int(live_end[nid])
            _leaf_may_block(q, px_arr, py_arr, rs, le) || continue
            for k in rs:le
                pxk = px_arr[k]
                pyk = py_arr[k]
                (pxk >= min_tx && pxk <= max_tx && pyk >= min_ty && pyk <= max_ty) ||
                    continue
                v = Int(perm[k])
                (v == pr || v == ear || v == nx) && continue
                if _pt_in_tri(q, pxk, pyk)
                    idx.blocker_hint = Int32(v)
                    return true
                end
            end
        else
            l_id = Int(nd.left)
            r_id = Int(nd.right)
            hit_l = l_id != 0 && ac[l_id] > Int32(0) && _box_hits(q, nodes[l_id])
            hit_r = r_id != 0 && ac[r_id] > Int32(0) && _box_hits(q, nodes[r_id])
            if hit_l & hit_r
                go_r = (nd.split_dim == Int8(1) ? b_x : b_y) > nd.split_val
                if go_r
                    push!(stack, Int32(l_id))
                    push!(stack, Int32(r_id))
                else
                    push!(stack, Int32(r_id))
                    push!(stack, Int32(l_id))
                end
            elseif hit_l
                push!(stack, Int32(l_id))
            elseif hit_r
                push!(stack, Int32(r_id))
            end
        end
    end

    return false
end

function _ear_blocked(
    idx::_KdTreeReflexIndex{T},
    poly,
    kernel::AbstractKernel,
    pr::Int,
    ear::Int,
    nx::Int,
    a::Point,
    b::Point,
    c::Point,
)::Bool where {T}
    ft = idx.tree
    nodes = ft.nodes
    ac = idx.active_count
    has_tree = !isempty(nodes) && (@inbounds ac[1] > Int32(0))
    has_ov = !isempty(idx.overflow)
    (has_tree || has_ov) || return false

    min_tx = T(min(a.x, b.x, c.x))
    max_tx = T(max(a.x, b.x, c.x))
    min_ty = T(min(a.y, b.y, c.y))
    max_ty = T(max(a.y, b.y, c.y))

    bh = Int(idx.blocker_hint)
    @inbounds if bh != 0 && bh != pr && bh != ear && bh != nx && idx.is_reflex[bh]
        pt = poly[bh]
        if pt.x >= min_tx && pt.x <= max_tx && pt.y >= min_ty && pt.y <= max_ty
            if _point_in_ccw_triangle(kernel, pt, a, b, c)
                return true
            end
        end
    end

    if has_ov
        @inbounds for v in idx.overflow
            if v != pr && v != ear && v != nx && idx.is_reflex[v]
                pt = poly[v]
                if pt.x >= min_tx && pt.x <= max_tx && pt.y >= min_ty && pt.y <= max_ty
                    if _point_in_ccw_triangle(kernel, pt, a, b, c)
                        idx.blocker_hint = Int32(v)
                        return true
                    end
                end
            end
        end
    end

    has_tree || return false

    start = @inbounds Int(ft.leaf_of[pr])
    if start == 0
        start = @inbounds Int(ft.leaf_of[nx])
    end

    u = _climb_to_subtree(ft, start, min_tx, max_tx, min_ty, max_ty)

    @inbounds ac[u] == Int32(0) && return false

    q = _make_tri_query(
        kernel,
        Point{2,T}(a.x, a.y),
        Point{2,T}(b.x, b.y),
        Point{2,T}(c.x, c.y),
        min_tx,
        max_tx,
        min_ty,
        max_ty,
    )
    @inbounds _box_hits(q, nodes[u]) || return false

    return _dfs_ear_blocked!(idx, q, pr, ear, nx, T(b.x), T(b.y), u)
end

# Ear clipping of a validated polygon with the unshrunk-cell Kd-tree reflex index; with at
# most 16 reflex vertices the plain reflex list is faster and is used instead.
function _ear_clipping_u_kd_tree(vertices, n::Int, kernel::AbstractKernel)
    return _ear_clip_with_index(vertices, n, kernel) do is_reflex, reflex_cnt
        reflex_cnt <= 16 && return _build_plain_reflex_index(is_reflex)
        return _build_kdtree_reflex_index(vertices, Vector{Bool}(is_reflex), reflex_cnt)
    end
end

"""
    triangulate_ear_cut_u_kd_tree(vertices; kernel=FloatingKernel()) -> PolygonTriangulation

Triangulate a simple counterclockwise polygon with `n >= 3` vertices using ear clipping
with an unshrunk-cell flat 2D `Kd_tree` (`BBT.FlatUTree2D`) spatial index over active
reflex vertices, allowing bottom-up parent climbing to terminate on exact unshrunk split
cells without a depth cap.

This is the index [`triangulate_ear_cut`](@ref) uses on large non-monotone polygons; this
function runs it on every input, without the small-polygon and monotone fast paths. For
inputs with at most 16 reflex vertices (including convex polygons), it falls back to the
direct reflex-list ear clipper.

# Preconditions
- `vertices` must form a simple counterclockwise polygon with `n >= 3` vertices.
"""
function triangulate_ear_cut_u_kd_tree(vertices; kernel::AbstractKernel = FloatingKernel())
    n = _validate_ccw_polygon(kernel, vertices)
    return _ear_clipping_u_kd_tree(vertices, n, kernel)
end

"""
    triangulate_ear_cut_kd_tree(vertices; kernel=FloatingKernel()) -> PolygonTriangulation

Deprecated: use [`triangulate_ear_cut_u_kd_tree`](@ref), which this now calls. The
shrinking-box `BBT.FlatTree2D` index it used was never faster than the unshrunk-cell one.
"""
function triangulate_ear_cut_kd_tree(vertices; kernel::AbstractKernel = FloatingKernel())
    Base.depwarn(
        "`triangulate_ear_cut_kd_tree` is deprecated, use `triangulate_ear_cut_u_kd_tree` instead.",
        :triangulate_ear_cut_kd_tree,
    )
    return triangulate_ear_cut_u_kd_tree(vertices; kernel = kernel)
end

"""
    triangulate_ear_cut_z_order(vertices; kernel=FloatingKernel()) -> PolygonTriangulation

Deprecated: use [`triangulate_ear_cut`](@ref), which this now calls. The Morton Z-order
reflex index was slower than the Kd-tree index on every tested input.
"""
function triangulate_ear_cut_z_order(vertices; kernel::AbstractKernel = FloatingKernel())
    Base.depwarn(
        "`triangulate_ear_cut_z_order` is deprecated, use `triangulate_ear_cut` instead.",
        :triangulate_ear_cut_z_order,
    )
    return triangulate_ear_cut(vertices; kernel = kernel)
end

export triangulate_ear_cut_kd_tree,
    triangulate_ear_cut_u_kd_tree, triangulate_ear_cut_z_order
