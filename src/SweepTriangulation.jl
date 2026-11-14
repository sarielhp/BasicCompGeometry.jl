# Deterministic worst-case O(n log n) simple polygon triangulation via plane-sweep
# monotone partitioning and linear-time monotone polygon triangulation.
# Reference: M. de Berg, O. Cheong, M. van Kreveld, and M. Overmars,
# Computational Geometry: Algorithms and Applications (3rd ed.), Chapter 3;
# M. R. Garey, D. S. Johnson, F. P. Preparata, and R. E. Tarjan (1978).

# Preallocated AVL tree storing active lower-boundary edges (going left-to-right under
# compare_xy, with polygon interior above) ordered from bottom to top at the sweep line.
mutable struct _SweepAVLTree
    left::Vector{Int}
    right::Vector{Int}
    parent::Vector{Int}
    height::Vector{Int}
    key::Vector{Int}    # edge index currently stored in tree node slot
    slot::Vector{Int}   # tree node slot currently holding edge e (0 if inactive)
    root::Int
    alloc::Int
end

function _SweepAVLTree(n::Int)
    return _SweepAVLTree(
        zeros(Int, n),
        zeros(Int, n),
        zeros(Int, n),
        zeros(Int, n),
        zeros(Int, n),
        zeros(Int, n),
        0,
        0,
    )
end

@inline _avl_height(t::_SweepAVLTree, u::Int) = u == 0 ? 0 : t.height[u]

@inline function _avl_update_height!(t::_SweepAVLTree, u::Int)
    t.height[u] = 1 + max(_avl_height(t, t.left[u]), _avl_height(t, t.right[u]))
    return nothing
end

@inline function _avl_replace_child!(t::_SweepAVLTree, p::Int, old_c::Int, new_c::Int)
    if p == 0
        t.root = new_c
    elseif t.left[p] == old_c
        t.left[p] = new_c
    else
        t.right[p] = new_c
    end
    new_c != 0 && (t.parent[new_c] = p)
    return nothing
end

function _avl_rotate_left!(t::_SweepAVLTree, x::Int)::Int
    y = t.right[x]
    p = t.parent[x]
    b = t.left[y]

    t.right[x] = b
    b != 0 && (t.parent[b] = x)

    t.left[y] = x
    t.parent[x] = y

    _avl_replace_child!(t, p, x, y)
    _avl_update_height!(t, x)
    _avl_update_height!(t, y)
    return y
end

function _avl_rotate_right!(t::_SweepAVLTree, y::Int)::Int
    x = t.left[y]
    p = t.parent[y]
    b = t.right[x]

    t.left[y] = b
    b != 0 && (t.parent[b] = y)

    t.right[x] = y
    t.parent[y] = x

    _avl_replace_child!(t, p, y, x)
    _avl_update_height!(t, y)
    _avl_update_height!(t, x)
    return x
end

function _avl_rebalance_up!(t::_SweepAVLTree, start_node::Int)
    u = start_node
    while u != 0
        _avl_update_height!(t, u)
        bf = _avl_height(t, t.left[u]) - _avl_height(t, t.right[u])
        if bf > 1
            lc = t.left[u]
            if _avl_height(t, t.left[lc]) < _avl_height(t, t.right[lc])
                _avl_rotate_left!(t, lc)
            end
            u = _avl_rotate_right!(t, u)
        elseif bf < -1
            rc = t.right[u]
            if _avl_height(t, t.right[rc]) < _avl_height(t, t.left[rc])
                _avl_rotate_right!(t, rc)
            end
            u = _avl_rotate_left!(t, u)
        end
        u = t.parent[u]
    end
    return nothing
end

# Return true iff new lower-boundary edge `e_new = (v -> q_new)` starting at sweep event `v`
# lies above active lower-boundary edge `e_old = (p_old -> q_old)` just to the right of `v`
# under the symbolic shear `(x', y') = (x + εy, y)`.
@inline function _sweep_edge_above(
    kernel::AbstractKernel,
    pts,
    next_v::Vector{Int},
    e_new::Int,
    v::Int,
    e_old::Int,
)::Bool
    q_new = next_v[e_new]
    q_old = next_v[e_old]
    # `e_old != v` always holds: `_avl_insert!` rejects an edge that is already in the tree.
    o = orient_sign(kernel, pts[e_old], pts[q_old], pts[v])
    o != 0 && return o > 0
    return orient_sign(kernel, pts[e_old], pts[q_old], pts[q_new]) > 0
end

function _avl_insert!(
    t::_SweepAVLTree,
    kernel::AbstractKernel,
    pts,
    next_v::Vector{Int},
    e::Int,
    v::Int,
)
    t.slot[e] == 0 ||
        throw(ArgumentError("polygon edges intersect or duplicate lower-boundary edge $e"))
    idx = (t.alloc += 1)
    t.key[idx] = e
    t.slot[e] = idx
    t.height[idx] = 1
    t.left[idx] = 0
    t.right[idx] = 0

    if t.root == 0
        t.root = idx
        t.parent[idx] = 0
        return nothing
    end

    cur = t.root
    par = 0
    go_right = false
    while cur != 0
        par = cur
        go_right = _sweep_edge_above(kernel, pts, next_v, e, v, t.key[cur])
        cur = go_right ? t.right[cur] : t.left[cur]
    end

    t.parent[idx] = par
    if go_right
        t.right[par] = idx
    else
        t.left[par] = idx
    end
    _avl_rebalance_up!(t, par)
    return nothing
end

# Delete edge `e` by handle without geometric comparisons at shared right endpoints.
function _avl_delete!(t::_SweepAVLTree, e::Int)
    z = t.slot[e]
    z == 0 && throw(
        ArgumentError(
            "polygon edges intersect or do not form a simple counterclockwise polygon (edge $e not active)",
        ),
    )
    t.slot[e] = 0

    if t.left[z] != 0 && t.right[z] != 0
        y = t.right[z]
        while t.left[y] != 0
            y = t.left[y]
        end
        ey = t.key[y]
        t.key[z] = ey
        t.slot[ey] = z
        z = y
    end

    child = t.left[z] != 0 ? t.left[z] : t.right[z]
    pz = t.parent[z]
    _avl_replace_child!(t, pz, z, child)
    t.left[z] = 0
    t.right[z] = 0
    t.parent[z] = 0
    t.height[z] = 0
    pz != 0 && _avl_rebalance_up!(t, pz)
    return nothing
end

# Return the highest active lower-boundary edge strictly below vertex `v`, or `0` if none.
function _avl_below(
    t::_SweepAVLTree,
    kernel::AbstractKernel,
    pts,
    next_v::Vector{Int},
    v::Int,
)::Int
    best = 0
    cur = t.root
    while cur != 0
        e = t.key[cur]
        if orient_sign(kernel, pts[e], pts[next_v[e]], pts[v]) > 0
            best = e
            cur = t.right[cur]
        else
            cur = t.left[cur]
        end
    end
    return best
end

@inline function _sweep_emit_diagonal!(
    diags::Vector{Tuple{Int,Int}},
    n::Int,
    u::Int,
    v::Int,
)
    d = abs(u - v)
    if d > 1 && d < n - 1
        push!(diags, u < v ? (u, v) : (v, u))
    end
    return nothing
end

@inline function _sweep_connect_if_merge!(
    diags::Vector{Tuple{Int,Int}},
    n::Int,
    v::Int,
    helper_v::Int,
    is_merge::BitVector,
)
    (helper_v != 0 && is_merge[helper_v]) && _sweep_emit_diagonal!(diags, n, v, helper_v)
    return nothing
end

function _require_edge_below(
    t::_SweepAVLTree,
    kernel::AbstractKernel,
    pts,
    next_v::Vector{Int},
    v::Int,
)::Int
    e_j = _avl_below(t, kernel, pts, next_v, v)
    e_j == 0 && throw(
        ArgumentError(
            "polygon edges intersect or do not form a simple counterclockwise polygon",
        ),
    )
    return e_j
end

# Partition a simple CCW polygon into compare_xy-monotone pieces via O(n log n) left-to-right
# plane sweep under symbolic shear `(x', y') = (x + εy, y)`.
function _sweep_partition_diagonals(kernel::AbstractKernel, pts)
    n = length(pts)
    next_v = [mod1(i + 1, n) for i in 1:n]

    edge_is_lower = falses(n)
    for s in 1:n
        c = compare_xy(kernel, pts[s], pts[next_v[s]])
        c == 0 && throw(ArgumentError("polygon contains duplicate adjacent vertex $s"))
        edge_is_lower[s] = c < 0
    end

    # Use MergeSort for deterministic worst-case O(n log n) event ordering.
    events = collect(1:n)
    sort!(events; alg = MergeSort, lt = (a, b) -> compare_xy(kernel, pts[a], pts[b]) < 0)
    for i in 1:(n-1)
        compare_xy(kernel, pts[events[i]], pts[events[i+1]]) == 0 &&
            throw(ArgumentError("polygon contains duplicate vertex $(events[i])"))
    end

    tree = _SweepAVLTree(n)
    helper = zeros(Int, n)
    is_merge = falses(n)
    diags = Tuple{Int,Int}[]
    sizehint!(diags, n - 3)

    for v in events
        e_in = v == 1 ? n : v - 1
        prev_left = edge_is_lower[e_in]
        next_right = edge_is_lower[v]

        if !prev_left && next_right
            # Both incident edges go strictly right: START (convex) or SPLIT (reflex).
            o = orient_sign(kernel, pts[e_in], pts[v], pts[next_v[v]])
            o == 0 && throw(
                ArgumentError("polygon has degenerate overlapping edges at vertex $v"),
            )
            if o < 0
                e_j = _require_edge_below(tree, kernel, pts, next_v, v)
                _sweep_emit_diagonal!(diags, n, v, helper[e_j])
                helper[e_j] = v
            end
            _avl_insert!(tree, kernel, pts, next_v, v, v)
            helper[v] = v
        elseif prev_left && !next_right
            # Both incident edges come from the left: END (convex) or MERGE (reflex).
            o = orient_sign(kernel, pts[e_in], pts[v], pts[next_v[v]])
            o == 0 && throw(
                ArgumentError("polygon has degenerate overlapping edges at vertex $v"),
            )
            _sweep_connect_if_merge!(diags, n, v, helper[e_in], is_merge)
            _avl_delete!(tree, e_in)
            if o < 0
                is_merge[v] = true
                e_j = _require_edge_below(tree, kernel, pts, next_v, v)
                _sweep_connect_if_merge!(diags, n, v, helper[e_j], is_merge)
                helper[e_j] = v
            end
        elseif prev_left && next_right
            # Regular lower-chain vertex (interior lies above v).
            _sweep_connect_if_merge!(diags, n, v, helper[e_in], is_merge)
            _avl_delete!(tree, e_in)
            _avl_insert!(tree, kernel, pts, next_v, v, v)
            helper[v] = v
        else
            # Regular upper-chain vertex (interior lies below v).
            e_j = _require_edge_below(tree, kernel, pts, next_v, v)
            _sweep_connect_if_merge!(diags, n, v, helper[e_j], is_merge)
            helper[e_j] = v
        end
    end

    tree.root == 0 || throw(
        ArgumentError(
            "polygon edges intersect or do not form a simple counterclockwise polygon",
        ),
    )
    return diags
end

"""
    triangulate_sweep(vertices; kernel=FloatingKernel()) -> PolygonTriangulation

Triangulate a simple counterclockwise polygon in deterministic worst-case `O(n log n)` time
and `O(n)` space using plane-sweep monotone partitioning (Garey–Johnson–Preparata–Tarjan /
de Berg et al.) followed by linear-time monotone polygon triangulation.

The sweep processes vertices in strict `compare_xy` order (symbolic shear `(x + εy, y)`),
maintains active lower-boundary edges in an AVL tree ordered bottom-to-top with merge-vertex
helpers, partitions the polygon into `compare_xy`-monotone pieces in `O(n log n)` time, and
triangulates each piece in `O(k)` time.

# Preconditions
- `vertices` must form a simple counterclockwise polygon with `n >= 3` vertices.
"""
function triangulate_sweep(vertices; kernel::AbstractKernel = FloatingKernel())
    n = _validate_ccw_polygon(kernel, vertices; check_duplicates = false)
    if n == 3
        return PolygonTriangulation(Tuple{Int,Int}[], NTuple{3,Int}[(1, 2, 3)])
    end

    partition_diags = _sweep_partition_diagonals(kernel, vertices)
    return _triangulate_from_partition(kernel, vertices, partition_diags)
end

export triangulate_sweep
