# Seidel's randomized incremental polygon triangulation in expected O(n log* n) time.
# Reference: R. Seidel, "A simple and fast incremental randomized algorithm for computing
# trapezoidal decompositions and for triangulating polygons", Computational Geometry:
# Theory and Applications 1 (1991), pp. 51-64.

const _SEIDEL_DAG_LEAF = Int8(0)
const _SEIDEL_DAG_X = Int8(1)
const _SEIDEL_DAG_Y = Int8(2)

mutable struct _SeidelDAGNode
    kind::Int8
    key::Int     # vertex index (_SEIDEL_DAG_X), segment index (_SEIDEL_DAG_Y), or trap id (_SEIDEL_DAG_LEAF)
    left::Int    # child for < x (_SEIDEL_DAG_X) or above (_SEIDEL_DAG_Y)
    right::Int   # child for >= x (_SEIDEL_DAG_X) or below (_SEIDEL_DAG_Y)
end

mutable struct _SeidelTrap
    leftp::Int    # left vertex index (0 = -∞ sentinel)
    rightp::Int   # right vertex index (-1 = +∞ sentinel)
    top::Int      # ceiling segment index (0 = +∞ sentinel)
    bot::Int      # floor segment index (-1 = -∞ sentinel)
    l0::Int       # upper-left neighbor trap id (0 = none)
    l1::Int       # lower-left neighbor trap id (0 = none)
    r0::Int       # upper-right neighbor trap id (0 = none)
    r1::Int       # lower-right neighbor trap id (0 = none)
    sink::Int     # index of corresponding node in query DAG
    alive::Bool
end

struct _SeidelMap{K<:AbstractKernel,V}
    kernel::K
    pts::V
    n::Int
    seg_p::Vector{Int}          # left endpoint (under compare_xy) of segment s
    seg_q::Vector{Int}          # right endpoint (under compare_xy) of segment s
    seg_ltr::Vector{Bool}       # true if CCW edge s -> mod1(s+1, n) goes left-to-right
    traps::Vector{_SeidelTrap}
    dag::Vector{_SeidelDAGNode}
    edge_start_node::Vector{Int}
    edge_inserted::BitVector
    vertex_inserted::BitVector
end

@noinline function _seidel_inconsistent(what::AbstractString)
    throw(
        ArgumentError(
            "Seidel triangulation: inconsistent predicate results ($what); the polygon is not simple, " *
            "or is numerically degenerate for the kernel (try kernel = ExactKernel() or algorithm = :sweep)",
        ),
    )
end

# Advance a trapezoid walk to `nxt`, failing loudly on a missing neighbor or a runaway walk.
@inline function _seidel_step(nxt::Int, steps::Int, limit::Int)::Int
    nxt == 0 && _seidel_inconsistent("walk left the trapezoidal map")
    steps > limit && _seidel_inconsistent("walk did not terminate")
    return nxt
end

# Compare vertex indices where 0 = -∞ and -1 = +∞ under the symbolic shear
# (x', y') = (x + εy, y) for infinitesimal ε > 0 (exact compare_xy order).
@inline function _seidel_cmp_pt(kernel::AbstractKernel, pts, a::Int, b::Int)::Int
    a == b && return 0
    (a == 0 || b == -1) && return -1
    (a == -1 || b == 0) && return 1
    return compare_xy(kernel, pts[a], pts[b])
end

@inline function _seidel_replace_right!(traps::Vector{_SeidelTrap}, nb_id::Int, old_t::Int, new_t::Int)
    nb_id == 0 && return nothing
    nb = traps[nb_id]
    nb.r0 == old_t && (nb.r0 = new_t)
    nb.r1 == old_t && (nb.r1 = new_t)
    return nothing
end

@inline function _seidel_replace_left!(traps::Vector{_SeidelTrap}, nb_id::Int, old_t::Int, new_t::Int)
    nb_id == 0 && return nothing
    nb = traps[nb_id]
    nb.l0 == old_t && (nb.l0 = new_t)
    nb.l1 == old_t && (nb.l1 = new_t)
    return nothing
end

function _seidel_init_map(kernel::AbstractKernel, pts)
    n = length(pts)
    seg_p = Vector{Int}(undef, n)
    seg_q = Vector{Int}(undef, n)
    seg_ltr = Vector{Bool}(undef, n)
    for s in 1:n
        nxt = mod1(s + 1, n)
        if compare_xy(kernel, pts[s], pts[nxt]) < 0
            seg_p[s] = s
            seg_q[s] = nxt
            seg_ltr[s] = true
        else
            seg_p[s] = nxt
            seg_q[s] = s
            seg_ltr[s] = false
        end
    end

    dag = _SeidelDAGNode[_SeidelDAGNode(_SEIDEL_DAG_LEAF, 1, 0, 0)]
    traps = _SeidelTrap[_SeidelTrap(0, -1, 0, -1, 0, 0, 0, 0, 1, true)]
    return _SeidelMap(
        kernel,
        pts,
        n,
        seg_p,
        seg_q,
        seg_ltr,
        traps,
        dag,
        ones(Int, n),
        falses(n),
        falses(n),
    )
end

# Locate the leaf trapezoid containing `pts[p]` slightly along direction `p -> q` (p <_{xy} q).
function _seidel_locate_start(m::_SeidelMap, start_node::Int, p::Int, q::Int)::Int
    kernel = m.kernel
    pts = m.pts
    dag = m.dag
    seg_p = m.seg_p
    seg_q = m.seg_q

    node_idx = start_node
    while true
        node = dag[node_idx]
        if node.kind == _SEIDEL_DAG_LEAF
            return node.key
        elseif node.kind == _SEIDEL_DAG_X
            c = compare_xy(kernel, pts[p], pts[node.key])
            node_idx = c < 0 ? node.left : node.right
        else
            seg = node.key
            a = seg_p[seg]
            b = seg_q[seg]
            o = orient_sign(kernel, pts[a], pts[b], pts[p])
            if o == 0
                o = orient_sign(kernel, pts[a], pts[b], pts[q])
            end
            node_idx = o > 0 ? node.left : node.right
        end
    end
end

# Locate any leaf trapezoid whose closure contains vertex `v` (used once per round at v = 1).
function _seidel_locate_vertex_trap(m::_SeidelMap, v::Int)::Int
    kernel = m.kernel
    pts = m.pts
    dag = m.dag
    seg_p = m.seg_p
    seg_q = m.seg_q

    node_idx = 1
    while true
        node = dag[node_idx]
        if node.kind == _SEIDEL_DAG_LEAF
            return node.key
        elseif node.kind == _SEIDEL_DAG_X
            c = compare_xy(kernel, pts[v], pts[node.key])
            node_idx = c < 0 ? node.left : node.right
        else
            seg = node.key
            o = orient_sign(kernel, pts[seg_p[seg]], pts[seg_q[seg]], pts[v])
            node_idx = o >= 0 ? node.left : node.right
        end
    end
end

function _seidel_collect_hit!(m::_SeidelMap, hit::Vector{Int}, t0::Int, p::Int, q::Int)
    empty!(hit)
    push!(hit, t0)
    kernel = m.kernel
    pts = m.pts
    traps = m.traps

    cur = t0
    limit = length(traps)
    while true
        t = traps[cur]
        rp = t.rightp
        _seidel_cmp_pt(kernel, pts, q, rp) <= 0 && break
        if t.r0 != 0 && t.r1 != 0
            o = orient_sign(kernel, pts[p], pts[q], pts[rp])
            cur = o < 0 ? t.r0 : t.r1
        else
            cur = t.r0 != 0 ? t.r0 : t.r1
        end
        cur = _seidel_step(cur, length(hit), limit)
        push!(hit, cur)
    end
    return hit
end

@inline function _seidel_new_trap!(
    traps::Vector{_SeidelTrap},
    dag::Vector{_SeidelDAGNode},
    leftp::Int,
    rightp::Int,
    top::Int,
    bot::Int;
    l0::Int = 0,
    l1::Int = 0,
    r0::Int = 0,
    r1::Int = 0,
)::Int
    t_id = length(traps) + 1
    sink = length(dag) + 1
    push!(dag, _SeidelDAGNode(_SEIDEL_DAG_LEAF, t_id, 0, 0))
    push!(traps, _SeidelTrap(leftp, rightp, top, bot, l0, l1, r0, r1, sink, true))
    return t_id
end

function _seidel_split_vertical!(m::_SeidelMap, t_id::Int, v::Int)::Tuple{Int,Int}
    traps = m.traps
    dag = m.dag
    T = traps[t_id]
    T.alive = false

    t_leftp = T.leftp
    t_rightp = T.rightp
    t_top = T.top
    t_bot = T.bot
    t_l0 = T.l0
    t_l1 = T.l1
    t_r0 = T.r0
    t_r1 = T.r1
    t_sink = T.sink

    tl_id = _seidel_new_trap!(traps, dag, t_leftp, v, t_top, t_bot; l0 = t_l0, l1 = t_l1)
    tr_id = _seidel_new_trap!(traps, dag, v, t_rightp, t_top, t_bot; l0 = tl_id, r0 = t_r0, r1 = t_r1)
    traps[tl_id].r0 = tr_id

    _seidel_replace_right!(traps, t_l0, t_id, tl_id)
    _seidel_replace_right!(traps, t_l1, t_id, tl_id)
    _seidel_replace_left!(traps, t_r0, t_id, tr_id)
    _seidel_replace_left!(traps, t_r1, t_id, tr_id)

    node = dag[t_sink]
    node.kind = _SEIDEL_DAG_X
    node.key = v
    node.left = traps[tl_id].sink
    node.right = traps[tr_id].sink

    m.vertex_inserted[v] = true
    return (tl_id, tr_id)
end

function _seidel_wire_left_endpoint!(
    m::_SeidelMap,
    t_id::Int,
    t_top::Int,
    t_bot::Int,
    t_l0::Int,
    t_l1::Int,
    u_id::Int,
    l_id::Int,
    p::Int,
)
    traps = m.traps
    if t_l0 != 0 && t_l1 != 0
        traps[u_id].l0 = t_l0
        traps[l_id].l0 = t_l1
        _seidel_replace_right!(traps, t_l0, t_id, u_id)
        _seidel_replace_right!(traps, t_l1, t_id, l_id)
        return nothing
    end

    lnb = t_l0 != 0 ? t_l0 : t_l1
    if t_bot > 0 && m.seg_p[t_bot] == p
        traps[u_id].l0 = lnb
        _seidel_replace_right!(traps, lnb, t_id, u_id)
    elseif t_top > 0 && m.seg_p[t_top] == p
        traps[l_id].l0 = lnb
        _seidel_replace_right!(traps, lnb, t_id, l_id)
    else
        traps[u_id].l0 = lnb
        traps[l_id].l0 = lnb
        if lnb != 0
            traps[lnb].r0 = u_id
            traps[lnb].r1 = l_id
        end
    end
    return nothing
end

function _seidel_wire_right_endpoint!(
    m::_SeidelMap,
    t_id::Int,
    t_top::Int,
    t_bot::Int,
    t_r0::Int,
    t_r1::Int,
    above_cur::Int,
    below_cur::Int,
    q::Int,
)
    traps = m.traps
    if t_r0 != 0 && t_r1 != 0
        traps[above_cur].r0 = t_r0
        traps[below_cur].r0 = t_r1
        _seidel_replace_left!(traps, t_r0, t_id, above_cur)
        _seidel_replace_left!(traps, t_r1, t_id, below_cur)
        return nothing
    end

    rnb = t_r0 != 0 ? t_r0 : t_r1
    if t_bot > 0 && m.seg_q[t_bot] == q
        traps[above_cur].r0 = rnb
        _seidel_replace_left!(traps, rnb, t_id, above_cur)
    elseif t_top > 0 && m.seg_q[t_top] == q
        traps[below_cur].r0 = rnb
        _seidel_replace_left!(traps, rnb, t_id, below_cur)
    else
        traps[above_cur].r0 = rnb
        traps[below_cur].r0 = rnb
        if rnb != 0
            traps[rnb].l0 = above_cur
            traps[rnb].l1 = below_cur
        end
    end
    return nothing
end

function _seidel_add_segment!(m::_SeidelMap, s::Int, hit::Vector{Int})
    p = m.seg_p[s]
    q = m.seg_q[s]
    t0 = _seidel_locate_start(m, m.edge_start_node[s], p, q)
    _seidel_collect_hit!(m, hit, t0, p, q)

    if !m.vertex_inserted[p]
        _, tr_id = _seidel_split_vertical!(m, hit[1], p)
        hit[1] = tr_id
    end
    if !m.vertex_inserted[q]
        tl_id, _ = _seidel_split_vertical!(m, hit[end], q)
        hit[end] = tl_id
    end

    traps = m.traps
    dag = m.dag
    pts = m.pts
    kernel = m.kernel

    above_cur = 0
    below_cur = 0
    k_len = length(hit)

    for k in 1:k_len
        t_id = hit[k]
        T = traps[t_id]
        T.alive = false

        t_leftp = T.leftp
        t_top = T.top
        t_bot = T.bot
        t_l0 = T.l0
        t_l1 = T.l1
        t_r0 = T.r0
        t_r1 = T.r1

        if k == 1
            u_id = _seidel_new_trap!(traps, dag, p, 0, t_top, s)
            l_id = _seidel_new_trap!(traps, dag, p, 0, s, t_bot)
            _seidel_wire_left_endpoint!(m, t_id, t_top, t_bot, t_l0, t_l1, u_id, l_id, p)
            above_cur = u_id
            below_cur = l_id
        else
            prev_id = hit[k - 1]
            prev_T = traps[prev_id]
            v = t_leftp
            other_r = (prev_T.r0 == t_id) ? prev_T.r1 : prev_T.r0
            other_l = (t_l0 == prev_id) ? t_l1 : t_l0

            if orient_sign(kernel, pts[p], pts[q], pts[v]) < 0
                traps[below_cur].rightp = v
                u_id = above_cur
                l_id = _seidel_new_trap!(traps, dag, v, 0, s, t_bot; l0 = below_cur, l1 = other_l)

                traps[below_cur].r0 = l_id
                traps[below_cur].r1 = other_r
                _seidel_replace_left!(traps, other_r, prev_id, below_cur)
                _seidel_replace_right!(traps, other_l, t_id, l_id)
                below_cur = l_id
            else
                traps[above_cur].rightp = v
                l_id = below_cur
                u_id = _seidel_new_trap!(traps, dag, v, 0, t_top, s; l0 = other_l, l1 = above_cur)

                traps[above_cur].r0 = other_r
                traps[above_cur].r1 = u_id
                _seidel_replace_left!(traps, other_r, prev_id, above_cur)
                _seidel_replace_right!(traps, other_l, t_id, u_id)
                above_cur = u_id
            end
        end

        if k == k_len
            traps[above_cur].rightp = q
            traps[below_cur].rightp = q
            _seidel_wire_right_endpoint!(m, t_id, t_top, t_bot, t_r0, t_r1, above_cur, below_cur, q)
        end

        node = dag[T.sink]
        node.kind = _SEIDEL_DAG_Y
        node.key = s
        node.left = traps[u_id].sink
        node.right = traps[l_id].sink
    end

    m.edge_inserted[s] = true
    return nothing
end

# Walk along uninserted edge from `v` (inside or on boundary of `start_tr`) to `nxt`,
# returning a trapezoid whose closure contains `nxt`.
function _seidel_walk_uninserted(m::_SeidelMap, start_tr::Int, v::Int, nxt::Int)::Int
    kernel = m.kernel
    pts = m.pts
    traps = m.traps
    cur = start_tr
    steps = 0
    limit = length(traps)

    if compare_xy(kernel, pts[v], pts[nxt]) < 0
        while _seidel_cmp_pt(kernel, pts, nxt, traps[cur].rightp) > 0
            t = traps[cur]
            if t.r0 != 0 && t.r1 != 0
                cur = orient_sign(kernel, pts[v], pts[nxt], pts[t.rightp]) < 0 ? t.r0 : t.r1
            else
                cur = t.r0 != 0 ? t.r0 : t.r1
            end
            cur = _seidel_step(cur, steps += 1, limit)
        end
    else
        while _seidel_cmp_pt(kernel, pts, nxt, traps[cur].leftp) < 0
            t = traps[cur]
            if t.l0 != 0 && t.l1 != 0
                cur = orient_sign(kernel, pts[nxt], pts[v], pts[t.leftp]) < 0 ? t.l0 : t.l1
            else
                cur = t.l0 != 0 ? t.l0 : t.l1
            end
            cur = _seidel_step(cur, steps += 1, limit)
        end
    end
    return cur
end

# Find a trapezoid adjacent to already-inserted edge `e` touching vertex `v` (at most 2 hops).
function _seidel_find_edge_trap_at_vertex(m::_SeidelMap, tr::Int, v::Int, e::Int)::Int
    traps = m.traps
    t0 = traps[tr]
    (t0.bot == e || t0.top == e) && return tr

    for nb1 in (t0.l0, t0.l1, t0.r0, t0.r1)
        nb1 == 0 && continue
        t1 = traps[nb1]
        (t1.leftp == v || t1.rightp == v) || continue
        (t1.bot == e || t1.top == e) && return nb1
        for nb2 in (t1.l0, t1.l1, t1.r0, t1.r1)
            (nb2 == 0 || nb2 == tr) && continue
            t2 = traps[nb2]
            (t2.leftp == v || t2.rightp == v) || continue
            (t2.bot == e || t2.top == e) && return nb2
        end
    end
    _seidel_inconsistent("no trapezoid adjacent to inserted edge $e at vertex $v")
end

function _seidel_walk_inserted(m::_SeidelMap, tr::Int, v::Int, nxt::Int, e::Int)::Int
    traps = m.traps
    cur = _seidel_find_edge_trap_at_vertex(m, tr, v, e)
    on_bot = (traps[cur].bot == e)
    steps = 0
    limit = length(traps)
    if compare_xy(m.kernel, m.pts[v], m.pts[nxt]) < 0
        while traps[cur].rightp != nxt
            t = traps[cur]
            if on_bot
                cur = (t.r0 != 0 && traps[t.r0].bot == e) ? t.r0 : t.r1
                cur = _seidel_step(cur, steps += 1, limit)
            else
                cur = (t.r0 != 0 && traps[t.r0].top == e) ? t.r0 : t.r1
                cur = _seidel_step(cur, steps += 1, limit)
            end
        end
    else
        while traps[cur].leftp != nxt
            t = traps[cur]
            if on_bot
                cur = (t.l0 != 0 && traps[t.l0].bot == e) ? t.l0 : t.l1
                cur = _seidel_step(cur, steps += 1, limit)
            else
                cur = (t.l0 != 0 && traps[t.l0].top == e) ? t.l0 : t.l1
                cur = _seidel_step(cur, steps += 1, limit)
            end
        end
    end
    return cur
end

# Trace the entire polygon boundary in the current trapezoidal map in O(n) expected time,
# refreshing `edge_start_node[e]` for every uninserted polygon edge `e`.
function _seidel_trace_polygon!(m::_SeidelMap)
    n = m.n
    pts = m.pts
    kernel = m.kernel
    traps = m.traps

    tr = _seidel_locate_vertex_trap(m, 1)
    for i in 1:n
        nxt = mod1(i + 1, n)
        if m.edge_inserted[i]
            tr = _seidel_walk_inserted(m, tr, i, nxt, i)
            continue
        end

        start_tr = tr
        if m.vertex_inserted[i]
            # Only e_prev is inserted; choose the trapezoid at vertex i containing edge i
            u = mod1(i - 1, n)
            t = traps[tr]
            if compare_xy(kernel, pts[u], pts[i]) < 0
                # e_prev came from the left into i
                t_right = t.leftp == i ? tr : (t.r0 != 0 ? t.r0 : t.r1)
                if compare_xy(kernel, pts[i], pts[nxt]) < 0
                    start_tr = t_right
                else
                    tr_obj = traps[t_right]
                    start_tr = orient_sign(kernel, pts[u], pts[i], pts[nxt]) > 0 ? tr_obj.l0 : tr_obj.l1
                end
            else
                # e_prev went right from i
                t_left = t.rightp == i ? tr : (t.l0 != 0 ? t.l0 : t.l1)
                if compare_xy(kernel, pts[nxt], pts[i]) < 0
                    start_tr = t_left
                else
                    tl_obj = traps[t_left]
                    start_tr = orient_sign(kernel, pts[i], pts[u], pts[nxt]) > 0 ? tl_obj.r0 : tl_obj.r1
                end
            end
        end

        ltr = (m.seg_p[i] == i)
        if ltr
            m.edge_start_node[i] = traps[start_tr].sink
        end
        tr = _seidel_walk_uninserted(m, start_tr, i, nxt)
        if !ltr
            m.edge_start_node[i] = traps[tr].sink
        end
    end
    return nothing
end

# Phase schedule N(h) = floor(n / log^(h)(n)) yielding O(log* n) phases.
function _seidel_phase_bounds(n::Int)
    bounds = Int[]
    d = log2(Float64(n))
    prev = 0
    while d > 1.5
        b = clamp(floor(Int, n / d), 1, n - 1)
        if b > prev
            push!(bounds, b)
            prev = b
        end
        d = log2(d)
    end
    push!(bounds, n)
    return bounds
end

function _seidel_build_trapezoidal_map(
    kernel::AbstractKernel,
    pts,
    rng::AbstractRNG,
)
    n = length(pts)
    m = _seidel_init_map(kernel, pts)
    perm = randperm(rng, n)
    bounds = _seidel_phase_bounds(n)
    hit = Int[]
    sizehint!(hit, 32)

    start_idx = 1
    num_phases = length(bounds)
    for phase in 1:num_phases
        stop_idx = bounds[phase]
        for idx in start_idx:stop_idx
            _seidel_add_segment!(m, perm[idx], hit)
        end
        if phase < num_phases
            _seidel_trace_polygon!(m)
        end
        start_idx = stop_idx + 1
    end
    return m
end

# Extract unique internal diagonals (leftp, rightp) from interior trapezoids in O(n) time.
function _seidel_extract_diagonals(m::_SeidelMap)
    n = m.n
    diags = Tuple{Int,Int}[]
    for T in m.traps
        T.alive || continue
        (T.top >= 1 && T.bot >= 1) || continue
        is_interior = m.seg_ltr[T.bot] && !m.seg_ltr[T.top]
        is_interior || continue
        u = T.leftp
        v = T.rightp
        d = abs(u - v)
        if d > 1 && d < n - 1
            push!(diags, u < v ? (u, v) : (v, u))
        end
    end
    return diags
end

# Partition polygon by non-crossing diagonals into CCW monotone pieces in O(n) time
# using cyclic vertex-index ordering of half-edges (no angular sorting required).
function _seidel_extract_pieces(n::Int, diags::Vector{Tuple{Int,Int}})
    d = length(diags)
    d == 0 && return [collect(1:n)]

    total_he = 2 * n + 2 * d
    he_from = Vector{Int}(undef, total_he)
    he_to = Vector{Int}(undef, total_he)
    he_twin = Vector{Int}(undef, total_he)

    for i in 1:n
        nxt = mod1(i + 1, n)
        he_from[i] = i
        he_to[i] = nxt
        he_twin[i] = n + i
        he_from[n + i] = nxt
        he_to[n + i] = i
        he_twin[n + i] = i
    end

    diag_by_dest = [Int[] for _ in 1:n]
    for k in 1:d
        u, v = diags[k]
        h1 = 2 * n + 2 * k - 1
        h2 = h1 + 1
        he_from[h1] = u
        he_to[h1] = v
        he_twin[h1] = h2
        he_from[h2] = v
        he_to[h2] = u
        he_twin[h2] = h1
        push!(diag_by_dest[v], h1)
        push!(diag_by_dest[u], h2)
    end

    adj_sorted = [Int[] for _ in 1:n]
    for dest in 1:n
        for h in diag_by_dest[dest]
            push!(adj_sorted[he_from[h]], h)
        end
    end

    next_he = Vector{Int}(undef, total_he)
    ccw_list = Int[]
    for v in 1:n
        empty!(ccw_list)
        push!(ccw_list, v) # forward CCW boundary edge v -> mod1(v + 1, n)
        for h in adj_sorted[v]
            he_to[h] > v && push!(ccw_list, h)
        end
        for h in adj_sorted[v]
            he_to[h] < v && push!(ccw_list, h)
        end
        rev_in = n + mod1(v - 1, n) # reverse edge v -> mod1(v - 1, n)
        push!(ccw_list, rev_in)

        m_v = length(ccw_list)
        for idx in 1:m_v
            h_out = ccw_list[idx]
            h_pred = idx == 1 ? ccw_list[m_v] : ccw_list[idx - 1]
            next_he[he_twin[h_out]] = h_pred
        end
    end

    visited = falses(total_he)
    pieces = Vector{Int}[]
    for start_h in 1:total_he
        (n < start_h <= 2 * n) && continue # skip exterior reverse boundary edges
        visited[start_h] && continue
        pc = Int[]
        cur = start_h
        while !visited[cur]
            visited[cur] = true
            push!(pc, he_from[cur])
            cur = next_he[cur]
        end
        push!(pieces, pc)
    end
    return pieces
end

# Triangulate a CCW piece monotone under `compare_xy` in O(k) time, emitting
# original-index diagonals and CCW triangles.
function _seidel_triangulate_monotone_piece!(
    kernel::AbstractKernel,
    pts,
    pc::Vector{Int},
    all_diagonals::Vector{Tuple{Int,Int}},
    all_triangles::Vector{NTuple{3,Int}},
)
    k = length(pc)
    k < 3 && return nothing
    if k == 3
        push!(all_triangles, (pc[1], pc[2], pc[3]))
        return nothing
    end

    lo = 1
    hi = 1
    for i in 2:k
        compare_xy(kernel, pts[pc[i]], pts[pc[lo]]) < 0 && (lo = i)
        compare_xy(kernel, pts[pc[i]], pts[pc[hi]]) > 0 && (hi = i)
    end

    c_minus = Int[lo]
    curr = lo
    while curr != hi
        curr = mod1(curr - 1, k)
        push!(c_minus, curr)
    end

    c_plus = Int[lo]
    curr = lo
    while curr != hi
        curr = mod1(curr + 1, k)
        push!(c_plus, curr)
    end

    if orient_sign(kernel, pts[pc[lo]], pts[pc[c_plus[2]]], pts[pc[c_minus[2]]]) > 0
        upper, lower = c_minus, c_plus
    else
        upper, lower = c_plus, c_minus
    end

    is_upper = falses(k)
    for v in @view upper[2:end-1]
        is_upper[v] = true
    end

    order = Int[]
    sizehint!(order, k)
    i = 1
    j = 1
    while i <= length(upper) || j <= length(lower)
        if j > length(lower) || (
            i <= length(upper) &&
            compare_xy(kernel, pts[pc[upper[i]]], pts[pc[lower[j]]]) <= 0
        )
            v = upper[i]
            i += 1
        else
            v = lower[j]
            j += 1
        end
        (isempty(order) || last(order) != v) && push!(order, v)
    end

    @inline function emit_diag!(a::Int, b::Int)
        if a != b && !_is_boundary_edge(a, b, k)
            oa = pc[a]
            ob = pc[b]
            push!(all_diagonals, oa < ob ? (oa, ob) : (ob, oa))
        end
    end

    # The orientation of every emitted triangle follows from the chain structure (the lower
    # chain runs left-to-right and the upper chain right-to-left in CCW order), so callers pass
    # the vertices already in counterclockwise order and no predicate is re-evaluated here.
    @inline emit_tri!(a::Int, b::Int, c::Int) = push!(all_triangles, (pc[a], pc[b], pc[c]))

    stack = Int[order[1], order[2]]
    for idx in 3:k-1
        v = order[idx]
        if is_upper[v] != is_upper[last(stack)]
            while length(stack) >= 2
                s_top = pop!(stack)
                s_next = last(stack)
                is_upper[v] ? emit_tri!(v, s_next, s_top) : emit_tri!(v, s_top, s_next)
                emit_diag!(v, s_top)
            end
            pop!(stack)
            push!(stack, order[idx - 1], v)
        else
            q = pop!(stack)
            while !isempty(stack)
                a = last(stack)
                turn = orient_sign(kernel, pts[pc[a]], pts[pc[q]], pts[pc[v]])
                (is_upper[v] ? turn < 0 : turn > 0) || break
                is_upper[v] ? emit_tri!(a, v, q) : emit_tri!(a, q, v)
                emit_diag!(a, v)
                q = pop!(stack)
            end
            push!(stack, q, v)
        end
    end

    v_last = last(order)
    while length(stack) >= 2
        s_top = pop!(stack)
        s_next = last(stack)
        is_upper[s_top] ? emit_tri!(v_last, s_top, s_next) : emit_tri!(v_last, s_next, s_top)
        emit_diag!(v_last, s_top)
    end
    return nothing
end

function _triangulate_from_partition(
    kernel::AbstractKernel,
    vertices,
    partition_diags::Vector{Tuple{Int,Int}},
)
    n = length(vertices)
    pieces = _seidel_extract_pieces(n, partition_diags)

    all_diagonals = copy(partition_diags)
    sizehint!(all_diagonals, n - 3)
    all_triangles = NTuple{3,Int}[]
    sizehint!(all_triangles, n - 2)

    for pc in pieces
        _seidel_triangulate_monotone_piece!(kernel, vertices, pc, all_diagonals, all_triangles)
    end

    return PolygonTriangulation(all_diagonals, all_triangles)
end

"""
    triangulate_seidel(vertices; kernel=FloatingKernel(), rng=Random.default_rng()) -> PolygonTriangulation

Triangulate a simple counterclockwise polygon using Seidel's randomized incremental algorithm
in expected `O(n log* n)` time.

Builds the trapezoidal decomposition of the polygon edges incrementally in `O(log* n)`
phases interleaved with `O(n)` polygon boundary tracing, decomposes the interior into
monotone pieces via trapezoid diagonals in `O(n)` time, and triangulates each piece in
linear time `O(k)`.

# Preconditions
- `vertices` must form a simple counterclockwise polygon with `n >= 3` vertices.
"""
function triangulate_seidel(
    vertices;
    kernel::AbstractKernel = FloatingKernel(),
    rng::AbstractRNG = Random.default_rng(),
)
    n = _validate_ccw_polygon(kernel, vertices)
    if n == 3
        return PolygonTriangulation(Tuple{Int,Int}[], NTuple{3,Int}[(1, 2, 3)])
    end

    trap_map = _seidel_build_trapezoidal_map(kernel, vertices, rng)
    trap_diags = _seidel_extract_diagonals(trap_map)
    return _triangulate_from_partition(kernel, vertices, trap_diags)
end

export triangulate_seidel
