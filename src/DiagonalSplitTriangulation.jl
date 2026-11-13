# Diagonal-split triangulation (prototype): recursively split the polygon by diagonals
# between reflex vertices and fan the resulting pieces that have at most one reflex vertex.
# A random pair of reflex vertices (half the piece's reflex count apart) is tried first;
# when it fails, a reflex endpoint is resolved by ray shooting (first boundary hit, then the
# angularly first reflex vertex in the swept triangle). Every diagonal is validated by a cone
# test at both endpoints and a Kd-tree over the original edges.
# See `plans/diagonal_split.md`.

# ---------------------------------------------------------------------------------------
# Edge Kd-tree. Built over edge midpoints; node and edge boxes are outward-rounded Float64
# boxes covering the true edge extents, so box pruning is conservative for any coordinate
# type. Edge `e` joins vertices `e` and `mod1(e + 1, n)`. The tree holds only the original
# edges: a candidate diagonal of a piece that crosses an earlier diagonal must also cross an
# original edge to get back into the piece, so diagonals never need to be indexed.
# ---------------------------------------------------------------------------------------

struct _EdgeTree
    left::Vector{Int32}
    right::Vector{Int32}
    rs::Vector{Int32}
    re::Vector{Int32}
    x0::Vector{Float64}
    x1::Vector{Float64}
    y0::Vector{Float64}
    y1::Vector{Float64}
    eid::Vector{Int32}       # slot -> edge id
    ex0::Vector{Float64}     # per-slot edge boxes
    ex1::Vector{Float64}
    ey0::Vector{Float64}
    ey1::Vector{Float64}
    tol::Float64             # absolute slack for the Float64 line-side test
    stack::Vector{Int32}
end

@inline _dn(x) = prevfloat(Float64(x))
@inline _up(x) = nextfloat(Float64(x))

function _build_edge_tree(poly, n::Int; leaf_size::Int = 8)
    px = Vector{Float64}(undef, n)
    py = Vector{Float64}(undef, n)
    cmax = 0.0
    @inbounds for i in 1:n
        a = poly[i]
        b = poly[i == n ? 1 : i + 1]
        ax = Float64(a.x); ay = Float64(a.y); bx = Float64(b.x); by = Float64(b.y)
        px[i] = 0.5 * ax + 0.5 * bx
        py[i] = 0.5 * ay + 0.5 * by
        cmax = max(cmax, abs(ax), abs(ay))
    end
    ids = Int32[Int32(i) for i in 1:n]
    ft = Kd_tree._build_flat_tree!(Kd_tree.FlatNode2D{Float64}, px, py, ids, n, leaf_size)
    nodes = ft.nodes
    nn = length(nodes)
    left = Vector{Int32}(undef, nn); right = Vector{Int32}(undef, nn)
    rs = Vector{Int32}(undef, nn); re = Vector{Int32}(undef, nn)
    x0 = fill(Inf, nn); x1 = fill(-Inf, nn); y0 = fill(Inf, nn); y1 = fill(-Inf, nn)
    m = length(ft.perm)
    eid = copy(ft.perm)
    ex0 = Vector{Float64}(undef, m); ex1 = similar(ex0); ey0 = similar(ex0); ey1 = similar(ex0)
    @inbounds for k in 1:m
        e = Int(eid[k])
        a = poly[e]
        b = poly[e == n ? 1 : e + 1]
        ex0[k] = _dn(min(a.x, b.x)); ex1[k] = _up(max(a.x, b.x))
        ey0[k] = _dn(min(a.y, b.y)); ey1[k] = _up(max(a.y, b.y))
    end
    # Children always have larger ids than their parent: refit bottom-up in reverse order.
    @inbounds for id in nn:-1:1
        nd = nodes[id]
        left[id] = nd.left; right[id] = nd.right
        rs[id] = nd.r_start; re[id] = nd.r_end
        if nd.is_leaf
            left[id] = 0; right[id] = 0
            for k in nd.r_start:nd.r_end
                x0[id] = min(x0[id], ex0[k]); x1[id] = max(x1[id], ex1[k])
                y0[id] = min(y0[id], ey0[k]); y1[id] = max(y1[id], ey1[k])
            end
        else
            for c in (nd.left, nd.right)
                x0[id] = min(x0[id], x0[c]); x1[id] = max(x1[id], x1[c])
                y0[id] = min(y0[id], y0[c]); y1[id] = max(y1[id], y1[c])
            end
        end
    end
    # The computed line-side value has absolute error far below 1e-12 * C^2 for
    # coordinates (and box corners) bounded by C, including Float64 conversion error.
    c = _up(cmax) + 1.0e-300
    tol = 1.0e-12 * (2c) * (2c) + 1.0e-300
    stack = Int32[]
    sizehint!(stack, 64)
    return _EdgeTree(left, right, rs, re, x0, x1, y0, y1, eid, ex0, ex1, ey0, ey1, tol, stack)
end

# Orientation sign used by all decisions of this algorithm. Under the plain Float64 kernel,
# values within Shewchuk's static error bound (orient2dfast, errboundA ~ 3.3e-16) are
# reported as 0. Every caller treats 0 conservatively (a diagonal touching an edge, a vertex
# being reflex, a point outside an open cone), so rounding can only reject a valid diagonal
# or delay a fan, never accept an invalid one.
@inline _cons_orient(kernel, a, b, c) = orient_sign(kernel, a, b, c)
@inline function _cons_orient(::FloatingKernel{Float64}, a, b, c)
    ax = Float64(a.x); ay = Float64(a.y)
    l = (Float64(b.x) - ax) * (Float64(c.y) - ay)
    r = (Float64(b.y) - ay) * (Float64(c.x) - ax)
    det = l - r
    bound = 1.0e-15 * (abs(l) + abs(r))
    return det > bound ? 1 : (det < -bound ? -1 : 0)
end

# `p` lies on the closed segment [a, b] (conservatively: near-collinear counts).
@inline _cons_on_seg(kernel, p, a, b) =
    _cons_orient(kernel, a, b, p) == 0 &&
    min(a.x, b.x) <= p.x <= max(a.x, b.x) && min(a.y, b.y) <= p.y <= max(a.y, b.y)

# Closed segments [p1,p2] and [q1,q2] share a point (conservatively under Float64).
# For (near-)collinear segments, bounding-box overlap decides: exact for truly collinear
# segments, conservative otherwise.
@inline function _closed_segs_meet(kernel, p1, p2, q1, q2)::Bool
    o1 = _cons_orient(kernel, p1, p2, q1)
    o2 = _cons_orient(kernel, p1, p2, q2)
    o1 * o2 > 0 && return false
    if o1 == 0 && o2 == 0
        return max(p1.x, p2.x) >= min(q1.x, q2.x) && max(q1.x, q2.x) >= min(p1.x, p2.x) &&
               max(p1.y, p2.y) >= min(q1.y, q2.y) && max(q1.y, q2.y) >= min(p1.y, p2.y)
    end
    o3 = _cons_orient(kernel, q1, q2, p1)
    o4 = _cons_orient(kernel, q1, q2, p2)
    return o3 * o4 <= 0
end

# Does segment (s1, s2) block the candidate diagonal (a, b)? Segments sharing an endpoint
# with it block only if their other endpoint lies on it.
@inline function _seg_blocks(kernel, poly, A, B, a::Int, b::Int, s1::Int, s2::Int)::Bool
    if s1 == a || s1 == b || s2 == a || s2 == b
        far = (s1 == a || s1 == b) ? s2 : s1
        (far == a || far == b) && return false
        return _cons_on_seg(kernel, poly[far], A, B)
    end
    return _closed_segs_meet(kernel, A, B, poly[s1], poly[s2])
end

# Does some original edge meet the closed segment between vertices `a` and `b` other than at
# `a` or `b` themselves? Edges incident to `a` or `b` count only if their other endpoint lies
# on the segment.
function _edge_blocked(et::_EdgeTree, poly, kernel, n::Int, a::Int, b::Int, stats)::Bool
    A = poly[a]
    B = poly[b]
    axf = Float64(A.x); ayf = Float64(A.y); bxf = Float64(B.x); byf = Float64(B.y)
    sx0 = _dn(min(A.x, B.x)); sx1 = _up(max(A.x, B.x))
    sy0 = _dn(min(A.y, B.y)); sy1 = _up(max(A.y, B.y))
    dx = bxf - axf
    dy = byf - ayf
    tol = et.tol
    stack = et.stack
    empty!(stack)
    push!(stack, Int32(1))
    @inbounds while !isempty(stack)
        id = Int(pop!(stack))
        stats[3] += 1
        (et.x1[id] < sx0 || et.x0[id] > sx1 || et.y1[id] < sy0 || et.y0[id] > sy1) && continue
        bx0 = et.x0[id]; bx1 = et.x1[id]; by0 = et.y0[id]; by1 = et.y1[id]
        s1 = dx * (by0 - ayf) - dy * (bx0 - axf)
        s2 = dx * (by0 - ayf) - dy * (bx1 - axf)
        s3 = dx * (by1 - ayf) - dy * (bx0 - axf)
        s4 = dx * (by1 - ayf) - dy * (bx1 - axf)
        ((s1 > tol) & (s2 > tol) & (s3 > tol) & (s4 > tol)) && continue
        ((s1 < -tol) & (s2 < -tol) & (s3 < -tol) & (s4 < -tol)) && continue
        l = et.left[id]
        if l == 0
            for k in et.rs[id]:et.re[id]
                (et.ex1[k] < sx0 || et.ex0[k] > sx1 || et.ey1[k] < sy0 || et.ey0[k] > sy1) &&
                    continue
                e = Int(et.eid[k])
                e2 = e == n ? 1 : e + 1
                stats[4] += 1
                _seg_blocks(kernel, poly, A, B, a, b, e, e2) && return true
            end
        else
            push!(stack, l)
            push!(stack, et.right[id])
        end
    end
    return false
end

# ---------------------------------------------------------------------------------------
# Ray shooting and the angular sweep (plain Float64 heuristics; every diagonal they propose
# is validated exactly by `_diag_ok`).
# ---------------------------------------------------------------------------------------

# Parameter interval of the ray A + t D (t >= 0) inside box [x0,x1] x [y0,y1]: entry t, or
# Inf if the ray misses the box.
@inline function _ray_box_entry(ax, ay, dx, dy, x0, x1, y0, y1)::Float64
    tlo = 0.0
    thi = Inf
    if dx == 0.0
        (ax < x0 || ax > x1) && return Inf
    else
        t1 = (x0 - ax) / dx; t2 = (x1 - ax) / dx
        tlo = max(tlo, min(t1, t2)); thi = min(thi, max(t1, t2))
    end
    if dy == 0.0
        (ay < y0 || ay > y1) && return Inf
    else
        t1 = (y0 - ay) / dy; t2 = (y1 - ay) / dy
        tlo = max(tlo, min(t1, t2)); thi = min(thi, max(t1, t2))
    end
    return tlo <= thi ? tlo : Inf
end

# Ray parameter at which A + t D meets segment [P, Q], or Inf.
@inline function _ray_seg_t(ax, ay, dx, dy, px, py, qx, qy)::Float64
    ex = qx - px; ey = qy - py
    den = dx * ey - dy * ex
    den == 0.0 && return Inf
    wx = px - ax; wy = py - ay
    t = (wx * ey - wy * ex) / den
    u = (wx * dy - wy * dx) / den
    (t > 0.0 && -1.0e-12 <= u <= 1.0 + 1.0e-12) || return Inf
    return t
end

# First original edge hit by the ray from vertex `a` in direction D, ignoring edges incident
# to `a` and only accepting hits before `best_t`. Returns (s1, s2, t); s1 == 0 if none.
function _ray_first_hit(et::_EdgeTree, poly, n::Int, a::Int, ax, ay, dx, dy, stats,
                        best_t::Float64 = Inf, b1::Int = 0, b2::Int = 0)
    stack = et.stack
    empty!(stack)
    push!(stack, Int32(1))
    @inbounds while !isempty(stack)
        id = Int(pop!(stack))
        stats[3] += 1
        _ray_box_entry(ax, ay, dx, dy, et.x0[id], et.x1[id], et.y0[id], et.y1[id]) < best_t ||
            continue
        l = et.left[id]
        if l == 0
            for k in et.rs[id]:et.re[id]
                _ray_box_entry(ax, ay, dx, dy, et.ex0[k], et.ex1[k], et.ey0[k], et.ey1[k]) <
                    best_t || continue
                e = Int(et.eid[k])
                e2 = e == n ? 1 : e + 1
                (e == a || e2 == a) && continue
                P = poly[e]; Q = poly[e2]
                t = _ray_seg_t(ax, ay, dx, dy, Float64(P.x), Float64(P.y), Float64(Q.x), Float64(Q.y))
                t < best_t && (best_t = t; b1 = e; b2 = e2)
            end
        else
            push!(stack, l)
            push!(stack, et.right[id])
        end
    end
    return b1, b2, best_t
end

@inline _o2(ax, ay, bx, by, cx, cy) = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)

# All four corners of the node box have sg * orient(P, Q, corner) < 0.
@inline _box_all_neg(sg, px, py, qx, qy, nd) =
    sg * _o2(px, py, qx, qy, nd.x0, nd.y0) < 0.0 && sg * _o2(px, py, qx, qy, nd.x1, nd.y0) < 0.0 &&
    sg * _o2(px, py, qx, qy, nd.x0, nd.y1) < 0.0 && sg * _o2(px, py, qx, qy, nd.x1, nd.y1) < 0.0

# Among the reflex vertices (and `xid`) in the closed triangle (A, H, X), the one reached
# first when rotating the ray A->H towards X (ties: nearest to A). If the ray's first hit
# was on segment (.., X), that vertex is visible from `a`.
function _first_swept_reflex(rt, rstack::Vector{Int32}, a::Int, ax, ay, hx, hy, xx, xy, xid::Int)
    sg = _o2(ax, ay, hx, hy, xx, xy)
    sg == 0.0 && return xid
    sg = sign(sg)
    best = xid; bx = xx; by = xy
    isempty(rt.nodes) && return best
    tx0 = min(ax, hx, xx); tx1 = max(ax, hx, xx); ty0 = min(ay, hy, xy); ty1 = max(ay, hy, xy)
    nodes = rt.nodes
    empty!(rstack)
    push!(rstack, Int32(1))
    @inbounds while !isempty(rstack)
        nd = nodes[Int(pop!(rstack))]
        (nd.x1 < tx0 || nd.x0 > tx1 || nd.y1 < ty0 || nd.y0 > ty1) && continue
        # prune: whole box strictly outside one triangle side, or strictly beyond the best ray
        (_box_all_neg(sg, ax, ay, hx, hy, nd) || _box_all_neg(sg, hx, hy, xx, xy, nd) ||
         _box_all_neg(sg, xx, xy, ax, ay, nd) || _box_all_neg(-sg, ax, ay, bx, by, nd)) && continue
        if nd.is_leaf
            for k in nd.r_start:nd.r_end
                v = Int(rt.perm[k])
                v == a && continue
                yx = rt.px[k]; yy = rt.py[k]
                (sg * _o2(ax, ay, hx, hy, yx, yy) >= 0.0 && sg * _o2(hx, hy, xx, xy, yx, yy) >= 0.0 &&
                 sg * _o2(xx, xy, ax, ay, yx, yy) >= 0.0) || continue
                o = sg * _o2(ax, ay, bx, by, yx, yy)
                if o < 0.0 ||
                   (o == 0.0 && (yx - ax)^2 + (yy - ay)^2 < (bx - ax)^2 + (by - ay)^2)
                    best = v; bx = yx; by = yy
                end
            end
        else
            push!(rstack, nd.left)
            push!(rstack, nd.right)
        end
    end
    return best
end

# ---------------------------------------------------------------------------------------
# Pieces and the halving search
# ---------------------------------------------------------------------------------------

mutable struct _DiagSplitState{P,K,R}
    poly::P
    kernel::K
    n::Int
    et::_EdgeTree
    rt::Kd_tree.FlatTree2D{Float64}   # static tree over the input's reflex vertices
    rstack::Vector{Int32}
    rng::R
    diagonals::Vector{Tuple{Int,Int}}
    triangles::Vector{NTuple{3,Int}}
    # stats: 1 diagonal tests, 2 failed tests, 3 edge-tree nodes, 4 exact edge tests,
    #        5 splits, 6 fallback pieces, 7 fanned pieces, 8 tests made by galloping,
    #        9 ray resolutions tried, 10 ray resolutions that gave a diagonal,
    #        11 total size of the pieces that were split
    stats::Vector{Int}
end

# Open interior cone at `a` (piece neighbors `u` before, `w` after; ccw) contains `x`.
# A (near-)straight angle at `a` uses the stricter convex rule.
@inline function _in_open_cone(kernel, u, a, w, x)::Bool
    s1 = _cons_orient(kernel, a, w, x) > 0
    s2 = _cons_orient(kernel, u, a, x) > 0
    return _cons_orient(kernel, u, a, w) < 0 ? (s1 | s2) : (s1 & s2)
end

@inline _is_reflex_at(kernel, poly, u::Int, v::Int, w::Int) =
    _cons_orient(kernel, poly[u], poly[v], poly[w]) <= 0

function _diag_ok(st::_DiagSplitState, verts::Vector{Int32}, p::Int, q::Int)::Bool
    m = length(verts)
    poly = st.poly
    kernel = st.kernel
    a = Int(verts[p]); b = Int(verts[q])
    st.stats[1] += 1
    ua = Int(verts[p == 1 ? m : p - 1]); wa = Int(verts[p == m ? 1 : p + 1])
    ub = Int(verts[q == 1 ? m : q - 1]); wb = Int(verts[q == m ? 1 : q + 1])
    ok = _in_open_cone(kernel, poly[ua], poly[a], poly[wa], poly[b]) &&
         _in_open_cone(kernel, poly[ub], poly[b], poly[wb], poly[a]) &&
         !_edge_blocked(st.et, poly, kernel, st.n, a, b, st.stats)
    ok || (st.stats[2] += 1)
    return ok
end

# Middle reflex position strictly between p and q (cyclically forward), or 0.
function _middle_reflex(refl::Vector{Int32}, m::Int, p::Int, q::Int)::Int
    if p < q
        lo = searchsortedfirst(refl, p + 1)
        hi = searchsortedlast(refl, q - 1)
        cnt = hi - lo + 1
        cnt <= 0 && return 0
        return Int(refl[lo + (cnt - 1) ÷ 2])
    else
        lo1 = searchsortedfirst(refl, p + 1)
        c1 = length(refl) - lo1 + 1
        c2 = searchsortedlast(refl, q - 1)
        cnt = c1 + c2
        cnt <= 0 && return 0
        t = (cnt - 1) ÷ 2
        return Int(t < c1 ? refl[lo1 + t] : refl[t - c1 + 1])
    end
end

# Exponential (galloping) search along a chain of length `d` starting at position `p`
# (the diagonal to `p + d` is known to fail): test `p + 2, p + 4, p + 8, ...` while valid,
# then binary-search between the last valid and first failed offset. Returns the farthest
# valid diagonal found, or (0, 0) if `p + 2` already fails.
function _gallop_chain(st::_DiagSplitState, verts::Vector{Int32}, p::Int, d::Int)
    m = length(verts)
    best = 0
    t = 2
    while t < d
        st.stats[8] += 1
        _diag_ok(st, verts, p, mod1(p + t, m)) || break
        best = t
        t *= 2
    end
    best == 0 && return (0, 0)
    lo = best
    hi = min(2best, d)
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        st.stats[8] += 1
        if _diag_ok(st, verts, p, mod1(p + mid, m))
            lo = mid
        else
            hi = mid
        end
    end
    return (p, mod1(p + lo, m))
end

# Search the chain p -> q (cyclically forward positions) for a valid diagonal.
function _try_chain(st::_DiagSplitState, verts::Vector{Int32}, refl::Vector{Int32}, p::Int, q::Int)
    m = length(verts)
    d = mod(q - p, m)
    d <= 1 && return (0, 0)
    if d != m - 1 && _diag_ok(st, verts, p, q)
        return (p, q)
    end
    mid = _middle_reflex(refl, m, p, q)
    if mid == 0
        # Convex chain p -> q: gallop forward from p for the farthest visible vertex.
        r = _gallop_chain(st, verts, p, d)
        r[1] != 0 && return r
        mid = mod1(p + d ÷ 2, m)
    end
    r = _try_chain(st, verts, refl, p, mid)
    r[1] != 0 && return r
    return _try_chain(st, verts, refl, mid, q)
end

# Resolve the reflex vertex at position `p` of the piece: shoot a ray along the interior
# bisector, take the first boundary segment hit, and propose the angularly first reflex
# vertex of the swept triangle (trying both endpoints of the hit segment). Returns the
# validated diagonal as positions, or (0, 0).
function _ray_resolve(st::_DiagSplitState, verts::Vector{Int32}, p::Int)
    m = length(verts)
    poly = st.poly
    st.stats[9] += 1
    a = Int(verts[p])
    U = poly[Int(verts[p == 1 ? m : p - 1])]; A = poly[a]; W = poly[Int(verts[p == m ? 1 : p + 1])]
    ax = Float64(A.x); ay = Float64(A.y)
    ux = Float64(U.x) - ax; uy = Float64(U.y) - ay
    wx = Float64(W.x) - ax; wy = Float64(W.y) - ay
    lu = hypot(ux, uy); lw = hypot(wx, wy)
    dx = -(ux / lu + wx / lw); dy = -(uy / lu + wy / lw)
    if hypot(dx, dy) < 1.0e-9   # straight angle: left normal of U -> W
        dx = -(wy - uy); dy = wx - ux
    end
    # The piece's own diagonals (consecutive piece vertices that are not consecutive in the
    # input) bound it too: scan them, then let the edge tree look only for nearer hits.
    n = st.n
    best_t = Inf; s1 = 0; s2 = 0
    @inbounds for t in 1:m
        v1 = Int(verts[t]); v2 = Int(verts[t == m ? 1 : t + 1])
        (v2 == v1 + 1 || (v1 == n && v2 == 1)) && continue
        (v1 == a || v2 == a) && continue
        P = poly[v1]; Q = poly[v2]
        tt = _ray_seg_t(ax, ay, dx, dy, Float64(P.x), Float64(P.y), Float64(Q.x), Float64(Q.y))
        tt < best_t && (best_t = tt; s1 = v1; s2 = v2)
    end
    s1, s2, t = _ray_first_hit(st.et, poly, n, a, ax, ay, dx, dy, st.stats, best_t, s1, s2)
    s1 == 0 && return (0, 0)
    hx = ax + t * dx; hy = ay + t * dy
    P1 = poly[s1]; P2 = poly[s2]
    d1 = (Float64(P1.x) - hx)^2 + (Float64(P1.y) - hy)^2
    d2 = (Float64(P2.x) - hx)^2 + (Float64(P2.y) - hy)^2
    for x in (d1 <= d2 ? (s1, s2) : (s2, s1))
        X = poly[x]
        z = _first_swept_reflex(st.rt, st.rstack, a, ax, ay, hx, hy, Float64(X.x), Float64(X.y), x)
        q = findfirst(==(Int32(z)), verts)
        q === nothing && continue
        d = mod(q - p, m)
        (d <= 1 || d == m - 1) && continue
        if _diag_ok(st, verts, p, q)
            st.stats[10] += 1
            return (p, q)
        end
    end
    return (0, 0)
end

function _fan_piece!(st::_DiagSplitState, verts::Vector{Int32}, s::Int)
    m = length(verts)
    v = Int(verts[s])
    for t in 1:m-2
        w1 = Int(verts[mod1(s + t, m)])
        w2 = Int(verts[mod1(s + t + 1, m)])
        push!(st.triangles, (v, w1, w2))
        t >= 2 && push!(st.diagonals, v < w1 ? (v, w1) : (w1, v))
    end
    st.stats[7] += 1
    return nothing
end

function _fallback_piece!(st::_DiagSplitState, verts::Vector{Int32})
    st.stats[6] += 1
    sub = [st.poly[v] for v in verts]
    res = triangulate_ear_cut_u_kd_tree(sub; kernel = st.kernel)
    for (i, j, k) in res.triangles
        push!(st.triangles, (Int(verts[i]), Int(verts[j]), Int(verts[k])))
    end
    for (i, j) in res.diagonals
        a = Int(verts[i]); b = Int(verts[j])
        push!(st.diagonals, a < b ? (a, b) : (b, a))
    end
    return nothing
end

# Split piece at positions i < j into verts[i:j] and verts[j:m] ∪ verts[1:i].
function _split_piece(st::_DiagSplitState, verts::Vector{Int32}, refl::Vector{Int32}, i::Int, j::Int)
    m = length(verts)
    kernel = st.kernel
    poly = st.poly
    v1 = verts[i:j]
    v2 = Vector{Int32}(undef, m - j + 1 + i)
    copyto!(v2, 1, verts, j, m - j + 1)
    copyto!(v2, m - j + 2, verts, 1, i)
    m1 = length(v1); m2 = length(v2)

    r1 = Int32[]
    _is_reflex_at(kernel, poly, Int(v1[m1]), Int(v1[1]), Int(v1[2])) && push!(r1, Int32(1))
    lo = searchsortedfirst(refl, i + 1); hi = searchsortedlast(refl, j - 1)
    for t in lo:hi
        push!(r1, refl[t] - Int32(i - 1))
    end
    _is_reflex_at(kernel, poly, Int(v1[m1 - 1]), Int(v1[m1]), Int(v1[1])) && push!(r1, Int32(m1))

    r2 = Int32[]
    _is_reflex_at(kernel, poly, Int(v2[m2]), Int(v2[1]), Int(v2[2])) && push!(r2, Int32(1))
    for t in searchsortedfirst(refl, j + 1):length(refl)
        push!(r2, refl[t] - Int32(j - 1))
    end
    for t in 1:searchsortedlast(refl, i - 1)
        push!(r2, refl[t] + Int32(m - j + 1))
    end
    _is_reflex_at(kernel, poly, Int(v2[m2 - 1]), Int(v2[m2]), Int(v2[1])) && push!(r2, Int32(m2))

    a = Int(verts[i]); b = Int(verts[j])
    push!(st.diagonals, a < b ? (a, b) : (b, a))
    st.stats[5] += 1
    return (v1, r1), (v2, r2)
end

function _process_pieces!(st::_DiagSplitState, verts0::Vector{Int32}, refl0::Vector{Int32})
    stack = [(verts0, refl0)]
    while !isempty(stack)
        verts, refl = pop!(stack)
        m = length(verts)
        if m == 3
            push!(st.triangles, (Int(verts[1]), Int(verts[2]), Int(verts[3])))
            continue
        end
        k = length(refl)
        if k == 0
            _fan_piece!(st, verts, 1)
            continue
        elseif k == 1
            _fan_piece!(st, verts, Int(refl[1]))
            continue
        end
        ia = rand(st.rng, 1:k)
        ib = mod1(ia + k ÷ 2, k)
        p = Int(refl[ia]); q = Int(refl[ib])
        d = mod(q - p, m)
        r = (d > 1 && d < m - 1 && _diag_ok(st, verts, p, q)) ? (p, q) : (0, 0)
        r[1] == 0 && (r = _ray_resolve(st, verts, p))
        r[1] == 0 && (r = _ray_resolve(st, verts, q))
        r[1] == 0 && (r = _try_chain(st, verts, refl, p, q))
        r[1] == 0 && (r = _try_chain(st, verts, refl, q, p))
        if r[1] == 0
            _fallback_piece!(st, verts)
            continue
        end
        i, j = minmax(r[1], r[2])
        st.stats[11] += m
        piece1, piece2 = _split_piece(st, verts, refl, i, j)
        push!(stack, piece1)
        push!(stack, piece2)
    end
    return nothing
end

function _triangulate_diagonal_split_impl(vertices, kernel::AbstractKernel, rng::AbstractRNG)
    n = _validate_ccw_polygon(kernel, vertices)
    stats = zeros(Int, 11)
    if n == 3
        return PolygonTriangulation(Tuple{Int,Int}[], NTuple{3,Int}[(1, 2, 3)]), stats
    end
    refl = Int32[]
    @inbounds for i in 1:n
        _is_reflex_at(kernel, vertices, i == 1 ? n : i - 1, i, i == n ? 1 : i + 1) &&
            push!(refl, Int32(i))
    end
    et = _build_edge_tree(vertices, n)
    rpx = Float64[Float64(vertices[v].x) for v in refl]
    rpy = Float64[Float64(vertices[v].y) for v in refl]
    rt = Kd_tree._build_flat_tree!(Kd_tree.FlatNode2D{Float64}, rpx, rpy, copy(refl), n, 8)
    st = _DiagSplitState(vertices, kernel, n, et, rt, Int32[], rng, Tuple{Int,Int}[],
                         NTuple{3,Int}[], stats)
    sizehint!(st.diagonals, n - 3)
    sizehint!(st.triangles, n - 2)
    _process_pieces!(st, Int32[Int32(i) for i in 1:n], refl)
    return PolygonTriangulation(st.diagonals, st.triangles), stats
end

"""
    BasicCompGeometry.triangulate_diagonal_split(vertices; kernel=FloatingKernel(), rng=Random.default_rng()) -> PolygonTriangulation

Triangulate a simple counterclockwise polygon with `n >= 3` vertices by recursively splitting
it along diagonals from reflex vertices, and fanning pieces with at most one reflex vertex.
Each split first tries a random pair of reflex vertices; if that fails, a reflex endpoint is
resolved by ray shooting (first boundary hit, then the angularly first reflex vertex of the
swept triangle); then a halving search; and finally u_kd ear clipping of the piece. All
diagonals are validated by cone tests and a Kd-tree over the edges.
Experimental.

# Preconditions
- `vertices` must form a simple counterclockwise polygon with `n >= 3` vertices.
"""
triangulate_diagonal_split(
    vertices;
    kernel::AbstractKernel = FloatingKernel(),
    rng::AbstractRNG = Random.default_rng(),
) = first(_triangulate_diagonal_split_impl(vertices, kernel, rng))

