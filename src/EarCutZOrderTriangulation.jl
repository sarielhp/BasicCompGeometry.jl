# Ear-clipping polygon triangulation using a Morton Z-order spatial index over
# active reflex vertices to prune point-in-triangle checks locally.

@inline function _quantize_coord(v::Float64, min_v::Float64, inv_span::Float64)::UInt32
    t = (v - min_v) * inv_span
    if !(t > 0.0)
        return UInt32(0)
    elseif t >= 65535.0
        return UInt32(65535)
    end
    return unsafe_trunc(UInt32, floor(t))
end

@inline function _spread_bits16(x::UInt32)::UInt32
    x = (x | (x << 8)) & 0x00ff00ff
    x = (x | (x << 4)) & 0x0f0f0f0f
    x = (x | (x << 2)) & 0x33333333
    x = (x | (x << 1)) & 0x55555555
    return x
end

@inline _morton_code(ix::UInt32, iy::UInt32)::UInt32 =
    _spread_bits16(ix) | (_spread_bits16(iy) << 1)

function _compute_z_bounds(vertices)
    p1 = vertices[1]
    min_x = Float64(p1.x)
    max_x = min_x
    min_y = Float64(p1.y)
    max_y = min_y
    @inbounds for i in 2:length(vertices)
        p = vertices[i]
        x = Float64(p.x)
        y = Float64(p.y)
        x < min_x && (min_x = x)
        x > max_x && (max_x = x)
        y < min_y && (min_y = y)
        y > max_y && (max_y = y)
    end
    span_x = max_x - min_x
    span_y = max_y - min_y
    inv_x = span_x > 0.0 ? (65535.0 / span_x) : 0.0
    inv_y = span_y > 0.0 ? (65535.0 / span_y) : 0.0
    ok = isfinite(min_x) && isfinite(max_x) && isfinite(min_y) && isfinite(max_y) &&
         isfinite(inv_x) && isfinite(inv_y) && inv_x > 0.0 && inv_y > 0.0
    return ok, min_x, min_y, inv_x, inv_y
end

mutable struct _ZOrderReflexIndex <: _AbstractReflexIndex
    is_reflex::BitVector
    z_code::Vector{UInt32}
    rz_prev::Vector{Int}
    rz_next::Vector{Int}
    rz_head::Int
    rz_tail::Int
    in_rz::BitVector
    b_first::Vector{Int}
    l0::Vector{UInt64}
    l1::UInt64
    shift::Int
    min_x::Float64
    min_y::Float64
    inv_x::Float64
    inv_y::Float64
end

@inline function _z_set_bit!(idx::_ZOrderReflexIndex, b::Int)
    w = ((b - 1) >>> 6) + 1
    idx.l0[w] |= UInt64(1) << ((b - 1) & 63)
    idx.l1 |= UInt64(1) << (w - 1)
    return nothing
end

@inline function _z_clear_bit!(idx::_ZOrderReflexIndex, b::Int)
    w = ((b - 1) >>> 6) + 1
    idx.l0[w] &= ~(UInt64(1) << ((b - 1) & 63))
    if idx.l0[w] == 0
        idx.l1 &= ~(UInt64(1) << (w - 1))
    end
    return nothing
end

@inline function _z_next_bucket(idx::_ZOrderReflexIndex, b_start::Int)::Int
    B = length(idx.b_first)
    b_start > B && return 0
    b_start < 1 && (b_start = 1)
    w = ((b_start - 1) >>> 6) + 1
    off = (b_start - 1) & 63
    masked = idx.l0[w] & (~UInt64(0) << off)
    if masked != 0
        return ((w - 1) << 6) + trailing_zeros(masked) + 1
    end
    w >= length(idx.l0) && return 0
    hi_l1 = idx.l1 & (~UInt64(0) << w)
    hi_l1 == 0 && return 0
    w2 = trailing_zeros(hi_l1) + 1
    return ((w2 - 1) << 6) + trailing_zeros(idx.l0[w2]) + 1
end

@inline function _z_prev_bucket(idx::_ZOrderReflexIndex, b_end::Int)::Int
    b_end < 1 && return 0
    B = length(idx.b_first)
    b_end > B && (b_end = B)
    w = ((b_end - 1) >>> 6) + 1
    off = (b_end - 1) & 63
    masked = idx.l0[w] & ((~UInt64(0)) >>> (63 - off))
    if masked != 0
        return ((w - 1) << 6) + (64 - leading_zeros(masked))
    end
    w <= 1 && return 0
    lo_l1 = idx.l1 & ((~UInt64(0)) >>> (65 - w))
    lo_l1 == 0 && return 0
    w2 = 64 - leading_zeros(lo_l1)
    return ((w2 - 1) << 6) + (64 - leading_zeros(idx.l0[w2]))
end

const _MORTON_X_MASK = 0x55555555
const _MORTON_Y_MASK = 0xaaaaaaaa

@inline function _in_morton_box(z::UInt32, min_zx::UInt32, max_zx::UInt32, min_zy::UInt32, max_zy::UInt32)::Bool
    zx = z & _MORTON_X_MASK
    zy = z & _MORTON_Y_MASK
    return min_zx <= zx <= max_zx && min_zy <= zy <= max_zy
end

# Tropf-Herzog (1981) BIGMIN: smallest Morton code >= z inside 2D box [lo, hi].
function _morton_bigmin(z::UInt32, lo::UInt32, hi::UInt32)::Tuple{Bool,UInt32}
    res = UInt32(0)
    found = false
    for k in 31:-1:0
        axis = isodd(k) ? _MORTON_Y_MASK : _MORTON_X_MASK
        mask = UInt32(1) << k
        low_both = mask - UInt32(1)
        low_axis = axis & low_both
        zb = (z & mask) != 0
        lb = (lo & mask) != 0
        hb = (hi & mask) != 0
        if !zb && !lb && hb
            res = (z & ~low_both) | mask | (lo & (low_both & ~axis))
            found = true
            hi |= low_axis
        elseif !zb && lb && hb
            return true, (z & ~low_both) | mask | (lo & low_both)
        elseif zb && !lb && !hb
            return found, res
        elseif zb && !lb && hb
            lo &= ~low_axis
        end
    end
    return true, z
end

# Tropf-Herzog (1981) LITMAX: largest Morton code <= z inside 2D box [lo, hi].
function _morton_litmax(z::UInt32, lo::UInt32, hi::UInt32)::Tuple{Bool,UInt32}
    res = UInt32(0)
    found = false
    for k in 31:-1:0
        axis = isodd(k) ? _MORTON_Y_MASK : _MORTON_X_MASK
        mask = UInt32(1) << k
        low_both = mask - UInt32(1)
        low_inc = low_both | mask
        low_axis = axis & low_both
        zb = (z & mask) != 0
        lb = (lo & mask) != 0
        hb = (hi & mask) != 0
        if !zb && !lb && hb
            hi |= low_axis
        elseif !zb && lb && hb
            return found, res
        elseif zb && !lb && !hb
            return true, (z & ~low_inc) | (hi & low_both)
        elseif zb && !lb && hb
            res = (z & ~low_inc) | (hi & (low_both & ~axis)) | low_axis
            found = true
            lo &= ~low_axis
        end
    end
    return true, z
end

@inline function _advance_q_in_box(idx::_ZOrderReflexIndex, q::Int, min_z::UInt32, max_z::UInt32)::Int
    zq = idx.z_code[q]
    ok, bm = _morton_bigmin(zq, min_z, max_z)
    (!ok || bm > max_z) && return 0
    bm <= zq && return idx.rz_next[q]
    b_cur = Int(zq >> idx.shift) + 1
    b_dst = Int(bm >> idx.shift) + 1
    if b_dst > b_cur
        nb = _z_next_bucket(idx, b_dst)
        return nb != 0 ? idx.b_first[nb] : 0
    end
    while q != 0 && idx.z_code[q] < bm
        q = idx.rz_next[q]
    end
    return q
end

@inline function _retreat_p_in_box(idx::_ZOrderReflexIndex, p::Int, min_z::UInt32, max_z::UInt32)::Int
    zp = idx.z_code[p]
    ok, lm = _morton_litmax(zp, min_z, max_z)
    (!ok || lm < min_z) && return 0
    lm >= zp && return idx.rz_prev[p]
    b_cur = Int(zp >> idx.shift) + 1
    b_dst = Int(lm >> idx.shift) + 1
    if b_dst < b_cur
        pb = _z_prev_bucket(idx, b_dst)
        pb == 0 && return 0
        nb = _z_next_bucket(idx, pb + 1)
        return nb != 0 ? idx.rz_prev[idx.b_first[nb]] : idx.rz_tail
    end
    while p != 0 && idx.z_code[p] > lm
        p = idx.rz_prev[p]
    end
    return p
end

function _build_zorder_reflex_index(
    vertices,
    is_reflex::BitVector,
    reflex_cnt::Int,
    min_x::Float64,
    min_y::Float64,
    inv_x::Float64,
    inv_y::Float64,
)
    n = length(vertices)
    z_code = Vector{UInt32}(undef, n)
    packed = Vector{UInt64}(undef, reflex_cnt)
    k = 0
    @inbounds for i in 1:n
        p = vertices[i]
        ix = _quantize_coord(Float64(p.x), min_x, inv_x)
        iy = _quantize_coord(Float64(p.y), min_y, inv_y)
        zc = _morton_code(ix, iy)
        z_code[i] = zc
        if is_reflex[i]
            k += 1
            packed[k] = (UInt64(zc) << 32) | UInt64(i)
        end
    end
    sort!(packed)

    B = prevpow(2, clamp(reflex_cnt ÷ 2, 64, 4096))
    shift = 32 - trailing_zeros(B)
    rz_prev = zeros(Int, n)
    rz_next = zeros(Int, n)
    in_rz = falses(n)
    b_first = zeros(Int, B)
    l0 = zeros(UInt64, B >>> 6)
    l1 = UInt64(0)

    rz_head = 0
    rz_tail = 0
    @inbounds for pos in 1:reflex_cnt
        v = Int(packed[pos] & UInt64(0xffffffff))
        in_rz[v] = true
        rz_prev[v] = rz_tail
        if rz_tail != 0
            rz_next[rz_tail] = v
        else
            rz_head = v
        end
        rz_tail = v
        b = Int(z_code[v] >> shift) + 1
        if b_first[b] == 0
            b_first[b] = v
            w = ((b - 1) >>> 6) + 1
            l0[w] |= UInt64(1) << ((b - 1) & 63)
            l1 |= UInt64(1) << (w - 1)
        end
    end

    return _ZOrderReflexIndex(
        is_reflex,
        z_code,
        rz_prev,
        rz_next,
        rz_head,
        rz_tail,
        in_rz,
        b_first,
        l0,
        l1,
        shift,
        min_x,
        min_y,
        inv_x,
        inv_y,
    )
end

@inline function _remove_reflex!(idx::_ZOrderReflexIndex, i::Int)
    idx.is_reflex[i] || return nothing
    idx.is_reflex[i] = false
    idx.in_rz[i] || return nothing
    idx.in_rz[i] = false

    p = idx.rz_prev[i]
    q = idx.rz_next[i]
    if p != 0
        idx.rz_next[p] = q
    else
        idx.rz_head = q
    end
    if q != 0
        idx.rz_prev[q] = p
    else
        idx.rz_tail = p
    end

    b = Int(idx.z_code[i] >> idx.shift) + 1
    if idx.b_first[b] == i
        if q != 0 && (Int(idx.z_code[q] >> idx.shift) + 1) == b
            idx.b_first[b] = q
        else
            idx.b_first[b] = 0
            _z_clear_bit!(idx, b)
        end
    end
    idx.rz_prev[i] = 0
    idx.rz_next[i] = 0
    return nothing
end

function _add_reflex!(idx::_ZOrderReflexIndex, v::Int)
    idx.is_reflex[v] && return nothing
    idx.is_reflex[v] = true
    idx.in_rz[v] && return nothing

    zv = idx.z_code[v]
    b = Int(zv >> idx.shift) + 1
    p = 0
    q = 0
    if idx.b_first[b] != 0
        q = idx.b_first[b]
        while q != 0 && (Int(idx.z_code[q] >> idx.shift) + 1) == b && idx.z_code[q] < zv
            q = idx.rz_next[q]
        end
        p = (q != 0) ? idx.rz_prev[q] : idx.rz_tail
    else
        nb = _z_next_bucket(idx, b + 1)
        if nb != 0
            q = idx.b_first[nb]
            p = idx.rz_prev[q]
        else
            p = idx.rz_tail
            q = 0
        end
    end

    idx.rz_prev[v] = p
    idx.rz_next[v] = q
    if p != 0
        idx.rz_next[p] = v
    else
        idx.rz_head = v
    end
    if q != 0
        idx.rz_prev[q] = v
    else
        idx.rz_tail = v
    end
    idx.in_rz[v] = true

    if idx.b_first[b] == 0 || q == idx.b_first[b]
        idx.b_first[b] = v
        _z_set_bit!(idx, b)
    end
    return nothing
end

@inline function _ear_candidate_blocks(
    poly,
    kernel::AbstractKernel,
    v::Int,
    pr::Int,
    ear::Int,
    nx::Int,
    a::Point,
    b::Point,
    c::Point,
    min_tx,
    max_tx,
    min_ty,
    max_ty,
)::Bool
    (v != pr && v != ear && v != nx) || return false
    pt = poly[v]
    (pt.x >= min_tx && pt.x <= max_tx && pt.y >= min_ty && pt.y <= max_ty) || return false
    return _point_in_ccw_triangle(kernel, pt, a, b, c)
end

@inline function _step_ear_scan(
    idx::_ZOrderReflexIndex,
    poly,
    kernel::AbstractKernel,
    v::Int,
    miss::Int,
    forward::Bool,
    pr::Int,
    ear::Int,
    nx::Int,
    a::Point,
    b::Point,
    c::Point,
    min_tx,
    max_tx,
    min_ty,
    max_ty,
    min_z::UInt32,
    max_z::UInt32,
    min_zx::UInt32,
    max_zx::UInt32,
    min_zy::UInt32,
    max_zy::UInt32,
)::Tuple{Bool,Int,Int}
    zv = idx.z_code[v]
    in_range = forward ? (zv >= min_z) : (zv <= max_z)
    if in_range && _in_morton_box(zv, min_zx, max_zx, min_zy, max_zy)
        _ear_candidate_blocks(poly, kernel, v, pr, ear, nx, a, b, c, min_tx, max_tx, min_ty, max_ty) &&
            return true, 0, 0
        return false, (forward ? idx.rz_next[v] : idx.rz_prev[v]), 0
    end
    miss += 1
    if miss >= 3
        jumped = forward ? _advance_q_in_box(idx, v, min_z, max_z) : _retreat_p_in_box(idx, v, min_z, max_z)
        return false, jumped, 0
    end
    return false, (forward ? idx.rz_next[v] : idx.rz_prev[v]), miss
end

function _ear_blocked(
    idx::_ZOrderReflexIndex,
    poly,
    kernel::AbstractKernel,
    pr::Int,
    ear::Int,
    nx::Int,
    a::Point,
    b::Point,
    c::Point,
)::Bool
    min_tx = min(a.x, b.x, c.x)
    max_tx = max(a.x, b.x, c.x)
    min_ty = min(a.y, b.y, c.y)
    max_ty = max(a.y, b.y, c.y)

    # Pad the discrete Morton box by 1 grid cell so Rational/BigInt conversion rounding
    # cannot exclude an interior reflex vertex on a cell boundary.
    qx0 = _quantize_coord(Float64(min_tx), idx.min_x, idx.inv_x)
    qx1 = _quantize_coord(Float64(max_tx), idx.min_x, idx.inv_x)
    qy0 = _quantize_coord(Float64(min_ty), idx.min_y, idx.inv_y)
    qy1 = _quantize_coord(Float64(max_ty), idx.min_y, idx.inv_y)
    qx0 > UInt32(0) && (qx0 -= UInt32(1))
    qy0 > UInt32(0) && (qy0 -= UInt32(1))
    qx1 < UInt32(65535) && (qx1 += UInt32(1))
    qy1 < UInt32(65535) && (qy1 += UInt32(1))
    min_z = _morton_code(qx0, qy0)
    max_z = _morton_code(qx1, qy1)
    min_zx = min_z & _MORTON_X_MASK
    max_zx = max_z & _MORTON_X_MASK
    min_zy = min_z & _MORTON_Y_MASK
    max_zy = max_z & _MORTON_Y_MASK

    shift = idx.shift
    b_lo = Int(min_z >> shift) + 1
    b_hi = Int(max_z >> shift) + 1
    b_ear = Int(idx.z_code[ear] >> shift) + 1

    b_first_lo = _z_next_bucket(idx, b_lo)
    if b_first_lo == 0 || b_first_lo > b_hi || idx.z_code[idx.b_first[b_first_lo]] > max_z
        return false
    end

    p = 0
    q = 0
    b_mid = _z_next_bucket(idx, b_ear)
    if b_mid != 0 && b_mid <= b_hi && idx.z_code[idx.b_first[b_mid]] <= max_z
        q = idx.b_first[b_mid]
        p = idx.rz_prev[q]
    else
        b_after = _z_next_bucket(idx, b_hi + 1)
        p = b_after != 0 ? idx.rz_prev[idx.b_first[b_after]] : idx.rz_tail
    end

    miss_p = 0
    miss_q = 0
    while p != 0 && idx.z_code[p] >= min_z && q != 0 && idx.z_code[q] <= max_z
        hit, p, miss_p = _step_ear_scan(
            idx, poly, kernel, p, miss_p, false, pr, ear, nx, a, b, c,
            min_tx, max_tx, min_ty, max_ty, min_z, max_z, min_zx, max_zx, min_zy, max_zy,
        )
        hit && return true
        hit, q, miss_q = _step_ear_scan(
            idx, poly, kernel, q, miss_q, true, pr, ear, nx, a, b, c,
            min_tx, max_tx, min_ty, max_ty, min_z, max_z, min_zx, max_zx, min_zy, max_zy,
        )
        hit && return true
    end
    while p != 0 && idx.z_code[p] >= min_z
        hit, p, miss_p = _step_ear_scan(
            idx, poly, kernel, p, miss_p, false, pr, ear, nx, a, b, c,
            min_tx, max_tx, min_ty, max_ty, min_z, max_z, min_zx, max_zx, min_zy, max_zy,
        )
        hit && return true
    end
    while q != 0 && idx.z_code[q] <= max_z
        hit, q, miss_q = _step_ear_scan(
            idx, poly, kernel, q, miss_q, true, pr, ear, nx, a, b, c,
            min_tx, max_tx, min_ty, max_ty, min_z, max_z, min_zx, max_zx, min_zy, max_zy,
        )
        hit && return true
    end
    return false
end

"""
    triangulate_ear_cut_z_order(vertices; kernel=FloatingKernel()) -> PolygonTriangulation

Triangulate a simple counterclockwise polygon with `n >= 3` vertices using Z-order ear clipping
with a Morton spatial index over active reflex vertices.

For inputs with at most 16 reflex vertices (including convex polygons) or non-finite
`Float64` bounding spans, falls back automatically to the direct reflex-list ear clipper.
Worst-case running time is `O(n r)`, where `r` is the reflex vertex count. Candidate
reflex vertices are pruned only by the ear's bounding box on a fixed 16-bit Morton grid over
the polygon's bounding box, so it is fast when ears stay small and well spread, but long
thin ears (e.g. random star polygons) or a few far-away vertices make it degrade toward
quadratic; prefer `triangulate_ear_cut_kd_tree` or `triangulate_sweep` for such inputs.

# Preconditions
- `vertices` must form a simple counterclockwise polygon with `n >= 3` vertices.
"""
function triangulate_ear_cut_z_order(
    vertices;
    kernel::AbstractKernel = FloatingKernel(),
)
    n = _validate_ccw_polygon(kernel, vertices)
    return _ear_clip_with_index(vertices, n, kernel) do is_reflex, reflex_cnt
        reflex_cnt <= 16 && return _build_plain_reflex_index(is_reflex)
        ok, min_x, min_y, inv_x, inv_y = _compute_z_bounds(vertices)
        ok || return _build_plain_reflex_index(is_reflex)
        return _build_zorder_reflex_index(vertices, is_reflex, reflex_cnt, min_x, min_y, inv_x, inv_y)
    end
end

export triangulate_ear_cut_z_order
