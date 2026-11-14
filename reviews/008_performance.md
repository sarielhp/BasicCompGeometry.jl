# Systems Code Review Report #008 (Lens: performance)

- **Date**: 2026-11-14
- **Auditor**: Gemini 3.8 Flash (Tier 0 Workhorse) via `review-audit`
- **Focus Lens**: `performance`
- **Model Tier**: `TIER0`
- **Backend**: `gemini`
- **Scope**: `src/PntSeqs.jl src/PntSeqHierarchy.jl src/PntSeqHausdorff.jl src/ChainIndex.jl src/LocalIncidence.jl src/Bundles.jl src/Arrangements.jl`
- **Status**: Action Required

---

Auditing via Gemini Flash (bws run) [Profile: performance]...
[SEVERITY]: Major
[LOCATION]: [Arrangements.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/Arrangements.jl#L446-L475) (`vertical_decomposition`)
[ROOT CAUSE]: Line 446 initializes `new_active = _ActiveTrap[]` using the unparameterized `UnionAll` type `_ActiveTrap` instead of `_ActiveTrap{R_type}`, and line 474 reassigns `active_traps = new_active`. Inside `vertical_decomposition`, `@code_warntype` confirms `active_traps` becomes `Union{Vector{_ActiveTrap{R_type}}, Vector{_ActiveTrap}}`, `new_active` is `Vector{_ActiveTrap}`, and loop elements `at` are typed abstractly as `_ActiveTrap`. Consequently, from slab `k = 2` onward, every call to `_vd_find_active`, `_vd_can_extend` (`at.bottom == bot_b && at.top == top_b`), `at.x_right = x_r`, and `_vd_make_trap` performs dynamic dispatch and boxes struct fields in the `O(n^3)`-iteration sweep.
[FAILURE TRACE]: On `160` non-parallel `Line{2,Float64}` inputs (`12,654` vertices, `38,244` trapezoids), `vertical_decomposition(ls, bb)` performs `72.10 M` allocations (`2.795 GiB`, `0.617 s`). Measuring `_vd_find_active` directly on `80` active trapezoids allocates `12,912` bytes per call with `Vector{_ActiveTrap}` versus `112` bytes per call with `Vector{_ActiveTrap{Float64}}`.
[REMEDIATION]:
```julia
# Preallocate two concrete buffers outside the slab loop and swap/empty! them:
active_traps = _ActiveTrap{R_type}[]
new_active = _ActiveTrap{R_type}[]
# Inside the slab loop:
empty!(new_active)
...
active_traps, new_active = new_active, active_traps
```

[SEVERITY]: Major
[LOCATION]: [Arrangements.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/Arrangements.jl#L514-L533) (`_tc_link_neighbors!`, called by `trapezoid_chains` at line 554)
[ROOT CAUSE]: `_tc_link_neighbors!` links adjacent trapezoids using an all-pairs nested loop `for i = 1:n, j = 1:n` over `n = length(traps)` trapezoids, resulting in `O(n^2)` time (`O(k^4)` for `k` input lines since `n = O(k^2)`). It also repeats `any(xv -> abs(xv - ti.x_right) <= tol, vert_line_x)` inside the inner `j` loop even though `ti.x_right` depends only on `i`.
[FAILURE TRACE]: On `80` lines (`9,592` trapezoids), `trapezoid_chains(decomp)` runs in `0.037 s`; doubling to `160` lines (`38,244` trapezoids, `4x` traps) increases runtime to `1.060 s` (`28.6x` slower), consistent with `O(n^2) = O(k^4)` scaling.
[REMEDIATION]:
```julia
function _tc_link_neighbors!(right_neighbor, left_neighbor, traps, vert_line_x, tol)
    n = length(traps)
    order = sortperm(traps, by = t -> t.x_left)
    x_lefts = [traps[idx].x_left for idx in order]
    sort!(vert_line_x)
    for i = 1:n
        ti = traps[i]
        xr = ti.x_right
        if !isempty(vert_line_x)
            k = searchsortedfirst(vert_line_x, xr - tol)
            (k <= length(vert_line_x) && abs(vert_line_x[k] - xr) <= tol) && continue
        end
        lo = searchsortedfirst(x_lefts, xr - tol)
        hi = searchsortedlast(x_lefts, xr + tol)
        for pos = lo:hi
            j = order[pos]
            i == j && continue
            tj = traps[j]
            y_bot = max(ti.corners[2][2], tj.corners[1][2])
            y_top = min(ti.corners[3][2], tj.corners[4][2])
            if y_top - y_bot > tol
                right_neighbor[i] = j
                left_neighbor[j] = i
            end
        end
    end
end
```

[SEVERITY]: Major
[LOCATION]: [Arrangements.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/Arrangements.jl#L290-L325) (`_vd_walls`, called by `vertical_decomposition` at line 430)
[ROOT CAUSE]: Inside `for v in in_view_verts` (`O(k^2)` vertices for `k` lines): (1) `findfirst(x -> abs(x - v[1]) <= tol, unique_x)` linearly scans the already-sorted `unique_x` vector (`O(k^2)` per vertex, `O(k^4)` total; similarly at line 320 for vertical lines); (2) `non_vert = filter(...)`, `y_vals = [...]`, `above_y = filter(...)`, and `below_y = filter(...)` allocate four `O(k)` temporary vectors on every vertex (`O(k^3)` allocations).
[FAILURE TRACE]: On `160` lines (`12,654` in-view vertices), `_vd_walls` alone executes `163.83 k` allocations totaling `130.032 MiB`.
[REMEDIATION]:
```julia
function _nearest_x_idx(unique_x, x_target, tol)
    idx = searchsortedfirst(unique_x, x_target)
    if idx <= length(unique_x) && abs(unique_x[idx] - x_target) <= tol
        return idx
    elseif idx > 1 && abs(unique_x[idx-1] - x_target) <= tol
        return idx - 1
    end
    return nothing
end

# Inside _vd_walls (pre-filter non_vert once, or check l.u[1] inline):
for v in in_view_verts
    k_idx = _nearest_x_idx(unique_x, v[1], tol)
    k_idx === nothing && continue
    x_val = unique_x[k_idx]
    wall_top = y_max
    wall_bot = y_min
    vy = v[2]
    for l in lines
        abs(l.u[1]) > 1e-14 || continue
        y = _line_y_at_x(l, x_val)
        if y > vy + tol
            wall_top = min(wall_top, y)
        elseif y < vy - tol
            wall_bot = max(wall_bot, y)
        end
    end
    push!(walls_at_event[k_idx], _VertWall{R}(wall_bot, wall_top))
end
```

[SEVERITY]: Major
[LOCATION]: [Bundles.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/Bundles.jl#L75-L85) (`_unique_points`, called by `vertices(...; unique=true)` and `src/Arrangements.jl:226` `_vd_in_view_vertices`)
[ROOT CAUSE]: `_unique_points` checks every incoming point `p` via a linear scan `any(q -> dist_sq(p, q) <= tol_sq, res)` over all previously kept points `res`, incurring `O(m^2)` distance checks for `m` points (`O(k^4)` when `m = O(k^2)` arrangement vertices).
[FAILURE TRACE]: Timing `_unique_points` on `5,000` vs `10,000` distinct `Point2F` points increases runtime from `3.88 ms` to `15.38 ms` (`4x` for `2x` points), directly bottlenecking `vertices(lines; unique=true)` and `vertical_decomposition` as `k` grows.
[REMEDIATION]:
```julia
function _unique_points(pts::Vector{Point{2,Float64}}; tol::Real = 1e-12)
    isempty(pts) && return pts
    tol_sq = tol^2
    order = sortperm(pts, by = p -> (p[1], p[2]))
    res = Point{2,Float64}[]
    for idx in order
        p = pts[idx]
        dup = false
        for j in length(res):-1:1
            q = res[j]
            (p[1] - q[1]) > tol && break
            if dist_sq(p, q) <= tol_sq
                dup = true
                break
            end
        end
        dup || push!(res, p)
    end
    return res
end
```

[SEVERITY]: Moderate
[LOCATION]: [Arrangements.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/Arrangements.jl#L351-L374), [Arrangements.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/Arrangements.jl#L390-L398), [Arrangements.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/Arrangements.jl#L436-L470) (`_vd_slab_lines`, `_vd_find_active`, `vertical_decomposition`)
[ROOT CAUSE]: In each of the `O(k^2)` slabs (`for k = 1:(length(unique_x)-1)`), `vertical_decomposition` calls `_vd_slab_lines` (allocating `slab_lines` and `unique_slab` and sorting `O(k)` lines) and allocates `bounds`, then for each of the `O(k)` `bounds` pairs linearly scans `active_traps` via `_vd_find_active` and shifts the array via `deleteat!(active_traps, matched_idx)` (`O(k^2)` per slab, `O(k^4)` total).
[FAILURE TRACE]: For `k` lines in general position, `length(unique_x)` is `O(k^2)` (`3,225` events for `80` lines), so every slab repeats `O(k)`-sized vector allocations and `O(k^2)` linear search + `deleteat!` operations.
[REMEDIATION]: Reuse preallocated scratch buffers (`slab_lines`, `unique_slab`, `bounds`) across slabs via `empty!`, and since both `bounds` and `active_traps` are ordered bottom-to-top by `y`, match them with a linear two-pointer scan instead of `_vd_find_active` + `deleteat!`.

[SEVERITY]: Moderate
[LOCATION]: [PntSeqs.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/PntSeqs.jl#L153-L177), [PntSeqs.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/PntSeqs.jl#L266-L273) (`at`, `sample_uniformly`)
[ROOT CAUSE]: `sample_uniformly(P, n)` calls `at(P, (i-1)/(n-1))` in a loop `for i = 1:n`, and `at` calls `prefix_lengths(pnt_seq)` (line 158) on every invocation, allocating an `O(m)` `Vector{Float64}` and summing all `m = length(P)` edges `n` times (`O(n * m)` time and allocations instead of `O(m + n)`). In addition, `at` is type-unstable when `T != Float64` because boundary branches return `pnt_seq.pnts[...]` (`Point{D,T}`) while interpolation returns `convex_comb(..., local_t::Float64)` (`Point{D,Float64}`).
[FAILURE TRACE]: Calling `sample_uniformly(p, 1000)` on a `1,000`-point `PntSeq2F` allocates `8,105,624` bytes (~`8.1 MB` for a `16 KB` output). `@code_warntype at(p32, 0.5f0)` on `PntSeq{2,Float32}` infers `Body::Union{SVector{2,Float32}, SVector{2,Float64}}`.
[REMEDIATION]:
```julia
function sample_uniformly(P::AbsPntSeq{D,T}, n::Int) where {D,T}
    n <= 1 && throw(ArgumentError("n must be > 1"))
    m = length(P)
    m == 0 && throw(ArgumentError("Empty point sequence"))
    OutT = promote_type(T, Float64)
    out = Vector{Point{D,OutT}}(undef, n)
    if m == 1
        fill!(out, Point{D,OutT}(P[1]))
        return PntSeq(out)
    end
    plens = prefix_lengths(P)
    total_len = plens[end]
    if total_len == 0.0
        fill!(out, Point{D,OutT}(P[1]))
        return PntSeq(out)
    end
    seg = 2
    for i = 1:n
        target = ((i - 1) / (n - 1)) * total_len
        while seg < m && plens[seg] < target
            seg += 1
        end
        denom = plens[seg] - plens[seg - 1]
        local_t = denom == 0.0 ? 0.0 : (target - plens[seg - 1]) / denom
        out[i] = Point{D,OutT}(convex_comb(P[seg - 1], P[seg], local_t))
    end
    return PntSeq(out)
end
```

[SEVERITY]: Moderate
[LOCATION]: [PntSeqHierarchy.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/PntSeqHierarchy.jl#L135-L147), [PntSeqHierarchy.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/PntSeqHierarchy.jl#L174-L187), [PntSeqHierarchy.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/PntSeqHierarchy.jl#L322-L330) (`SimplificationHierarchy`, `_fill!`, `_p_extract!`)
[ROOT CAUSE]: `SimplificationHierarchy{D,T,M<:HierarchyMeasure}` declares field `plt::Palette{T}` without the second type parameter `R<:SplitRule` of `Palette{T,R}`, making `plt` an abstract `UnionAll` field (`isconcretetype == false`) and forcing dynamic dispatch on `palette_extract(h.P, h.plt, ...)` in `hierarchy_approx` (line 461) and `hierarchy_build!` (line 493). Furthermore, `_fill!` (lines 136, 138) and `_p_extract!` (line 175) return the `Float64` literal `0.0` in base cases instead of `zero(T)`, introducing `Union{T, Float64}` type instability across recursive calls when `T != Float64`.
[FAILURE TRACE]: For `h = hierarchy_init(pts)` on `PntSeq2F`, `fieldtype(typeof(h), :plt)` is `Palette{Float64}` (`isconcretetype=false`), and `@code_warntype hierarchy_approx(h, 0.5, 25, 4.0)` shows `Base.getproperty(h, :plt)::Palette{Float64}` followed by a dynamically dispatched call to `palette_extract`. For `PntSeq{2,Float32}`, `@code_warntype _fill!` and `_p_extract!` infer `Body::Union{Tuple{Float32, Int64}, Tuple{Float64, Int64}}` and `Body::Union{Float32, Float64}`.
[REMEDIATION]:
```julia
mutable struct SimplificationHierarchy{D,T,R<:SplitRule,M<:HierarchyMeasure}
    len::T
    P::PntSeq{D,T}
    plt::Palette{T,R}
    polys::Vector{PntSeq{D,T}}
    widths::Vector{Float64}
    lock::ReentrantLock
    measure::M
end

# In _fill! (lines 136, 138): return (zero(eltype(errs)), 0)
# In _p_extract! (line 175): k == 0 && return zero(eltype(pal))
```

[SEVERITY]: Moderate
[LOCATION]: [PntSeqs.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/PntSeqs.jl#L61-L66) (`dist(::AbsPntSeq, i::Int, j::Int)`)
[ROOT CAUSE]: `dist(P::AbsPntSeq, i::Int, j::Int)` returns the `Float64` literal `0.0` when `i == j`, while `dist(P.pnts[i], P.pnts[j])` returns `float(T)` (e.g. `Float32` or `BigFloat`).
[FAILURE TRACE]: `@code_warntype dist(p32, 1, 2)` on `PntSeq{2,Float32}` infers `Body::Union{Float32, Float64}`.
[REMEDIATION]:
```julia
@inline function dist(P::AbsPntSeq{D,T}, i::Int, j::Int) where {D,T}
    if i == j
        return zero(float(T))
    end
    return dist(P.pnts[i], P.pnts[j])
end
```

[SEVERITY]: Moderate
[LOCATION]: [LocalIncidence.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/LocalIncidence.jl#L40-L55), [LocalIncidence.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/LocalIncidence.jl#L79-L93) (`vertical_ray_segment`, `first_positive_vertical_contact`)
[ROOT CAUSE]: `vertical_ray_segment` unconditionally widens all six input coordinates `(p, a, b)` to `BigInt` / `Rational{BigInt}` via `_exact` (lines 42-43) before checking the `x`-span condition `(px < ax || px > bx)` (line 50). During a linear scan in `first_positive_vertical_contact`, segments whose `x`-range does not contain `_cx(p)` still heap-allocate `BigInt` temporaries.
[FAILURE TRACE]: Calling `first_positive_vertical_contact(segs, (500, -1), 1)` on `1,000` integer segments allocates `346,288` bytes (~`346` bytes per segment) even though `999` segments do not span `x = 500`.
[REMEDIATION]:
```julia
function vertical_ray_segment(p, direction::Integer, a, b)
    s = _check_direction(direction)
    px_raw, ax_raw, bx_raw = _cx(p), _cx(a), _cx(b)
    if ax_raw == bx_raw
        px_raw == ax_raw || return NoHit()
        return _vertical_overlap(p, s, _exact(px_raw), _exact(_cy(p)), _exact(ax_raw), _exact(_cy(a)), _exact(_cy(b)))
    end
    lo_x, hi_x = ax_raw < bx_raw ? (ax_raw, bx_raw) : (bx_raw, ax_raw)
    (px_raw < lo_x || px_raw > hi_x) && return NoHit()
    px, py = _exact(px_raw), _exact(_cy(p))
    ax, ay, bx, by = _exact(ax_raw), _exact(_cy(a)), _exact(bx_raw), _exact(_cy(b))
    if ax > bx
        ax, ay, bx, by = bx, by, ax, ay
    end
    y = ay + (px - ax) * (by - ay) // (bx - ax)
    t = s * (y - py)
    t < 0 && return NoHit()
    return PointHit(_mkpoint(p, Rational{BigInt}(px), Rational{BigInt}(y)), Rational{BigInt}(t))
end
```
