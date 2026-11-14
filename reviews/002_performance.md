# Systems Code Review Report #002 (Lens: performance)

- **Date**: 2026-11-13
- **Auditor**: Gemini 3.8 Flash (Tier 0 Workhorse) via `review-audit`
- **Focus Lens**: `performance`
- **Model Tier**: `TIER0`
- **Backend**: `gemini`
- **Scope**: `src/EarCutKdTreeTriangulation.jl src/EarCutZOrderTriangulation.jl src/VisibilityPolygon.jl src/MinkowskiSum.jl`
- **Status**: Action Required

---

Auditing via Gemini Flash (bws run) [Profile: performance]...
[SEVERITY]: Moderate
[LOCATION]: src/VisibilityPolygon.jl:57-65, 98-113, 118-173
[ROOT CAUSE]: Variable capture and reassignment in `visibility_polygon` forces compiler variable boxing (`Core.Box`), leading to severe container-level type instability and allocation churn in the inner visibility sweep. Specifically:
1. `R` is reassigned (`R = Rational{BigInt}`) and captured by comprehension `poly_conv = [point(R(p.x), R(p.y)) for p in poly]`, placing `R` into a heap-allocated `Core.Box`. Consequently, `directions`, `unique_rays`, `raw_polygon`, `deduped`, and `final_vertices` are instantiated with non-concrete element types (`Array{Point{2, _A}} where _A`).
2. `u_start` and `u_end` are assigned across multiple branches and captured by closures `r -> compare_angle(kernel, r, u_start) == 0`, placing both vectors into `Core.Box`.
3. `best_t1` and `best_t2` are initialized to `Inf` (Float64) regardless of `R`, mutating types across edge iterations.

In the inner ray-sweep loop (`for (e_idx, (a, b)) in enumerate(edges)`), calls to `_ray_edge_intersect` dynamically dispatch twice per edge, allocating intermediate tuples and boxing `best_t1`/`best_t2`. For a 100-vertex polygon, 100 visibility queries churn 11.51 million heap allocations and 365 MiB of garbage.
[FAILURE TRACE]:
```julia
poly = [point(cos(t), sin(t)) for t in range(0, 2pi, length=101)[1:100]]
q = point(0.0, 0.0)
# 100 calls produce 11.51M allocations and 365 MiB allocation churn:
@time for _ in 1:100; visibility_polygon(poly, q); end
# @code_warntype reveals:
#   R@_41::CORE.BOX
#   u_start::CORE.BOX
#   u_end::CORE.BOX
#   Body::ARRAY{STATICARRAYSCORE.SVECTOR{2, _A}, 1} WHERE _A
```
[REMEDIATION]:
Avoid variable reassignment before closure capture, initialize ray distance sentinels with `typemax(R)`, and preserve concrete element types:
```julia
function visibility_polygon(
    poly::AbstractVector{Point{2, T_poly}},
    q::Point{2, T_q};
    kernel::AbstractKernel = FloatingKernel(),
) where {T_poly, T_q}
    n = length(poly)
    n >= 3 || throw(ArgumentError("polygon must have at least 3 vertices; got n = $n"))
    area2 = twice_area(poly; kernel = kernel)
    area2 > 0 || throw(ArgumentError("polygon must be oriented counterclockwise with positive area; got twice_area = $area2"))

    R_init = promote_type(T_poly, T_q)
    R = (kernel isa ExactKernel && !(R_init <: Rational)) ? Rational{BigInt} : R_init

    poly_conv = Point{2, R}[point(R(p.x), R(p.y)) for p in poly]
    q_conv = point(R(q.x), R(q.y))

    loc_type, loc_idx = _point_in_polygon(kernel, q_conv, poly_conv)
    loc_type != :exterior || throw(ArgumentError("vantage point q = $q must lie inside or on the boundary of the polygon"))

    directions = Point{2, R}[]
    sizehint!(directions, n)
    for p in poly_conv
        if p != q_conv
            push!(directions, p - q_conv)
        end
    end

    sort!(directions, lt = (d1, d2) -> angle_less(kernel, d1, d2))

    unique_rays = Point{2, R}[]
    sizehint!(unique_rays, length(directions))
    for d in directions
        if isempty(unique_rays) || compare_angle(kernel, d, last(unique_rays)) != 0
            push!(unique_rays, d)
        end
    end
    m = length(unique_rays)
    m >= 2 || throw(ArgumentError("polygon vertices must span at least 2 distinct directions from q"))

    edges = [(poly_conv[j], poly_conv[mod1(j + 1, n)]) for j in 1:n]

    start_ray_idx = 1
    end_ray_idx = 1
    if loc_type == :vertex
        cu = poly_conv[loc_idx]
        nx = poly_conv[mod1(loc_idx + 1, n)]
        pr = poly_conv[mod1(loc_idx - 1, n)]
        u_s = nx - cu
        u_e = pr - cu
        s_idx = findfirst(r -> compare_angle(kernel, r, u_s) == 0, unique_rays)
        e_idx = findfirst(r -> compare_angle(kernel, r, u_e) == 0, unique_rays)
        start_ray_idx = s_idx === nothing ? 1 : s_idx
        end_ray_idx = e_idx === nothing ? 1 : e_idx
    elseif loc_type == :edge
        a = poly_conv[loc_idx]
        b = poly_conv[mod1(loc_idx + 1, n)]
        u_s = b - q_conv
        u_e = a - q_conv
        s_idx = findfirst(r -> compare_angle(kernel, r, u_s) == 0, unique_rays)
        e_idx = findfirst(r -> compare_angle(kernel, r, u_e) == 0, unique_rays)
        start_ray_idx = s_idx === nothing ? 1 : s_idx
        end_ray_idx = e_idx === nothing ? 1 : e_idx
    end

    raw_polygon = Point{2, R}[]
    sizehint!(raw_polygon, 2m)
    curr_ray_idx = start_ray_idx

    inf_val = typemax(R)
    while true
        next_ray_idx = mod1(curr_ray_idx + 1, m)
        u1 = unique_rays[curr_ray_idx]
        u2 = unique_rays[next_ray_idx]

        best_edge_idx = 0
        best_t1 = inf_val
        best_t2 = inf_val
        best_pt1 = q_conv
        best_pt2 = q_conv

        for (e_idx, (a, b)) in enumerate(edges)
            hit1, t1, pt1 = _ray_edge_intersect(kernel, q_conv, u1, a, b, R)
            hit2, t2, pt2 = _ray_edge_intersect(kernel, q_conv, u2, a, b, R)

            if hit1 && hit2
                is_closer = false
                if best_edge_idx == 0
                    is_closer = true
                else
                    diff1 = t1 - best_t1
                    diff2 = t2 - best_t2
                    if kernel isa ExactKernel
                        if (diff1 < 0 && diff2 <= 0) || (diff1 <= 0 && diff2 < 0)
                            is_closer = true
                        end
                    else
                        if (diff1 < -1e-9 && diff2 < 1e-9) || (diff2 < -1e-9 && diff1 < 1e-9)
                            is_closer = true
                        end
                    end
                end
                if is_closer
                    best_edge_idx = e_idx
                    best_t1 = t1
                    best_t2 = t2
                    best_pt1 = pt1
                    best_pt2 = pt2
                end
            end
        end

        if best_edge_idx != 0
            push!(raw_polygon, best_pt1)
            push!(raw_polygon, best_pt2)
        end

        if loc_type != :interior && next_ray_idx == end_ray_idx
            break
        end
        curr_ray_idx = next_ray_idx
        if loc_type == :interior && curr_ray_idx == start_ray_idx
            break
        end
    end
    # ... remaining deduplication with Point{2, R}[]
```

---

[SEVERITY]: Moderate
[LOCATION]: src/VisibilityPolygon.jl:301-395
[ROOT CAUSE]: In `polygon_kernel`, `R` is reassigned (`R = Rational{BigInt}`) and captured by comprehension `poly_conv = [point(R(p.x), R(p.y)) for p in poly]`, forcing `R` into a `Core.Box`. As a consequence, `curr_poly`, `next_poly`, `deduped`, and `final_poly` are typed as `Array{Point{2, _A}} where _A` rather than concrete `Vector{Point{2, R}}`. 

In the Sutherland-Hodgman clipping loops (`for i in 1:n` and `for j in 1:m`), point accesses (`p_prev`, `p_curr`), edge arithmetic, and clipping intersections dynamically unbox and dynamically dispatch. Calling `polygon_kernel` 1,000 times on a 100-vertex polygon triggers 61.35 million allocations (1.776 GB).
[FAILURE TRACE]:
```julia
poly = [point(cos(t), sin(t)) for t in range(0, 2pi, length=101)[1:100]]
# 1,000 calls produce 61.35M allocations (1.776 GiB, 1.326s):
@time for _ in 1:1000; BasicCompGeometry.polygon_kernel(poly; kernel=FloatingKernel()); end
```
[REMEDIATION]:
Compute coordinate type `R` via single assignment without mutation prior to the generator/comprehension, and type-annotate intermediate polygon buffers:
```julia
function polygon_kernel(
    poly::AbstractVector{Point{2, T}};
    kernel::AbstractKernel = FloatingKernel(),
) where {T}
    n = length(poly)
    n >= 3 || throw(ArgumentError("polygon must have at least 3 vertices; got n = $n"))
    area2 = twice_area(poly; kernel = kernel)
    area2 > 0 || throw(ArgumentError("polygon must be oriented counterclockwise with positive area; got twice_area = $area2"))

    R = (kernel isa ExactKernel && !(T <: Rational)) ? Rational{BigInt} : T
    poly_conv = Point{2, R}[point(R(p.x), R(p.y)) for p in poly]

    min_x = minimum(p.x for p in poly_conv)
    max_x = maximum(p.x for p in poly_conv)
    min_y = minimum(p.y for p in poly_conv)
    max_y = maximum(p.y for p in poly_conv)

    curr_poly = Point{2, R}[
        point(min_x, min_y),
        point(max_x, min_y),
        point(max_x, max_y),
        point(min_x, max_y),
    ]

    tol = 1e-11

    for i in 1:n
        a = poly_conv[i]
        b = poly_conv[mod1(i + 1, n)]

        isempty(curr_poly) && break
        next_poly = Point{2, R}[]
        m = length(curr_poly)
        # ... Sutherland-Hodgman clipping preserving Point{2, R} ...
```

---

[SEVERITY]: Moderate
[LOCATION]: src/MinkowskiSum.jl:50-116
[ROOT CAUSE]: In `is_convex`, multiple variable reassignments prior to generator/comprehension closures cause Julia's lowering compiler to box `R`, `deduped`, and `idx_min` into `Core.Box`:
1. `R` is reassigned (`R = Rational{BigInt}`) and captured by `p_conv = [point(R(p.x), R(p.y)) for p in poly]`.
2. `deduped` is reassigned (`deduped = reverse(deduped)`).
3. `idx_min` is updated in a loop (`idx_min = i`).
4. Both `deduped` and `idx_min` are captured by `edges = [deduped[mod1(idx_min + idx, k)] - ... for idx in 1:k]`.

Because `deduped` and `idx_min` reside in `Core.Box`, vertex indexing `deduped[i]` in the orientation check loop (`for i in 1:k`) and minimum vertex search dynamically unbox on every iteration, dispatching dynamically to `orient_sign`, `dist`, and `compare_angle`. 1,000 calls on a 100-vertex polygon generate 2.31 million heap allocations (69.5 MiB).
[FAILURE TRACE]:
```julia
poly = [point(cos(t), sin(t)) for t in range(0, 2pi, length=101)[1:100]]
# 1,000 calls produce 2.31M allocations (69.5 MiB):
@time for _ in 1:1000; is_convex(poly); end
# @code_warntype on var"#is_convex#254" displays:
#   R@_19::CORE.BOX
#   deduped@_15::CORE.BOX
#   idx_min@_13::CORE.BOX
#   edges::ANY
#   nz_edges::ANY
```
[REMEDIATION]:
Assign `R` and `deduped` without variable reassignment, and compute `edges` using a typed comprehension `Point{2, R}[...]`:
```julia
function is_convex(
    poly::AbstractVector{Point{2, T}};
    kernel::AbstractKernel = FloatingKernel(),
) where {T}
    n = length(poly)
    n <= 2 && return true

    R = (kernel isa ExactKernel && !(T <: Rational)) ? Rational{BigInt} : T
    p_conv = Point{2, R}[point(R(p.x), R(p.y)) for p in poly]
    deduped_raw = _dedup_consecutive(p_conv, kernel)

    k = length(deduped_raw)
    k <= 2 && return true

    area2 = twice_area(deduped_raw; kernel = kernel)
    iszero(area2) && return false

    deduped = area2 < 0 ? reverse(deduped_raw) : deduped_raw

    tol = 1e-11
    for i in 1:k
        p = deduped[mod1(i - 1, k)]
        q = deduped[i]
        r = deduped[mod1(i + 1, k)]

        ts = orient_sign(kernel, p, q, r)
        if kernel isa ExactKernel
            ts < 0 && return false
        else
            if ts < 0
                d1 = dist(p, q)
                d2 = dist(q, r)
                cross = (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x)
                if cross < -tol * max(d1 * d2, 1.0)
                    return false
                end
            end
        end
    end

    idx_min = 1
    for i in 2:k
        if deduped[i].y < deduped[idx_min].y || (deduped[i].y == deduped[idx_min].y && deduped[i].x < deduped[idx_min].x)
            idx_min = i
        end
    end

    edges = Point{2, R}[deduped[mod1(idx_min + idx, k)] - deduped[mod1(idx_min + idx - 1, k)] for idx in 1:k]
    nz_edges = filter(e -> !(iszero(e.x) && iszero(e.y)), edges)
    length(nz_edges) < 3 && return false

    for idx in 1:(length(nz_edges) - 1)
        cmp = compare_angle(kernel, nz_edges[idx], nz_edges[idx + 1])
        if cmp > 0
            return false
        end
    end

    return true
end
```

---

[SEVERITY]: Moderate
[LOCATION]: src/MinkowskiSum.jl:212-263
[ROOT CAUSE]: In `_minkowski_convex`, an inner helper `function find_start(pts)` is defined directly inside the method body. Because `p1` and `p2` are reassigned (`p1 = reverse(p1)` and `p2 = reverse(p2)`), the closure capture mechanism of Julia lowers the local activation frame by boxing `p1` and `p2` into `Core.Box`.

Consequently:
1. `v1`, `v2`, `edges1`, and `edges2` are inferred as `Vector{Any}`.
2. In the cyclic two-pointer angular merge loop (`while i <= n || j <= m`), every vector element read (`v1[mod1(i, n)]`, `edges1[i]`), addition (`+`), and `compare_angle(kernel, e1, e2)` is dynamically dispatched on `Any` and boxed.
3. 1,000 invocations churn 762,000 allocations (46.18 MiB).
[FAILURE TRACE]:
```julia
poly1 = [point(cos(t), sin(t)) for t in range(0, 2pi, length=101)[1:100]]
poly2 = [point(cos(t), sin(t)) for t in range(0, 2pi, length=101)[1:100]]
# 1,000 calls produce 762k allocations (46.18 MiB, 25-36% GC time):
@time for _ in 1:1000; BasicCompGeometry._minkowski_convex(poly1, poly2, FloatingKernel()); end
# @code_warntype reveals:
#   p1@_25::CORE.BOX
#   p2@_24::CORE.BOX
#   edges1::ANY
#   edges2::ANY
#   v1::ANY
#   v2::ANY
#   cmp::Int64 (via dynamic dispatch in every merge loop step)
```
[REMEDIATION]:
Lift `_find_start_vertex` to a file-level `@inline` function, avoid mutating `p1`/`p2`, and construct typed vectors:
```julia
@inline function _find_start_vertex(pts)
    idx = 1
    for i in 2:length(pts)
        if pts[i].y < pts[idx].y || (pts[i].y == pts[idx].y && pts[i].x < pts[idx].x)
            idx = i
        end
    end
    return idx
end

function _minkowski_convex(
    poly1::AbstractVector{Point{2, T1}},
    poly2::AbstractVector{Point{2, T2}},
    kernel::AbstractKernel,
) where {T1, T2}
    R = promote_type(T1, T2)
    R_exact = (kernel isa ExactKernel && !(R <: Rational)) ? Rational{BigInt} : R

    p1_raw = Point{2, R_exact}[point(R_exact(p.x), R_exact(p.y)) for p in poly1]
    p2_raw = Point{2, R_exact}[point(R_exact(p.x), R_exact(p.y)) for p in poly2]

    p1_d = _dedup_consecutive(p1_raw, kernel)
    p2_d = _dedup_consecutive(p2_raw, kernel)

    n = length(p1_d)
    m = length(p2_d)

    if n == 0 || m == 0
        return Point{2, R_exact}[]
    elseif n == 1 && m == 1
        return Point{2, R_exact}[p1_d[1] + p2_d[1]]
    elseif n == 1
        return Point{2, R_exact}[q + p1_d[1] for q in p2_d]
    elseif m == 1
        return Point{2, R_exact}[p + p2_d[1] for p in p1_d]
    end

    p1 = (n >= 3 && twice_area(p1_d; kernel = kernel) < 0) ? reverse(p1_d) : p1_d
    p2 = (m >= 3 && twice_area(p2_d; kernel = kernel) < 0) ? reverse(p2_d) : p2_d

    i0 = _find_start_vertex(p1)
    j0 = _find_start_vertex(p2)

    v1 = Point{2, R_exact}[p1[mod1(i0 + k - 1, n)] for k in 1:n]
    v2 = Point{2, R_exact}[p2[mod1(j0 + k - 1, m)] for k in 1:m]

    edges1 = Point{2, R_exact}[v1[mod1(k + 1, n)] - v1[k] for k in 1:n]
    edges2 = Point{2, R_exact}[v2[mod1(k + 1, m)] - v2[k] for k in 1:m]

    result = Point{2, R_exact}[]
    sizehint!(result, n + m)
    # ... linear merge loop now executes with zero allocations and static dispatch ...
```

---

[SEVERITY]: Moderate
[LOCATION]: src/MinkowskiSum.jl:304-318
[ROOT CAUSE]: Specialized overloads for `Point` and `Segment` construct intermediate point vectors using `Point{2}[...]` without coordinate type specialization:
```julia
minkowski_sum(s1::Segment{2}, s2::Segment{2}; ...) =
    minkowski_sum(Point{2}[s1.p, s1.q], Point{2}[s2.p, s2.q]; ...)

minkowski_sum(p::Point{2}, poly::AbstractVector{<:Point{2}}; ...) =
    minkowski_sum(Point{2}[p], poly; ...)
```
Because `Point{2}` is an unspecialized type alias (`Point{2, T} where T`), `Point{2}[p]` creates `Vector{StaticArraysCore.SVector{2}}` whose element type is a `UnionAll`. This instantiates an abstractly typed container (`isconcretetype(eltype(v)) == false`), destroying type inference at the call site and returning `Body::ANY`.
[FAILURE TRACE]:
```julia
p = point(1.0, 2.0)
poly = [point(0.0, 0.0), point(1.0, 0.0), point(0.0, 1.0)]
# @code_warntype on minkowski_sum(p, poly) reports:
#   Body::ANY
```
[REMEDIATION]:
Parameterize the signatures and construct concretely typed vectors:
```julia
minkowski_sum(s1::Segment{2, T1}, s2::Segment{2, T2}; kernel::AbstractKernel = FloatingKernel(), kwargs...) where {T1, T2} =
    minkowski_sum(Point{2, T1}[s1.p, s1.q], Point{2, T2}[s2.p, s2.q]; kernel = kernel, kwargs...)

minkowski_sum(p::Point{2, T1}, poly::AbstractVector{Point{2, T2}}; kernel::AbstractKernel = FloatingKernel(), kwargs...) where {T1, T2} =
    minkowski_sum(Point{2, T1}[p], poly; kernel = kernel, kwargs...)

minkowski_sum(poly::AbstractVector{Point{2, T1}}, p::Point{2, T2}; kernel::AbstractKernel = FloatingKernel(), kwargs...) where {T1, T2} =
    minkowski_sum(poly, Point{2, T2}[p]; kernel = kernel, kwargs...)

minkowski_sum(s::Segment{2, T1}, poly::AbstractVector{Point{2, T2}}; kernel::AbstractKernel = FloatingKernel(), kwargs...) where {T1, T2} =
    minkowski_sum(Point{2, T1}[s.p, s.q], poly; kernel = kernel, kwargs...)

minkowski_sum(poly::AbstractVector{Point{2, T1}}, s::Segment{2, T2}; kernel::AbstractKernel = FloatingKernel(), kwargs...) where {T1, T2} =
    minkowski_sum(poly, Point{2, T2}[s.p, s.q]; kernel = kernel, kwargs...)
```

---

[SEVERITY]: Moderate
[LOCATION]: src/EarCutZOrderTriangulation.jl:186-196, 203-215
[ROOT CAUSE]: In `_advance_q_in_box` and `_retreat_p_in_box`, when target Morton code `bm` (or `lm`) falls into a different bucket than `b_cur`, the routine fetches `q = idx.b_first[nb]` (or the tail vertex of bucket `pb`) and returns it immediately:
```julia
    if b_dst > b_cur
        nb = _z_next_bucket(idx, b_dst)
        return nb != 0 ? idx.b_first[nb] : 0
    end
    while q != 0 && idx.z_code[q] < bm
        q = idx.rz_next[q]
    end
    return q
```
When `nb == b_dst`, `idx.b_first[nb]` is the first vertex in bucket `b_dst`, which often has `z_code < bm` (since `bm` can lie in the middle of bucket `b_dst`). Because it returns without advancing `q` to `>= bm`, `_step_ear_scan` encounters candidate vertices with Morton codes strictly outside the active query box, registering up to 3 spurious misses and executing extra bounding-box and in-triangle queries before a subsequent jump executes the inner `while` scan.
[FAILURE TRACE]:
Ear-clipping triangulation on polygons where reflex vertices span multiple buckets and `bm` falls strictly inside bucket `b_dst` beyond the bucket's first vertex.
[REMEDIATION]:
Unify the bucket jump and internal scan so the helper always returns a candidate satisfying `idx.z_code[q] >= bm` (and `idx.z_code[p] <= lm`):
```julia
@inline function _advance_q_in_box(idx::_ZOrderReflexIndex, q::Int, min_z::UInt32, max_z::UInt32)::Int
    zq = idx.z_code[q]
    ok, bm = _morton_bigmin(zq, min_z, max_z)
    (!ok || bm > max_z) && return 0
    bm <= zq && return idx.rz_next[q]
    b_cur = Int(zq >> idx.shift) + 1
    b_dst = Int(bm >> idx.shift) + 1
    if b_dst > b_cur
        nb = _z_next_bucket(idx, b_dst)
        nb == 0 && return 0
        q = idx.b_first[nb]
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
        p = nb != 0 ? idx.rz_prev[idx.b_first[nb]] : idx.rz_tail
    end
    while p != 0 && idx.z_code[p] > lm
        p = idx.rz_prev[p]
    end
    return p
end
```
