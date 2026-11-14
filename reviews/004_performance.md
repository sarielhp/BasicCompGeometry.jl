# Systems Code Review Report #004 (Lens: performance)

- **Date**: 2026-11-13
- **Auditor**: Gemini 3.8 Flash (Tier 0 Workhorse) via `review-audit`
- **Focus Lens**: `performance`
- **Model Tier**: `TIER0`
- **Backend**: `gemini`
- **Scope**: `src/EarCutKdTreeTriangulation.jl src/EarCutZOrderTriangulation.jl src/VisibilityPolygon.jl src/MinkowskiSum.jl`
- **Status**: Action Required

---

Auditing via Gemini Flash (bws run) [Profile: performance]...
[SEVERITY]: Major
[LOCATION]: [VisibilityPolygon.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/VisibilityPolygon.jl#L57-L173) (`visibility_polygon`, lines 57–65 and 98–173)
[ROOT CAUSE]: `R` is conditionally reassigned (`R = Rational{BigInt}` at line 61) and captured by the array comprehension `[point(R(p.x), R(p.y)) for p in poly]` at line 64, which lowers to a closure over `R`. Julia therefore wraps `R` in a `Core.Box` (`Any`), destroying type inference for `poly_conv` (`Vector`), `q_conv` (`Any`), `directions`, `unique_rays`, `edges`, `raw_polygon`, `deduped`, and `final_vertices`, and forcing runtime dynamic dispatch on every iteration of the `O(m * n)` ray-edge intersection loop (lines 118–173). Additionally, `u_start` and `u_end` are assigned in two `elseif` branches (lines 102–103 and 109–110) and captured by `findfirst` closures (lines 104–105, 111–112), creating two further `Core.Box` allocations.
[FAILURE TRACE]: Calling `visibility_polygon(pts_star, point(0.0, 0.0))` on a `Vector{Point{2, Float64}}` with `FloatingKernel()` yields `@code_warntype` slots `R@_41::Core.Box`, `u_start::Core.Box`, `u_end::Core.Box`, `poly_conv::Vector`, `q_conv::Any`, `best_t1::Any`, `t1::Any`, and return type `Vector{SVector{2, _A}} where _A`. On a 40-vertex star polygon it allocates 656,936 bytes; on a 200-vertex star polygon it allocates 15,027,824 bytes.
[REMEDIATION]:
```julia
@inline _vis_coord_type(::FloatingKernel, ::Type{T1}, ::Type{T2}) where {T1, T2} =
    float(promote_type(T1, T2))
@inline _vis_coord_type(::ExactKernel, ::Type{T1}, ::Type{T2}) where {T1, T2} =
    promote_type(T1, T2) <: Rational ? promote_type(T1, T2) : Rational{BigInt}

# Inside visibility_polygon, bind R once and pass to a typed helper / function barrier:
R = _vis_coord_type(kernel, eltype(eltype(poly)), eltype(q))
return _visibility_polygon_impl(poly, q, kernel, R)

# Inside _visibility_polygon_impl, replace boxed u_start/u_end captures with a helper:
@inline _find_ray(kernel, unique_rays, u) =
    findfirst(r -> compare_angle(kernel, r, u) == 0, unique_rays)::Int
```

[SEVERITY]: Major
[LOCATION]: [VisibilityPolygon.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/VisibilityPolygon.jl#L301-L395) (`polygon_kernel`, lines 301–306 and 322–395)
[ROOT CAUSE]: `R` is reassigned at line 304 (`R = Rational{BigInt}`) and captured by `[point(R(p.x), R(p.y)) for p in poly]` at line 306, boxing `R` in a `Core.Box`. This turns `poly_conv` into an untyped `Vector`, `curr_poly`/`next_poly`/`deduped`/`final_poly` into `Vector{SVector{2, _A}} where _A`, and all geometric predicates/intersections inside the `for i in 1:n; for j in 1:m` Sutherland–Hodgman clipping loop into dynamically dispatched `Any` operations. Furthermore, each outer iteration `1:n` allocates two fresh vectors (`next_poly` at line 327 and `deduped` at line 377) rather than reusing scratch buffers.
[FAILURE TRACE]: `@code_warntype` on `polygon_kernel(pts_star)` (`Vector{Point{2, Float64}}`, `FloatingKernel()`) reports `R@_21::Core.Box`, `poly_conv::Vector`, `ab::Any`, `d::Any`, `pt::Any`, and return type `Vector{SVector{2, _A}} where _A`. Measured allocations are 392,792 bytes for 40 vertices and 7,018,312 bytes for 200 vertices.
[REMEDIATION]:
```julia
R = _vis_coord_type(kernel, eltype(eltype(poly)), eltype(eltype(poly)))
return _polygon_kernel_impl(poly, kernel, R)

# Inside _polygon_kernel_impl(..., ::Type{R}) where {R}:
poly_conv = Point{2, R}[point(R(p.x), R(p.y)) for p in poly]
curr_poly = Point{2, R}[point(min_x, min_y), point(max_x, min_y), point(max_x, max_y), point(min_x, max_y)]
next_poly = Point{2, R}[]
deduped = Point{2, R}[]
for i in 1:n
    isempty(curr_poly) && break
    empty!(next_poly)
    # ... clip into next_poly ...
    empty!(deduped)
    # ... deduplicate next_poly into deduped ...
    curr_poly, deduped = deduped, curr_poly
end
```

[SEVERITY]: Major
[LOCATION]: [VisibilityPolygon.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/VisibilityPolygon.jl#L124-L127) and [VisibilityPolygon.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/VisibilityPolygon.jl#L454-L492) (`visibility_polygon` and `_ray_edge_intersect`)
[ROOT CAUSE]: `_ray_edge_intersect` has a value-dependent return type whenever `R` is `Float32` or `Int` under `FloatingKernel()`: miss branches (lines 459, 467, 491) return `(false, zero(R), q)::Tuple{Bool, R, Point{2, R}}`, while the hit branch (lines 483–488) computes `s_clamped = clamp(s, 0.0, 1.0)` using `Float64` literals and `t = num_t / denom`, returning `Tuple{Bool, Float64, Point{2, Float64}}` (for `R = Int`) or `Tuple{Bool, Float32, Point{2, Float64}}` (for `R = Float32`). In `visibility_polygon` lines 124–127, `best_t1` is initialized as `typemax(Rational{BigInt})` or `Inf` (`Float64`) regardless of `R`, introducing `Union{Rational{BigInt}, Rational{Int64}}` when `R = Rational{Int64}` with `ExactKernel()`, and causing `push!(raw_polygon, best_pt1)` to fail with `InexactError` when `R = Int` with `FloatingKernel()`.
[FAILURE TRACE]:
1. `@code_warntype BasicCompGeometry._ray_edge_intersect(FloatingKernel(), q32, u32, a32, b32, Float32)` returns `Union{Tuple{Bool, Float32, SVector{2, Float32}}, Tuple{Bool, Float32, SVector{2, Float64}}}`.
2. `visibility_polygon([point(0, 0), point(10, 0), point(10, 10), point(5, 3), point(0, 10)], point(2, 1))` throws `InexactError(:Int64, (Int64, 4.455445544554456))` (and `polygon_kernel` on the same integer polygon throws `InexactError(:Int64, (Int64, 2.8571428571428568))`).
[REMEDIATION]:
```julia
# Promote FloatingKernel coordinates to float(R) in visibility_polygon/polygon_kernel,
# initialize best_t1 = kernel isa ExactKernel ? typemax(R) : R(Inf),
# and keep Float32/Float64 arithmetic homogeneous in _ray_edge_intersect:
s_clamped = clamp(s, zero(t), one(t))
pt = point(R(a.x + s_clamped * ba_x), R(a.y + s_clamped * ba_y))
return (true, R(t), pt)
```

[SEVERITY]: Major
[LOCATION]: [MinkowskiSum.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/MinkowskiSum.jl#L181-L263) and [MinkowskiSum.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/MinkowskiSum.jl#L346-L383) (`_minkowski_convex` and `minkowski_sum_decomposition`)
[ROOT CAUSE]:
1. In `_minkowski_convex` (lines 191–236), `p1` and `p2` are reassigned (`p1 = _dedup_consecutive(p1, kernel)` at line 194 and `p1 = reverse(p1)` at line 213; likewise `p2` at lines 195 and 216) and then captured in the comprehensions `v1 = [p1[mod1(i0 + k - 1, n)] for k in 1:n]` and `v2 = [p2[mod1(j0 + k - 1, m)] for k in 1:m]`. Julia boxes `p1` and `p2` in `Core.Box`, turning `v1`, `v2`, `edges1`, `edges2`, `e1`, `e2`, and the return type of `_minkowski_convex` into `Any`. Every step of the `O(n + m)` merge loop (lines 242–263) runs via dynamic dispatch.
2. In `minkowski_sum_decomposition` (lines 352–358), `R` is conditionally reassigned (`R = Rational{BigInt}`) and captured by `[point(R(p.x), R(p.y)) for p in poly1]` / `poly2`, boxing `R` in `Core.Box`. `p1`/`p2` become untyped `Vector`, `pieces1`/`pieces2` and `results` hold abstract `Vector{SVector{2, _A}} where _A` elements, and the `O(n * m)` double loop (lines 377–381) invokes `_minkowski_convex` via dynamic dispatch on every triangle pair.
3. Each call to `_minkowski_convex` in that `O(n * m)` loop allocates up to 11 temporary `Vector`s (`p1`, `p2`, two `_dedup_consecutive` outputs, `v1`, `v2`, `edges1`, `edges2`, `result`, `deduped`, `final_poly`), even when `c1` and `c2` are already typed, CCW, 3-vertex triangles.
[FAILURE TRACE]: `@code_warntype` on `_minkowski_convex(sq, sq, FloatingKernel())` shows `p1@_25::Core.Box`, `p2@_24::Core.Box`, `v1::Any`, `edges1::Any`, `Body::Any`; `@code_warntype` on `minkowski_sum_decomposition(pts_star, pts_star)` shows `R@_17::Core.Box`, `p1::Vector`, `Body::Vector{Vector{SVector{2, _A}}} where _A`, and dynamically dispatched `_minkowski_convex(...)::Any`. Calling `minkowski_sum_decomposition(pts_star, pts_star)` on two 40-vertex star polygons allocates 5,105,480 bytes (`minkowski_sum(sq, sq)` on 4-vertex squares allocates 10,784 bytes).
[REMEDIATION]:
```julia
# 1. Eliminate Core.Box by separating typed conversion/normalization from comprehensions:
R = _minkowski_coord_type(kernel, _scalar_type(poly1), _scalar_type(poly2))
return _minkowski_sum_decomposition_impl(poly1, poly2, kernel, algorithm, R)

# 2. In _minkowski_convex, avoid reassigning variables captured by comprehensions (or use loops),
#    and index into p1/p2 directly with modular offsets i0/j0 instead of allocating v1, v2, edges1, edges2:
p1_ccw = (n >= 3 && twice_area(p1_dedup; kernel = kernel) < 0) ? reverse(p1_dedup) : p1_dedup
p2_ccw = (m >= 3 && twice_area(p2_dedup; kernel = kernel) < 0) ? reverse(p2_dedup) : p2_dedup
result = Vector{Point{2, R}}()
sizehint!(result, n + m)
```

[SEVERITY]: Moderate
[LOCATION]: [MinkowskiSum.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/MinkowskiSum.jl#L57-L115) (`is_convex`)
[ROOT CAUSE]: `is_convex` triggers three `Core.Box` captures: `R` (reassigned at line 59 and captured at line 61), `deduped` (reassigned at line 73 and captured at line 104), and `idx_min` (mutated at line 100 and captured at line 104). Consequently, `p_conv`, `edges`, `nz_edges`, `area2`, and `ts` are untyped (`Vector`/`Any`) and `twice_area`, `orient_sign`, and `compare_angle` dispatch dynamically. Additionally, lines 104–105 allocate two vectors (`edges` and `nz_edges`) solely to verify consecutive angle monotonicity in a single scan.
[FAILURE TRACE]: `@code_warntype` on `is_convex(sq)` (`Vector{Point{2, Float64}}`) exhibits `idx_min@_13::Core.Box`, `deduped@_15::Core.Box`, `R@_19::Core.Box`, `p_conv::Vector`, `edges::Any`, and `nz_edges::Any`.
[REMEDIATION]:
```julia
# Compute R before a typed helper _is_convex_impl(poly, kernel, R), bind deduped_ccw without
# reassigning deduped, and check edge angle monotonicity in-place without allocating edges/nz_edges:
prev_edge = Point{2, R}(zero(R), zero(R))
has_prev = false
nz_count = 0
for idx in 1:k
    e = deduped_ccw[mod1(idx_min + idx, k)] - deduped_ccw[mod1(idx_min + idx - 1, k)]
    (iszero(e.x) && iszero(e.y)) && continue
    nz_count += 1
    if has_prev && compare_angle(kernel, prev_edge, e) > 0
        return false
    end
    prev_edge = e
    has_prev = true
end
return nz_count >= 3
```

[SEVERITY]: Moderate
[LOCATION]: [MinkowskiSum.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/MinkowskiSum.jl#L304-L317) (`minkowski_sum` overloads for `Point{2}` and `Segment{2}`)
[ROOT CAUSE]: The convenience overloads wrap points and segments into `Point{2}[...]` vectors. Because `Point{2}` (`SVector{2, T} where T`) is an abstract `UnionAll` type rather than a concrete type, `Point{2}[...]` constructs an abstractly typed `Vector{Point{2}}` whose elements are heap-boxed, causing dynamic dispatch in downstream `is_convex` and `_minkowski_convex` calls.
[FAILURE TRACE]: `typeof(Point{2}[s1.p, s1.q])` evaluates to `Vector{SVector{2}}` (`isconcretetype(eltype(...)) == false`), and `minkowski_sum(s1, s2)` for two 2-point `Segment{2, Float64}` values allocates 8,832 bytes.
[REMEDIATION]:
```julia
@inline _pts_vec(s::Segment{2, T}) where {T} = Point{2, T}[s.p, s.q]
@inline _pts_vec(p::Point{2, T}) where {T} = Point{2, T}[p]

minkowski_sum(s1::Segment{2}, s2::Segment{2}; kernel::AbstractKernel = FloatingKernel(), kwargs...) =
    minkowski_sum(_pts_vec(s1), _pts_vec(s2); kernel = kernel, kwargs...)
minkowski_sum(p::Point{2}, poly::AbstractVector{<:Point{2}}; kernel::AbstractKernel = FloatingKernel(), kwargs...) =
    minkowski_sum(_pts_vec(p), poly; kernel = kernel, kwargs...)
# ... and symmetrically for the remaining Point/Segment overloads.
```
