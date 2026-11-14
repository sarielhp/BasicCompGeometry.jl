# Systems Code Review Report #003 (Lens: performance)

- **Date**: 2026-11-13
- **Auditor**: Claude / Codex (Tier 1 Standard) via `review-audit`
- **Focus Lens**: `performance`
- **Model Tier**: `STANDARD`
- **Backend**: `auto`
- **Scope**: `src/EarCutKdTreeTriangulation.jl src/EarCutZOrderTriangulation.jl src/VisibilityPolygon.jl src/MinkowskiSum.jl`
- **Status**: Action Required

---

Auditing via Codex CLI (gpt-5.6-sol, tier: standard) [Profile: performance]...
[SEVERITY]: Moderate  
[LOCATION]: `src/VisibilityPolygon.jl:57-65`, with additional boxed captures at `102-112`  
[ROOT CAUSE]: `R` is reassigned and then captured by the comprehension constructing `poly_conv`. Julia lowers it to `Core.Box`, so `poly_conv`, `q_conv`, edge tuples, intersection parameters, and intermediate points become abstractly inferred. The boundary-query `findfirst` closures similarly box branch-local `u_start` and `u_end`. Consequently, the quadratic sector/edge loop at lines 118-151 performs dynamic dispatch and boxes arithmetic results. For concrete `Vector{Point{2,Float64}}` input, `@code_warntype` reports `R::Core.Box`, `poly_conv::Vector`, `q_conv::Any`, and `t1`, `t2`, `pt1`, and `pt2` as `Any`.  
[FAILURE TRACE]: Use a regular 256-vertex `Float64` polygon with `q = point(0.0, 0.0)`. After compilation, one `visibility_polygon(poly, q)` call allocates approximately 24.5 MB. Measurements at 32, 64, 128, and 256 vertices were approximately 0.43, 1.61, 6.21, and 24.49 MB respectively, demonstrating allocation growth inside the sector-by-edge scan rather than only output-sized allocation. Boundary queries additionally exercise the boxed `u_start`/`u_end` closures.  
[REMEDIATION]:
```julia
@inline _construction_type(::ExactKernel, ::Type{T}) where {T} =
    T <: Rational ? T : Rational{BigInt}
@inline _construction_type(::AbstractKernel, ::Type{T}) where {T} = T

@inline function _find_ray(kernel, rays, target)
    @inbounds for i in eachindex(rays)
        compare_angle(kernel, rays[i], target) == 0 && return i
    end
    return nothing
end

T0 = promote_type(eltype(eltype(poly)), eltype(q))
R = _construction_type(kernel, T0)  # assigned exactly once

poly_conv = Point{2,R}[point(R(p.x), R(p.y)) for p in poly]
q_conv = point(R(q.x), R(q.y))

# For vertex/edge cases:
start_ray_idx = _find_ray(kernel, unique_rays, u_start)
end_ray_idx = _find_ray(kernel, unique_rays, u_end)
```

[SEVERITY]: Moderate  
[LOCATION]: `src/MinkowskiSum.jl:191-243`, especially reassignments at `194-195` and `212-216` followed by comprehensions at `207-209` and `232-233`  
[ROOT CAUSE]: `_minkowski_convex` repeatedly reassigns `p1` and `p2`, then captures those variables in comprehensions. Julia therefore represents both as `Core.Box`. This contaminates `n`, `m`, `v1`, `v2`, `edges1`, and `edges2`, and ultimately the advertised linear-time edge merge. With fully concrete `Vector{Point{2,Float64}}` arguments, `@code_warntype` reports `p1` and `p2` as boxed, the edge arrays and lengths as `Any`, and the function body as `Any`.  
[FAILURE TRACE]: Call `_minkowski_convex` through the public `minkowski_sum` API with two ordinary convex `Float64` polygons, such as two 1,000-vertex sampled circles. The instability occurs independently of heterogeneous input or custom kernels: indexing the rotated vertex arrays and comparing merged edges proceeds through values inferred as `Any`. A 1,000-by-1,000 call allocates roughly 560 KB despite an output containing only about 1,000–2,000 static points.  
[REMEDIATION]:
```julia
raw1 = Point{2,R}[point(R(p.x), R(p.y)) for p in poly1]
raw2 = Point{2,R}[point(R(p.x), R(p.y)) for p in poly2]
dedup1 = _dedup_consecutive(raw1, kernel)
dedup2 = _dedup_consecutive(raw2, kernel)

n = length(dedup1)
m = length(dedup2)

# Handle empty/singleton cases using dedup1 and dedup2 here.

p1 = n >= 3 && twice_area(dedup1; kernel=kernel) < 0 ?
     reverse(dedup1) : dedup1
p2 = m >= 3 && twice_area(dedup2; kernel=kernel) < 0 ?
     reverse(dedup2) : dedup2

# p1 and p2 are now assigned once and can safely be captured.
v1 = [p1[mod1(i0 + k - 1, n)] for k in 1:n]
v2 = [p2[mod1(j0 + k - 1, m)] for k in 1:m]
```

[SEVERITY]: Moderate  
[LOCATION]: `src/MinkowskiSum.jl:304-317`  
[ROOT CAUSE]: The point, segment, and polygon convenience overloads explicitly construct `Point{2}[...]`. Because `Point` aliases `SVector`, `Point{2}` is a non-concrete UnionAll type; these expressions create `Vector{SVector{2}}` rather than `Vector{SVector{2,T}}`. Iteration and coordinate access consequently require dynamic dispatch before `_minkowski_convex` can convert the values into a concrete representation.  
[FAILURE TRACE]: For two `Segment{2,Float64}` values, `Point{2}[s.p, s.q]` has element type `SVector{2}` and `isconcretetype(eltype(...)) == false`. After warm-up, `minkowski_sum(s1, s2)` allocated approximately 4,704 bytes, versus approximately 2,176 bytes when the identical endpoints were supplied through `[s.p, s.q]`, which infers `Vector{Point{2,Float64}}`. The same defect affects point–polygon and segment–polygon calls.  
[REMEDIATION]:
```julia
minkowski_sum(s1::Segment{2}, s2::Segment{2}; kernel=FloatingKernel(), kwargs...) =
    minkowski_sum([s1.p, s1.q], [s2.p, s2.q]; kernel=kernel, kwargs...)

minkowski_sum(p::Point{2}, poly::AbstractVector{<:Point{2}}; kernel=FloatingKernel(), kwargs...) =
    minkowski_sum([p], poly; kernel=kernel, kwargs...)

minkowski_sum(s::Segment{2}, poly::AbstractVector{<:Point{2}}; kernel=FloatingKernel(), kwargs...) =
    minkowski_sum([s.p, s.q], poly; kernel=kernel, kwargs...)
```
Apply the symmetric replacements to the reverse-argument overloads.

[SEVERITY]: Moderate  
[LOCATION]: `src/MinkowskiSum.jl:352-379`  
[ROOT CAUSE]: `minkowski_sum_decomposition` reassigns `R` before capturing it in the `p1` and `p2` comprehensions. Julia lowers `R` to `Core.Box`; both converted polygons become abstract `Vector`s, the convex-piece collections retain abstract element types, and the nested piece-pair loop invokes `_minkowski_convex` with an `Any`-inferred result. This affects even the fast path where both inputs are already convex.  
[FAILURE TRACE]: Pass two concrete `Vector{Point{2,Float64}}` convex polygons to `minkowski_sum_decomposition`. No triangulation is needed, but `@code_warntype` still reports `R::Core.Box`, `p1::Vector`, abstract `pieces1`/`pieces2`, and the call at line 379 as returning `Any`. For non-convex inputs, the same instability is repeated across every Cartesian pair of generated triangles.  
[REMEDIATION]:
```julia
T0 = promote_type(_scalar_type(poly1), _scalar_type(poly2))
R = _construction_type(kernel, T0)  # use the single-assignment helper above

p1 = Point{2,R}[point(R(p.x), R(p.y)) for p in poly1]
p2 = Point{2,R}[point(R(p.x), R(p.y)) for p in poly2]

pieces1::Vector{Vector{Point{2,R}}} =
    length(p1) <= 3 || is_convex(p1; kernel=kernel) ?
    [p1] :
    [[p1[i], p1[j], p1[k]] for (i, j, k) in
        triangulate_polygon(p1; kernel=kernel, algorithm=algorithm).triangles]
```
Construct `pieces2` identically. This keeps both the decomposition lists and the nested `_minkowski_convex` calls concretely typed.
