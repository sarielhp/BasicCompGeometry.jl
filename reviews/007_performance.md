# Systems Code Review Report #007 (Lens: performance)

- **Date**: 2026-11-13
- **Auditor**: Flagship Opus / Deep (Tier 2 Deep) via `review-audit`
- **Focus Lens**: `performance`
- **Model Tier**: `DEEP`
- **Backend**: `auto`
- **Scope**: `src/Points.jl src/Kernels.jl src/KernelPredicates.jl src/Segments.jl src/Planes.jl src/BBoxes.jl src/Transforms2D.jl src/VirtArray.jl`
- **Status**: Action Required

---

Auditing via Codex CLI (gpt-5.6-sol, tier: deep) [Profile: performance]...
[SEVERITY]: Moderate  
[LOCATION]: `src/Points.jl`, `dist_subspace(p, q, dir)` return expression  
[ROOT CAUSE]: The positive branch returns `sqrt(d_sq)`, while the nonpositive branch returns `zero(T)`. For integer or rational coordinates these have different types, so inference produces a value-dependent union in a predicate likely to be called repeatedly.  
[FAILURE TRACE]: With `Point{2,Int}`, a nonzero projected distance returns `Float64`, while identical points return `Int64(0)`. A loop accumulating `dist_subspace` results therefore requires union splitting or boxing in less-specialized consumers.  
[REMEDIATION]:
```julia
function dist_subspace(p::Point{D,T}, q::Point{D,T}, dir::Point{D,T}) where {D,T}
    diff = p - q
    proj_len = dot(diff, dir)
    d_sq = dot(diff, diff) - proj_len^2
    return sqrt(max(zero(d_sq), d_sq))
end
```

[SEVERITY]: Moderate  
[LOCATION]: `src/Points.jl`, `rand_gaussian(D::Int, ...)` and `rand_point(D::Int, ...)`  
[ROOT CAUSE]: The runtime value `D` is inserted into the result type `Point{D,T}`. Julia cannot generally infer a static-array type parameter from a nonconstant `Int`, so callers receive an abstract `SVector` result despite every individual invocation constructing a concrete point.  
[FAILURE TRACE]: Reading `D` from configuration and executing `p = rand_point(D); dist_sq(p, p)` leaves the dimension of `p` uninferred. Repeating this in a simulation introduces dynamic dispatch and prevents unrolling of static-vector operations.  
[REMEDIATION]:
```julia
rand_gaussian(::Val{D}, ::Type{T} = Float64) where {D,T<:Real} =
    Point{D,T}(ntuple(_ -> randn(T), Val(D)))

rand_point(::Val{D}, ::Type{T} = Float64) where {D,T<:Real} =
    rand(Point{D,T})

# Hot callers use rand_point(Val(3)); retain the Int overload only as
# explicitly dynamic convenience API if required.
```

[SEVERITY]: Moderate  
[LOCATION]: `src/Segments.jl`, `nn_point(s_p, s_q, qr)`  
[ROOT CAUSE]: The degenerate branch returns the original `Point{D,T}`, whereas the projection branch can promote coordinates through division and the `0.0`/`1.0` clamp bounds. The return type consequently depends on whether the segment happens to have zero length.  
[FAILURE TRACE]: For integer endpoints, `nn_point(point(0,0), point(0,0), point(1,1))` returns `Point{2,Int}`, while a nondegenerate integer segment returns a floating-point point. A workload containing both ordinary and zero-length segments forces union-typed nearest points in its inner loop. Float32 inputs can likewise be widened by the Float64 clamp constants.  
[REMEDIATION]:
```julia
function nn_point(
    s_p::Point{D,T},
    s_q::Point{D,T},
    qr::Point{D,T},
) where {D,T}
    R = typeof(one(T) / one(T))
    p = Point{D,R}(s_p)
    q = Point{D,R}(s_q)
    r = Point{D,R}(qr)

    delta = q - p
    d_sq = dot(delta, delta)
    iszero(d_sq) && return p

    t = dot(r - p, delta) / d_sq
    return convex_comb(p, q, clamp(t, zero(t), one(t)))
end
```

[SEVERITY]: Moderate  
[LOCATION]: `src/Segments.jl`, `bisection_point(seg, p, q)` parallel and nonparallel returns  
[ROOT CAUSE]: The parallel branch hard-codes `0.0` and returns `seg.p`; the nonparallel branch returns the arithmetic type of `t` and `at(seg,t)`. These tuple fields differ for integer, rational, and Float32 geometry, making the result dependent on the parallelism test.  
[FAILURE TRACE]: An integer segment parallel to the bisector returns `Tuple{Bool,Float64,Point{D,Int}}`; a nonparallel integer segment returns a floating-point intersection point. Processing both cases in one clipping or Voronoi loop prevents a single concrete tuple type.  
[REMEDIATION]:
```julia
function bisection_point(
    seg::Segment{D,T},
    p::Point{D,T},
    q::Point{D,T},
) where {D,T}
    R = typeof(one(T) / one(T))
    sp, sq = Point{D,R}(seg.p), Point{D,R}(seg.q)
    pr, qr = Point{D,R}(p), Point{D,R}(q)

    dir = qr - pr
    mid = (qr + pr) / R(2)
    pos = dot(dir, mid)
    denom = dot(dir, sq - sp)

    if abs(denom) < R(1e-12)
        return false, zero(R), sp
    end

    t = (pos - dot(dir, sp)) / denom
    return zero(R) <= t <= one(R), t, convex_comb(sp, sq, t)
end
```

[SEVERITY]: Moderate  
[LOCATION]: `src/Segments.jl`, entire `dist_segment_segment(...)` block constructing `c`, `m`, and `sol`  
[ROOT CAUSE]: Every distance query constructs a heap `Vector`, a heap `Matrix`, invokes generic `det` and `norm`, and runs a generic linear solve for a fixed 2×2 system. This creates avoidable allocation and linear-algebra dispatch overhead in a fundamental pairwise-distance operation.  
[FAILURE TRACE]: Computing distances between thousands of segment pairs—such as collision checking or nearest-neighbor filtering—allocates at least the RHS, matrix, and solve result for every nonparallel pair. The problem size is statically two, so none of these heap objects is algorithmically necessary.  
[REMEDIATION]:
```julia
A = dot(v_2, v_2)
B = dot(v_2, v_3)
C = dot(v_3, v_3)
r1 = -dot(v_1, v_2)
r2 = -dot(v_1, v_3)
Δ = A * C - B * B

scale = max(abs(A * C), abs(B * B), one(float(Δ)))
if abs(Δ) <= 1e-12 * scale
    # Existing endpoint fallback
else
    s = clamp((r1 * C - B * r2) / Δ, 0.0, 1.0)
    t = clamp((A * r2 - B * r1) / Δ, 0.0, 1.0)
end
```

[SEVERITY]: Moderate  
[LOCATION]: `src/Transforms2D.jl`, `apply_transform(m, p)` augmented-vector multiplication  
[ROOT CAUSE]: `[p[1], p[2], one(T)]` creates a heap vector, and generic matrix-vector multiplication creates another result vector merely to compute two scalar coordinates.  
[FAILURE TRACE]: Transforming a point cloud with `apply_transform.(Ref(m), points)` performs at least two temporary vector allocations per point. For one million points, this creates millions of short-lived heap objects despite the fixed 3×3 operation.  
[REMEDIATION]:
```julia
@inline function apply_transform(
    m::AbstractMatrix{T},
    p::Point{2,T},
)::Point{2,T} where {T}
    x = muladd(m[1, 1], p.x, muladd(m[1, 2], p.y, m[1, 3]))
    y = muladd(m[2, 1], p.x, muladd(m[2, 2], p.y, m[2, 3]))
    return Point{2,T}(x, y)
end
```

[SEVERITY]: Moderate  
[LOCATION]: `src/BBoxes.jl`, `union_bbox(objects::AbstractVector) = union_bbox(objects...)`  
[ROOT CAUSE]: A runtime vector is expanded into a vararg call. Large vectors therefore require construction and dispatch of a correspondingly large argument tuple, while the vararg implementation also materializes `boxes = map(bbox, objects)`. This adds O(n) temporary storage and can cause excessive specialization/dispatch overhead for what should be a streaming reduction.  
[FAILURE TRACE]: Calling `union_bbox` on a vector of tens of thousands of boxes or segments splats every element into a separate positional argument before reduction. Memory and compilation/dispatch costs grow with the number of objects, independently of the constant-sized result.  
[REMEDIATION]:
```julia
function union_bbox(objects::AbstractVector)
    isempty(objects) &&
        throw(ArgumentError("union_bbox requires at least one object"))

    first_box = bbox(first(objects))
    D = length(first_box.mini)
    T = mapreduce(o -> eltype(bbox(o).mini), promote_type, objects)
    result = BBox{D,T}()

    for object in objects
        bb = bbox(object)
        length(bb.mini) == D ||
            throw(DimensionMismatch("all objects must have the same dimension"))
        bb.f_init || continue
        bound!(result, Point{D,T}(bb.mini))
        bound!(result, Point{D,T}(bb.maxi))
    end
    return result
end
```
