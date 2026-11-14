# Systems Code Review Report #001 (Lens: performance)

- **Date**: 2026-11-13
- **Auditor**: Gemini 3.8 Flash (Tier 0 Workhorse) via `review-audit`
- **Focus Lens**: `performance`
- **Model Tier**: `TIER0`
- **Backend**: `gemini`
- **Scope**: `src/MonotoneTriangulation.jl src/SweepTriangulation.jl src/SeidelTriangulation.jl`
- **Status**: Action Required

---

Auditing via Gemini Flash (bws run) [Profile: performance]...
Successfully registered persistent workspace: /home/sariel/prog/26/BasicCompGeometry
[SEVERITY]: Moderate
[LOCATION]: src/MonotoneTriangulation.jl:173-175
[ROOT CAUSE]: In `triangulate_monotone`, when the sweep processes a vertex `v` belonging to the opposite chain from `last(stack)`, the stack is emptied by popping all elements. At line 174, instead of reusing the existing allocated buffer `stack`, a new heap vector is allocated via `stack = Int[order[j-1], v]`. This discards the existing buffer and forces repeated allocations and resizings during every chain alternation in the inner triangulation loop.
[FAILURE TRACE]: Triangulating any x-monotone polygon whose vertices alternate frequently between upper and lower boundary chains (such as a comb, staircase, or zigzag polygon with $k$ vertices):
```julia
# A polygon with k vertices alternating between upper (y=1) and lower (y=-1) chains
z = [point(0.0, 0.0); [point(2i - 1.0, 1.0) for i in 1:1000]; point(2001.0, 0.0); [point(2i, -1.0) for i in 1000:-1:1]]
triangulate_monotone(z) # Allocates hundreds of small Vector{Int} instances inside the inner loop
```
[REMEDIATION]:
```julia
            pop!(stack)
            empty!(stack)
            push!(stack, order[j-1], v)
```

---

[SEVERITY]: Moderate
[LOCATION]: src/MonotoneTriangulation.jl:155-161, 178-185
[ROOT CAUSE]: Redundant geometric predicate evaluation in the inner triangulation loop. At lines 179-181:
```julia
turn = orient(kernel, vertices[a], vertices[q], vertices[v])
(is_upper[v] ? turn < 0 : turn > 0) || break
push!(triangles, make_ccw(a, q, v))
```
The turn orientation `turn` is already computed by `orient(...)`. When `is_upper[v]` is true, `turn < 0`, so `(a, q, v)` forms a clockwise turn and CCW order is known to be `(a, v, q)`. When `is_upper[v]` is false, `turn > 0` and CCW order is `(a, q, v)`. Calling `make_ccw(a, q, v)` immediately re-evaluates `orient_sign(kernel, vertices[a], vertices[q], vertices[v])` on the exact same vertex triplet, doubling expensive kernel predicate evaluations (especially noticeable under `ExactKernel()`). Furthermore, for cross-chain triangles (lines 170 and 193), the CCW orientation is deterministically known from `is_upper[v]` without any predicate call.
[FAILURE TRACE]:
```julia
poly = [point(0, 0), point(1, 2), point(2, 3), point(3, 1), point(4, 0), point(2, -1)]
# Under ExactKernel(), every triangle formed within the chain evaluates orient_sign twice per triangle
triangulate_monotone(poly; kernel = ExactKernel())
```
[REMEDIATION]:
```julia
    # Orientation follows directly from chain classification without re-evaluating predicates
    @inline emit_tri!(i1::Int, i2::Int, i3::Int) = push!(triangles, (i1, i2, i3))

    # Cross-chain:
    is_upper[v] ? emit_tri!(v, s_next, s_top) : emit_tri!(v, s_top, s_next)

    # Same-chain:
    is_upper[v] ? emit_tri!(a, v, q) : emit_tri!(a, q, v)

    # Terminal:
    is_upper[s_top] ? emit_tri!(v_last, s_top, s_next) : emit_tri!(v_last, s_next, s_top)
```

---

[SEVERITY]: Moderate
[LOCATION]: src/MonotoneTriangulation.jl:421-432
[ROOT CAUSE]: `_check_distinct_vertices` allocates an intermediate `Vector{Tuple{T, T}}` via a comprehension `[(vertices[i][1], vertices[i][2]) for i in 1:n]` solely to pass to `sortperm`. For large polygons validated in `_validate_ccw_polygon` (called by `triangulate_ear_cut` and `triangulate_seidel`), this creates an unnecessary $O(n)$ array allocation and garbage collection overhead.
[FAILURE TRACE]:
```julia
poly = [point(cos(th), sin(th)) for th in range(0, 2pi, length=20000)[1:end-1]]
BasicCompGeometry._check_distinct_vertices(FloatingKernel(), poly) # Allocates intermediate tuple vector
```
[REMEDIATION]:
```julia
function _check_distinct_vertices(kernel::AbstractKernel, vertices)
    n = length(vertices)
    order = sortperm(vertices; alg = MergeSort, by = p -> (p[1], p[2]))
    for i in 1:n-1
        a = vertices[order[i]]
        b = vertices[order[i + 1]]
        (a[1] == b[1] && a[2] == b[2]) && throw(
            ArgumentError("polygon contains duplicate vertex $(order[i + 1]) (equal to vertex $(order[i]))"),
        )
    end
    return nothing
end
```

---

[SEVERITY]: Moderate
[LOCATION]: src/SweepTriangulation.jl:268, 289
[ROOT CAUSE]: In `_sweep_partition_diagonals`, an $n$-element vector `prev_v = [mod1(i - 1, n) for i in 1:n]` is allocated on the heap, but `prev_v[v]` is accessed only once per event at line 292 as `e_in = prev_v[v]`, which is an inline cyclic decrement `v == 1 ? n : v - 1`. Additionally, `diags = Tuple{Int,Int}[]` is initialized without `sizehint!(diags, n - 3)`, incurring repeated reallocations as partition diagonals are discovered and pushed.
[FAILURE TRACE]:
```julia
poly = [point(cos(th), sin(th)) for th in range(0, 2pi, length=10000)[1:end-1]]
# Triangulating via sweep allocates the redundant prev_v array and resizes diags multiple times
triangulate_sweep(poly)
```
[REMEDIATION]:
```julia
    @inline _prev_v(i::Int, n::Int) = i == 1 ? n : i - 1
    # Replace prev_v[v] with _prev_v(v, n) and sizehint diags:
    diags = Tuple{Int,Int}[]
    sizehint!(diags, n - 3)
```
