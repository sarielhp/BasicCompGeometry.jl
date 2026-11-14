# Systems Code Review Report #006 (Lens: performance)

- **Date**: 2026-11-13
- **Auditor**: Claude / Codex (Tier 1 Standard) via `review-audit`
- **Focus Lens**: `performance`
- **Model Tier**: `STANDARD`
- **Backend**: `auto`
- **Scope**: `src/ConvexHull.jl src/ConvexHull3D.jl src/MVBB.jl src/DiscreteHull.jl src/LongestConvexSubset.jl src/Halfplanes.jl`
- **Status**: Action Required

---

Auditing via Codex CLI (gpt-5.6-sol, tier: standard) [Profile: performance]...
[SEVERITY]: Moderate  
[LOCATION]: `src/LongestConvexSubset.jl:45-58` — `LCSContext` and `RotateOrder`  
[ROOT CAUSE]: Both structures erase the concrete point-sequence subtype by storing `PS::AbsPntSeq{2,T}`. `RotateOrder.lt` is invoked throughout every angular sort, while the DP repeatedly accesses `lct.PS`; Julia therefore cannot infer the concrete sequence representation from either structure. This introduces dynamic field access/dispatch in the dominant \(O(n^3 \log n)\) workflow.  
[FAILURE TRACE]: Call `compute_largest_convex_subset` with an ordinary `PntSeq` containing several hundred points. `pt_init` constructs a `RotateOrder`, and every sort comparison executes `o.PS[a]` through the abstract field. Subsequent DP iterations similarly retrieve points through `LCSContext.PS` millions of times.  
[REMEDIATION]:
```julia
mutable struct LCSContext{T,P<:AbsPntSeq{2,T}}
    PS::P
    turns::Vector{PointTurn}
    Lens::Vector{Int}
end

struct RotateOrder{T,P<:AbsPntSeq{2,T}} <: Base.Order.Ordering
    PS::P
    base::Point{2,T}
end

# Allow constructor inference to retain P.
lct = LCSContext(PS, pt, Int[])
ro = RotateOrder(PS, PS[loc])
```

[SEVERITY]: Moderate  
[LOCATION]: `src/LongestConvexSubset.jl:64-93, 334-340` — `pt_init`, `lcs_init`, and the target loop in `compute_largest_convex_subset`  
[ROOT CAUSE]: The angular ordering and \(O(n^2)\) `PointTurn` topology are independent of `target`, but both upper and lower contexts are reconstructed for every one of the \(n\) targets. This repeats all tail-vector allocation, initialization, and sorting \(2n\) times, producing \(\Theta(n^3)\) temporary integer/boolean storage traffic and repeating \(\Theta(n^3\log n)\) comparison work that can be precomputed once.  
[FAILURE TRACE]: For 1,000 input points, each `lcs_init` constructs approximately 500,000 tail entries across three arrays and sorts every tail. `compute_largest_convex_subset` performs that initialization 2,000 times even though only DP values depend on the target. Allocation and garbage collection dominate before the DP completes.  
[REMEDIATION]:
```julia
lct_upper = lcs_init(work_ps)
lct_lower = lcs_init(ps_reflected)

for target in 1:n
    # empty! retains the already allocated Lens capacity.
    empty!(lct_upper.Lens)
    empty!(lct_lower.Lens)
    foreach(pt -> pt.f_terminal = false, lct_upper.turns)
    foreach(pt -> pt.f_terminal = false, lct_lower.turns)

    compute_dp_table!(lct_upper, target)
    compute_dp_table!(lct_lower, target)

    best_start = _best_start(lct_upper, lct_lower)
    # ...
end
```
The existing DP overwrites the relevant `tvals` and `jump` entries, so the sorted tail topology can safely be reused.

[SEVERITY]: Moderate  
[LOCATION]: `src/MVBB.jl:246` and `src/MVBB.jl:352` — priority queues in `approx_diam` and `approx_diam_subspace`  
[ROOT CAUSE]: `PriorityQueue{Tuple{Node,Node},T}` uses the unparameterized `Node` UnionAll type. Both the key type and `Tuple{Node,Node}` are non-concrete, discarding `Node{D,T,S,V}` on every queue insertion and removal. Consequently, the dominant refinement loop cannot infer the node fields or specialize downstream calls reliably.  
[FAILURE TRACE]: Run either diameter routine on a large point cloud with a small approximation tolerance, causing thousands of queue operations. Each `first(pq)` returns nodes whose concrete parameters have been erased. Repository inspection confirms both `isconcretetype(Node)` and `isconcretetype(Tuple{Node,Node})` are false.  
[REMEDIATION]:
```julia
root = tree.root
N = typeof(root)
pq = PriorityQueue{Tuple{N,N},T}(Base.Order.Reverse)
pq[(root, root)] = max_dist(root.bb, root.bb)
```
Apply the same change to `approx_diam_subspace`.

[SEVERITY]: Moderate  
[LOCATION]: `src/ConvexHull3D.jl:245-283`, especially lines 272-273  
[ROOT CAUSE]: Every newly created horizon face concatenates two conflict lists into a fresh vector and then runs `unique!`, which also requires duplicate-tracking state. Neighboring horizon faces repeatedly copy heavily overlapping lists. This creates allocation and hashing churn inside the incremental hull’s central nested loop.  
[FAILURE TRACE]: Insert points sampled from a dense spherical shell. Each insertion exposes several faces, and their conflict lists contain many of the same unprocessed points. Every horizon edge executes `[f.visible_points; f_nonvis.visible_points]`, repeatedly copying those indices before visibility testing.  
[REMEDIATION]:
```julia
# Allocate once in convex_hull_3d and pass into _replace_visible_faces!.
candidate_seen = BitSet()

# For each new horizon face:
empty!(candidate_seen)
for list in (f.visible_points, f_nonvis.visible_points), pt in list
    pt in candidate_seen && continue
    push!(candidate_seen, pt)

    if pt != p_idx && is_visible(new_f, pnts[pt])
        push!(new_f.visible_points, pt)
        push!(point_conflicts[pt], nf_idx)
    end
end
```
This retains linear duplicate elimination while reusing scratch storage and avoiding one concatenated vector per horizon edge.
