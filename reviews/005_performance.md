# Systems Code Review Report #005 (Lens: performance)

- **Date**: 2026-11-13
- **Auditor**: Gemini 3.8 Flash (Tier 0 Workhorse) via `review-audit`
- **Focus Lens**: `performance`
- **Model Tier**: `TIER0`
- **Backend**: `gemini`
- **Scope**: `src/BBT.jl src/WSPD.jl src/Diameter.jl src/MetricSpace.jl src/RMQ.jl`
- **Status**: Action Required

---

Auditing via Gemini Flash (bws run) [Profile: performance]...
[SEVERITY]: Major
[LOCATION]: [src/RMQ.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/RMQ.jl#L36-L50) and [src/RMQ.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/RMQ.jl#L90-L100)
[ROOT CAUSE]: `_sparse_levels(keys, less, minima::Vector{Int})` has two compounding performance defects:
1. `half` is captured by the comprehension generator `[ _better_index(keys, less, previous[i], previous[i+half]) for i = 1:length(minima)-2half+1 ]` and reassigned (`half *= 2`) in the outer `while` loop, forcing Julia to wrap `half` in a `Core.Box` (`Any`) and infer `row::ANY` and an untyped `Base.Generator`.
2. Neither `build_rmq` nor `_sparse_levels` parameterizes `less::F` (and neither calls `less` directly, only forwarding it to `_better_index`), so Julia's specialization heuristic skips specializing `_sparse_levels` on `less`, resulting in dynamic dispatch of `_better_index(keys, less::Function, ...)` on every entry of every sparse-table level.
[FAILURE TRACE]: Calling `BasicCompGeometry._sparse_levels(keys, isless, collect(1:100_000))` on a 100,000-element `Vector{Float64}` allocates **82,009,008 bytes** in 0.058s (versus **11,754,400 bytes** in 0.024s when rewritten without the `Core.Box` comprehension and with `keys::V, less::F where {V,F}`). Eliminating only the `Core.Box` without parameterizing `less::F` still allocates **58,760,672 bytes** due to per-element dynamic dispatch.
[REMEDIATION]:
```julia
function _sparse_levels(keys::V, less::F, minima::Vector{Int}) where {V,F}
    isempty(minima) && return Vector{Int}[]
    levels = [minima]
    half = 1
    n = length(minima)
    while half <= n ÷ 2
        previous = last(levels)
        m = n - 2 * half + 1
        row = Vector{Int}(undef, m)
        @inbounds for i in 1:m
            row[i] = _better_index(keys, less, previous[i], previous[i+half])
        end
        push!(levels, row)
        half *= 2
    end
    return levels
end
```

[SEVERITY]: Major
[LOCATION]: [src/RMQ.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/RMQ.jl#L90-L100)
[ROOT CAUSE]: `build_rmq(keys::AbstractVector; less = isless, backend::Symbol = :blocked)` branches on the runtime `Symbol` value `backend` to return either `SparseRMQ` or `BlockRMQ`. Because `backend` is not constant-propagated across the keyword wrapper, `build_rmq` has a value-dependent union return type `Union{BlockRMQ{V,F}, SparseRMQ{V,F}}` (and does not parameterize `less::F`). Downstream `build_chain_index` in [src/ChainIndex.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/ChainIndex.jl#L77-L88) constructs `ChainIndex` from four `build_rmq` calls, inferring a 16-way union `ChainIndex{V, A, B, C, D} where {A<:Union{...}, B<:Union{...}, C<:Union{...}, D<:Union{...}}` that exceeds Julia's union-splitting limit.
[FAILURE TRACE]:
1. `@code_warntype BasicCompGeometry.build_rmq(rand(1000))` infers `Body::Union{BlockRMQ{Vector{Float64}, typeof(isless)}, SparseRMQ{Vector{Float64}, typeof(isless)}}`.
2. In a function that builds a `ChainIndex` via `build_chain_index(pts)` and queries `leftmost(w)` in a loop, every `O(1)` query performs dynamic dispatch and heap-allocates 32 bytes (`43,632` bytes total for 1,000 queries; `331,632` bytes total for 10,000 queries).
[REMEDIATION]:
```julia
@inline function _build_rmq(keys::V, less::F, ::Val{:blocked}) where {V,F}
    block_size = max(1, _floor_log2(max(1, length(keys))))
    masks, minima = _block_masks(keys, less, block_size)
    return BlockRMQ(keys, less, block_size, masks, _sparse_levels(keys, less, minima))
end

@inline function _build_rmq(keys::V, less::F, ::Val{:sparse}) where {V,F}
    return SparseRMQ(keys, less, _sparse_levels(keys, less, collect(eachindex(keys))))
end

Base.@constprop :aggressive function build_rmq(
    keys::V;
    less::F = isless,
    backend::Union{Symbol,Val} = Val(:blocked),
) where {V<:AbstractVector,F}
    Base.require_one_based_indexing(keys)
    b = backend isa Val ? backend : Val(backend)
    b isa Val{:blocked} || b isa Val{:sparse} ||
        throw(ArgumentError("unknown RMQ backend: $backend"))
    return _build_rmq(keys, less, b)
end
```

[SEVERITY]: Moderate
[LOCATION]: [src/BBT.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/BBT.jl#L94-L129)
[ROOT CAUSE]: `get_min_max_orig_index` uses dynamic reflection `getfield(node, field_name)` and `setfield!(node, field_name, val)` with a runtime `Symbol` (`:min_index` or `:max_index`) across a recursive traversal. Because `field_name` is not a compile-time constant across the recursion, `val = getfield(node, field_name)` infers as a 7-way `Union` over every field type of `Node{D,T,S,V}` (`Union{Nothing, UnitRange{Int}, Bool, S, Int, Node{...}, BBox{...}}`), `setfield!` boxes `Int` on the heap at every uncached node, and the leaf/unexpanded branch allocates a `Base.Generator` for `min_max_func(orig_index(tree.PS, i) for i in node.r)`. This is called on every popped pair in `WSPD.expand_bichromatic!` ([src/WSPD.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/WSPD.jl#L173-L174)).
[FAILURE TRACE]: Building and expanding a `BBT.Tree` over 2,000 `Point2F` points and calling `BBT.get_min_orig_index(tree, tree.root)` once (after warmup) allocates **35,072 bytes**, whereas direct field access (`node.min_index`) with a tight `for i in node.r` loop allocates **0 bytes**.
[REMEDIATION]:
```julia
function _range_ext_orig_index(tree::Tree, r::UnitRange{Int}, ::typeof(min))::Int
    best = typemax(Int)
    @inbounds for i in r
        oi = orig_index(tree.PS, i)
        oi < best && (best = oi)
    end
    return best
end

function _range_ext_orig_index(tree::Tree, r::UnitRange{Int}, ::typeof(max))::Int
    best = typemin(Int)
    @inbounds for i in r
        oi = orig_index(tree.PS, i)
        oi > best && (best = oi)
    end
    return best
end

function get_min_orig_index(tree::Tree{D,T,S,V}, node::Node{D,T,S,V})::Int where {D,T,S,V}
    val = node.min_index
    val !== nothing && return val
    idx = if !node.f_leaf && !isnothing(node.left) && !isnothing(node.right)
        min(get_min_orig_index(tree, node.left), get_min_orig_index(tree, node.right))
    else
        _range_ext_orig_index(tree, node.r, min)
    end
    node.min_index = idx
    return idx
end

function get_max_orig_index(tree::Tree{D,T,S,V}, node::Node{D,T,S,V})::Int where {D,T,S,V}
    val = node.max_index
    val !== nothing && return val
    idx = if !node.f_leaf && !isnothing(node.left) && !isnothing(node.right)
        max(get_max_orig_index(tree, node.left), get_max_orig_index(tree, node.right))
    else
        _range_ext_orig_index(tree, node.r, max)
    end
    node.max_index = idx
    return idx
end
```

[SEVERITY]: Moderate
[LOCATION]: [src/WSPD.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/WSPD.jl#L31-L38), [src/WSPD.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/WSPD.jl#L51-L58), and [src/WSPD.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/WSPD.jl#L130-L143)
[ROOT CAUSE]: Three type-parameter/type-stability defects in `WSPD`:
1. `WSPD.init(pnts::AbsPntSeq{D,T}, sep::S)` binds `S` from `sep` instead of from the distance type `S_tree` of `tree = BBT.Tree_init(pnts)` (`typeof(dist(pnts[1], pnts[1]))`). When `typeof(sep) != S_tree`, constructing `PD{D,T,S,V}` fails to convert `BBT.Tree{D,T,S_tree,V}` to `BBT.Tree{D,T,S,V}`.
2. `separation(pair::WSPair)` returns the `Float64` constant `MAX_SEP = 1e18` when `pair.dist <= 0.0` and `pair.max_sides_diam / pair.dist` (type `S`) otherwise, yielding a value-dependent return type `Union{Float32, Float64}` for `Float32` trees in `expand!` and `expand_bichromatic!`.
3. `mutable struct PD{D,T,S,V}` declares `pnts::AbsPntSeq{D,T}` using the abstract type `AbsPntSeq{D,T}` instead of a type parameter.
[FAILURE TRACE]:
- Calling `approx_diameter(P, 0.1)` ([src/Diameter.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/Diameter.jl#L30)) where `P = PntSeq([Point{2,Float32}(...)])` (or `WSPD.init(P_float64, 2)`) throws `MethodError: Cannot convert an object of type BBT.Tree{2, Float32, Float32, ...} to an object of type BBT.Tree{2, Float32, Float64, ...}`.
- `@code_warntype BasicCompGeometry.WSPD.separation(first(WSPD.init(ps32, 2.0f0).heap))` infers `Body::Union{Float32, Float64}`.
[REMEDIATION]:
```julia
@inline function separation(pair::WSPair{D,T,S,V})::S where {D,T,S,V}
    if pair.dist <= zero(S)
        return S(MAX_SEP)
    end
    return pair.max_sides_diam / pair.dist
end

mutable struct PD{D,T,S,V,P<:AbsPntSeq{D,T}}
    pnts::P
    finals::Vector{WSPair{D,T,S,V}}
    heap::BinaryHeap{WSPair{D,T,S,V},WSPDOrder}
    tree::BBT.Tree{D,T,S,V}
    sep::S
    hash_pairs::Dict{PairInt,Bool}
end

function init(pnts::P, sep::Real) where {D,T,P<:AbsPntSeq{D,T}}
    tree = BBT.Tree_init(pnts)
    S = typeof(tree.root.diam)
    V = typeof(Points(pnts))
    PairT = WSPair{D,T,S,V}
    heap = BinaryHeap{PairT}(WSPDOrder(), PairT[])
    W = PD{D,T,S,V,P}(pnts, PairT[], heap, tree, S(sep), Dict{PairInt,Bool}())
    push_pair!(W, tree.root, tree.root)
    return W
end
```
