# Algorithms and spatial structures

## Hulls and diameters

The 2D hull uses the monotone chain algorithm and returns a `PntSeq` in
counterclockwise order. The 3D hull uses an incremental algorithm and returns
vertex indices and triangular faces; it retains floating-point decisions.

`exact_diameter` scans all pairs in quadratic time and returns a distance.
Here, "exact" describes the search over all pairs; it does not select exact
predicate arithmetic. `approx_diameter(points, ε)` uses a WSPD for a
`(1 + ε)` approximation. The subspace variants project distances into the
subspace orthogonal to a supplied unit direction.

```@example algorithms
using BasicCompGeometry

points = PntSeq([point(0.0, 0.0), point(3.0, 0.0), point(3.0, 4.0), point(0.0, 4.0)])
hull = convex_hull(points; kernel=ExactKernel())
(hull_vertices=length(hull), diameter=exact_diameter(points),
 approximation=approx_diameter(points, 0.1))
```

See [the algorithm reference](reference/algorithms.md) for these operations and
the longest convex subset algorithms. `MVBB` supplies approximate 3D oriented
bounding boxes and a 2D rotating-calipers minimum-area rectangle operation;
see [the spatial reference](reference/spatial.md).

## Nearest neighbors with BBT

Construct a bounding box tree with `BBT.Tree_init`. Nodes can expand lazily;
`BBT.Tree_fully_expand` refines the entire tree.

```@example algorithms
tree = BBT.Tree_init(points)
BBT.Tree_fully_expand(tree)
query = point(2.8, 3.9)
distance, nearest, original_index = BBT.exact_naive_scan(tree, query)
(distance=distance, point=nearest, index=original_index)
```

Queries return `(distance, point, original_index)`:

| Query | Behavior |
|:------|:---------|
| `BBT.exact_naive_scan(tree, query)` | Linear scan of all points |
| `BBT.approx_nn(tree, query, c)` | Search with approximation factor `c ≥ 1` |
| `BBT.silly_nn(tree, query)` | Tree descent heuristic |
| `BBT.hybrid_nn(tree, query, c)` | Combines descent and approximate search |

The descent heuristic has no approximation-factor guarantee. Performance of
tree searches depends on the data and dimension. `BBT.Tree_draw` requires
`Cairo` and `Colors`; see [Cairo and LaTeX](cairo.md).

## Well-separated pairs

`WSPD.init(points, sep)` creates a construction state; `WSPD.expand!` returns
its finalized node pairs. The parameter is the maximum ratio of the larger
node-box diameter to the distance between the two boxes: smaller values
request stronger separation.

```@example algorithms
state = WSPD.init(points, 0.25)
pairs = WSPD.expand!(state)
(pair_count=length(pairs), representatives=WSPD.get_reps(state, first(pairs)))
```

`WSPD.top_refine!` refines the active pair with the highest diameter bound.
`expand_bichromatic!` restricts construction to pairs crossing an input-index
cut. Representative and index helpers are listed in
[the spatial reference](reference/spatial.md).

For polygonal-chain simplification and resampling, see
[Geometry and point sequences](geometry.md). For greedy permutations, see
[Finite metric spaces](metrics.md).
