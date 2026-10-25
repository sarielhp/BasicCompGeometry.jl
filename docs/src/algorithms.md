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

## Simplification hierarchies

`BasicCompGeometry.PntSeqHierarchy` stores simplifications of a curve at decreasing
widths. A palette records the error of a recursive subdivision of the curve, so that all
vertices above a threshold are extracted in time proportional to their number, and a
`SimplificationHierarchy` keeps a sorted set of simplifications for fast lookup by
width. The error being controlled is defined by a `HierarchyMeasure`; the default uses the
Hausdorff distance to the chord, and other measures, such as a Frechet distance, supply
their own error, simplification and certification functions. See
[Simplification hierarchies](reference/hierarchy.md).

```julia
import BasicCompGeometry as CG
H = CG.PntSeqHierarchy
curve = CG.rand_pnt_seq(2, Float64, 200)
h = H.hierarchy_init(curve)
H.hierarchy_build!(h)
simplified, width = H.hierarchy_approx(h, 0.3)
```

## Polygon triangulation

`triangulate_polygon` triangulates an arbitrary simple counterclockwise polygon with $n \ge 3$ vertices.
It returns a `PolygonTriangulation` (`MonotoneTriangulation` is a compatibility alias) containing
`n - 3` internal diagonals and `n - 2` counterclockwise triangles.

Supported algorithm options via keyword `algorithm`:
- `:ear_clipping` (default): Fast ear clipping with reflex-vertex acceleration running in $O(n r)$ time where $r$ is the count of reflex vertices.
- `:seidel`: Raimund Seidel's randomized incremental algorithm running in expected $O(n \log^* n)$ time, also available via `triangulate_seidel(vertices; kernel=FloatingKernel(), rng=Random.default_rng())`.

For strictly $x$-monotone or $y$-monotone polygons (checked via `is_x_monotone` and `is_y_monotone`),
triangulation executes in linear time `O(n)`.

```@example algorithms
# Non-monotone U-shaped polygon
u_poly = [
    point(0.0, 0.0),
    point(4.0, 0.0),
    point(4.0, 4.0),
    point(3.0, 4.0),
    point(3.0, 1.0),
    point(1.0, 1.0),
    point(1.0, 4.0),
    point(0.0, 4.0)
]
res = triangulate_polygon(u_poly; algorithm = :seidel)
(diagonals = res.diagonals, triangles = res.triangles)
```

`triangulate_monotone` decomposes a simple counterclockwise x-monotone polygon into
triangles in deterministic linear time `O(k)` using the monotone stack algorithm.
It returns a `PolygonTriangulation` with `k - 3` internal diagonals and `k - 2`
counterclockwise triangles:

```@example algorithms
poly = [
    point(0.0, 0.0),
    point(1.0, -1.0),
    point(3.0, -1.0),
    point(4.0, 0.0),
    point(2.0, 2.0)
]
res = triangulate_monotone(poly)
(diagonals = res.diagonals, triangles = res.triangles)
```

`monotone_chains` and `merge_chains` separate and order the upper and lower boundary
chains in `O(k)` time. See [the algorithm reference](reference/algorithms.md).
