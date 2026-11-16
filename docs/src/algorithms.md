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

Construct a bounding box tree with `BBT.Tree_init` (`Kd_tree` is an alias for
`BBT`). Nodes can expand lazily; `BBT.Tree_fully_expand` refines the entire tree,
and each node stores a `parent` link accessible via `parent(node)` or
`BBT.parent_node(node)`.

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
- `:sweep` (default): Deterministic plane-sweep monotone partitioning and triangulation running in worst-case $O(n \log n)$ time, also available via `triangulate_sweep(vertices; kernel=FloatingKernel())`.
- `:auto`: The small-polygon fast path of `:ear_cut` (reflex vertices in a 64-bit mask, adapted from [Triage](https://github.com/sarielhp/triage)) for $n \le 64$, where it is several times faster than `:sweep`, and `:sweep` for larger polygons. It is opt-in, so that the default result does not depend on the polygon size.
- `:seidel`: Raimund Seidel's randomized incremental algorithm running in expected $O(n \log^* n)$ time, also available via `triangulate_seidel(vertices; kernel=FloatingKernel(), rng=Random.default_rng())`.
- `:ear_cut` (also accepted as `:ear_clipping`, its pre-0.5.15 name): Ear clipping in $O(n r)$ worst-case time, where $r$ is the count of reflex vertices, and subquadratic in practice; also available via `triangulate_ear_cut(vertices; kernel=FloatingKernel())`. It takes the small-polygon path for $n \le 64$, a linear-time path when `is_x_monotone` or `is_y_monotone` holds, and otherwise indexes the active reflex vertices in an unshrunk-cell Kd-tree (a plain list when $r \le 16$).
- `:ear_cut_u_kd_tree`: The Kd-tree ear clipping of `:ear_cut` on every input, without its fast paths: a flat array-backed 2D `Kd_tree` (`BBT.FlatUTree2D`) over active reflex vertices with uncapped cell-containment parent climbing (worst case $O(n r)$, about $n^{1.2}$ on random star polygons), also available via `triangulate_ear_cut_u_kd_tree(vertices; kernel=FloatingKernel())`.
- `:ear_cut_z_order` and `:ear_cut_kd_tree` are deprecated aliases of `:ear_cut` and `:ear_cut_u_kd_tree` (the Morton Z-order and shrinking-box Kd-tree indices were never faster).

Every returned triangle has non-negative orientation under the given kernel, and with
`ExactKernel()` the triangles tile the polygon exactly; input that cannot be triangulated
that way (zero area, crossing edges, or a polygon numerically degenerate for the kernel)
raises an `ArgumentError` instead of returning clockwise or overlapping triangles. Weakly
simple polygons (a vertex touching another edge) are triangulated by the four ear-cutting
algorithms; `:sweep` and `:seidel` may reject some of them with an `ArgumentError`. For
numerically delicate input, use `kernel = ExactKernel()`.

When calling `triangulate_polygon(...; algorithm = :ear_cut)`, strictly $x$-monotone or $y$-monotone polygons (checked via `is_x_monotone` and `is_y_monotone`)
use a linear-time `O(n)` fast path; `triangulate_monotone` also runs directly in $O(n)$ time.

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

## Visibility polygons

`visibility_polygon` computes the star-shaped visibility polygon ``\mathrm{Vis}(q, P)`` consisting of all points inside a counterclockwise simple polygon ``P`` visible from a vantage point ``q \in P``. Under floating-point kernels a vantage point within
`1e-10` times the polygon's bounding-box extent of a vertex (or else an edge) is treated as
lying on it, and the visibility polygon of that boundary point is returned; near a reflex
vertex that can differ substantially from the exact answer for the unsnapped point. Use
`ExactKernel()` when the exact vantage point matters.

```@example algorithms
p_u = [
    point(0.0, 0.0),
    point(4.0, 0.0),
    point(4.0, 4.0),
    point(3.0, 4.0),
    point(3.0, 1.0),
    point(1.0, 1.0),
    point(1.0, 4.0),
    point(0.0, 4.0)
]

vis = visibility_polygon(p_u, point(2.0, 0.5))
(vertices_count = length(vis), twice_area = twice_area(vis))
```

Under `kernel = ExactKernel()`, exact rational coordinates are used and preserved throughout.
Point containment for polygons is verified via `is_inside(q, poly)`.

## Polygon kernel and star-shaped polygons

`polygon_kernel` computes the kernel of a simple polygon $\mathrm{Ker}(P) = \{ q \in P \mid \forall p \in P, \overline{qp} \subset P \}$, which is the convex region of points that can see the entire polygon. If the kernel is non-empty, the polygon is star-shaped, tested by `is_star_shaped(poly)`.

```@example algorithms
# L-shaped polygon: kernel is the corner square [0, 1] x [0, 1]
p_l = [
    point(0.0, 0.0), point(3.0, 0.0), point(3.0, 1.0),
    point(1.0, 1.0), point(1.0, 3.0), point(0.0, 3.0)
]
ker = polygon_kernel(p_l)
(star_shaped = is_star_shaped(p_l), kernel_area = twice_area(ker) / 2)
```

## Discrete hulls

`discrete_hull(C)` computes the convex hull of the lattice points of a convex
body `C` without enumerating them, by Har-Peled's output-sensitive algorithm.
For a disk of radius $r$ it takes $O(r^{2/3} \log r)$ time, against $O(r)$ for
the naive column scan `discrete_hull_columns`. A disk of radius $10^7$ takes
under 0.1 s. Disks and rational polygons are handled exactly; see
[the algorithm reference](reference/algorithms.md#Discrete-hulls-of-convex-bodies)
for other bodies, and `examples/discrete_hull_disk.jl` for a full example.

```@example algorithms
D = LatticeDisk((0.3, 0.7), 1000)
H = discrete_hull(D)
(vertices = length(H), ratio = length(H) / 1000^(2/3),
 matches_naive = H == discrete_hull_columns(D))
```

## Minkowski sum of polygons

The Minkowski sum of two sets $P$ and $Q$ is defined as $P \oplus Q = \{ p + q \mid p \in P, q \in Q \}$.

For two convex polygons with $n$ and $m$ vertices, `minkowski_sum(poly1, poly2)` computes the exact convex Minkowski sum in optimal $O(n + m)$ time by merging their directed edges in cyclic polar angle order. It supports both `FloatingKernel()` and `ExactKernel()`, as well as points and segments (producing translations, stadium shapes, and swept polygons).

Convexity is verified by `is_convex(poly)`. For non-convex polygons:
- `minkowski_sum(poly1, poly2; hull=true)` computes the convex hull $\mathrm{CH}(P \oplus Q) = \mathrm{CH}(P) \oplus \mathrm{CH}(Q)$ in $O(n \log n + m \log m)$ time.
- `minkowski_sum_decomposition(poly1, poly2; algorithm=:sweep)` decomposes non-convex polygons into convex pieces (triangles) via [`triangulate_polygon`](@ref) (defaulting to deterministic $O(n \log n)$ `:sweep` triangulation), returning a vector of convex polygons whose union is $P \oplus Q$.

```@example algorithms
# Minkowski sum of a square and a triangle:
sq = [point(0.0, 0.0), point(1.0, 0.0), point(1.0, 1.0), point(0.0, 1.0)]
tri = [point(0.0, 0.0), point(1.0, 0.0), point(0.0, 1.0)]
m_sum = minkowski_sum(sq, tri)
(vertices = length(m_sum), area = twice_area(m_sum) / 2)
```
