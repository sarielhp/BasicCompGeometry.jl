# Hull and diameter algorithms

See [Algorithms and spatial structures](../algorithms.md) for examples and
[Kernels and numerical scope](../kernels.md) for arithmetic limits.
Feature submodule operations use qualified names.

## 2D hull and diameter queries

| Call | Result |
|:-----|:-------|
| `convex_hull(points; kernel=FloatingKernel())` | Counterclockwise `PntSeq` hull, using monotone chain in `O(n log n)` time |
| `exact_diameter(points)` | Maximum pair distance from an `O(n²)` scan |
| `approx_diameter(points, ε)` | WSPD-based `(1 + ε)` diameter approximation |
| `exact_diameter_subspace(points, dir)` | Diameter after orthogonal projection to the complement of unit `dir` |
| `approx_diameter_subspace(points, ε, dir)` | Approximate projected diameter |

The hull accepts a vector of 2D points or an `AbsPntSeq`. Diameter operations
accept an `AbsPntSeq`. Request `ExactKernel()` for exact hull decisions; the
diameter routines use ordinary distance arithmetic.

```@autodocs
Modules = [BasicCompGeometry]
Pages = ["ConvexHull.jl", "Diameter.jl", "MonotoneTriangulation.jl", "VisibilityPolygon.jl"]
Filter = obj -> obj in Main.documented_api
```

## Visibility polygons

| Call | Result |
|:-----|:-------|
| `visibility_polygon(poly, q; kernel=FloatingKernel())` | Compute the star-shaped visibility polygon from vantage point `q` inside or on the boundary of simple polygon `poly` |
| `is_inside(q, poly; kernel=FloatingKernel())` | Test whether query point `q` lies in the interior or on the boundary of `poly` |

## Polygon triangulation

| Call | Result |
|:-----|:-------|
| `triangulate_polygon(vertices; kernel=FloatingKernel(), algorithm=:ear_clipping)` | Triangulate any simple CCW polygon into `n-3` diagonals and `n-2` CCW triangles (`:ear_clipping` or `:seidel`) |
| `triangulate_seidel(vertices; kernel=FloatingKernel(), rng=Random.default_rng())` | Seidel's randomized incremental triangulation in expected `O(n log* n)` time |
| `triangulate_monotone(vertices; kernel=FloatingKernel())` | Linear-time triangulation returning `k-3` internal diagonals and `k-2` CCW triangles |
| `is_x_monotone(vertices)` | Return `true` if polygon is strictly x-monotone |
| `is_y_monotone(vertices)` | Return `true` if polygon is strictly y-monotone |
| `monotone_chains(vertices; kernel=FloatingKernel())` | Upper and lower chains strictly increasing in x |
| `merge_chains(upper, lower)` | Linear-time merge of two chains sorted by x |
| `stack_diagonals(vertices; kernel=FloatingKernel())` | List of `k-3` internal diagonals for x-monotone polygon |
| `polygon_diagonals(vertices; kernel=FloatingKernel())` | List of `n-3` internal diagonals for any simple polygon |

## 3D convex hull

```@autodocs
Modules = [BasicCompGeometry.ConvexHull3D]
Private = false
```

## Longest convex subset

```@autodocs
Modules = [BasicCompGeometry.LongestConvexSubset]
Private = false
```
