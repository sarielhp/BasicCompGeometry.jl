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
Pages = ["ConvexHull.jl", "Diameter.jl", "MonotoneTriangulation.jl", "VisibilityPolygon.jl", "MinkowskiSum.jl"]
Filter = obj -> obj in Main.documented_api
```

## Discrete hulls of convex bodies

The discrete hull of a planar convex body $C$ is $\mathrm{CH}(C \cap \mathbb{Z}^2)$.
`discrete_hull` implements the output-sensitive algorithm of Har-Peled
(*An Output Sensitive Algorithm for Discrete Convex Hulls*, CGTA 10, 1998).
It reads $C$ only through `chord_range` and makes $O(|C'| \log \mathrm{diam}(C))$
calls, where $|C'|$ is the number of hull vertices. For a disk of radius $r$
this is $O(r^{2/3} \log r)$.

| Call | Result |
|:-----|:-------|
| `discrete_hull(C; seed=nothing, start=nothing)` | Hull vertices as `Point2I`, counterclockwise from the lowest-leftmost vertex, with no collinear vertices; empty if $C \cap \mathbb{Z}^2 = \emptyset$ |
| `chord_range(C, o, d)` | `nothing` if the line `o + t*d` misses `C`, else `(⌈t_lo⌉, ⌊t_hi⌋)` |
| `LatticeDisk(center, radius)` | Disk with exact rational center and radius; floats are rounded to multiples of `1/den` (default `2^16`) |
| `ConvexPolygonBody(vertices)` | Convex polygon with exact rational vertices in counterclockwise order |
| `FloatChordBody(f)` | Body given by a floating-point chord function `f(o, d)`; exact only up to the rounding of `f` |
| `find_edge_direction(C, p, L, U)` | Continued-fraction search for the next edge direction from `p` in the lattice cone spanned by `(L, U)` |
| `discrete_hull_start(C, seed)` | Lowest-leftmost lattice point of `C` by unbounded and binary search; needs lattice-connected `C ∩ Z²` |
| `discrete_hull_columns(C[, xs])` | Naive $O(r)$ hull of the lowest and highest lattice points on each column |
| `lattice_convex_hull(points)` | Exact convex hull of integer points |

For other bodies, subtype `AbsLatticeBody`, implement `chord_range`, and pass
either `seed` (a lattice point of `C`) or `start` (a vertex of the hull). The
`chord_range` results determine the hull, so the hull is exact when they are.

```@autodocs
Modules = [BasicCompGeometry]
Pages = ["DiscreteHull.jl"]
Filter = obj -> obj in Main.documented_api
```

## Minkowski sums and convexity

| Call | Result |
|:-----|:-------|
| `minkowski_sum(poly1, poly2; kernel=FloatingKernel(), hull=false)` | Compute the Minkowski sum $P \oplus Q$ in optimal $O(n + m)$ time for convex polygons, or its convex hull with `hull=true` |
| `minkowski_sum_decomposition(poly1, poly2; kernel=FloatingKernel())` | Compute convex pieces whose union forms the Minkowski sum of simple polygons via triangulation |
| `is_convex(poly; kernel=FloatingKernel())` | Test whether polygon `poly` is convex |

## Visibility polygons and kernels

| Call | Result |
|:-----|:-------|
| `visibility_polygon(poly, q; kernel=FloatingKernel())` | Compute the star-shaped visibility polygon from vantage point `q` inside or on the boundary of simple polygon `poly` |
| `polygon_kernel(poly; kernel=FloatingKernel())` | Compute the convex polygon forming the kernel $\mathrm{Ker}(P)$ of points that see the entire polygon |
| `is_star_shaped(poly; kernel=FloatingKernel())` | Test whether simple polygon `poly` has a non-empty kernel |
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
