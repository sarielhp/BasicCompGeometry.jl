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
Pages = ["ConvexHull.jl", "Diameter.jl"]
Filter = obj -> obj in Main.documented_api
```

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
