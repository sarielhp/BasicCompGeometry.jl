# Finite metric spaces

`MetricSpace` separates indexed distance calculations from point storage.
Its `AbsFMS` interface represents a finite collection by indices `1:n`:

- `size(space)` returns the number of points as an integer.
- `dist(space, i, j)` returns the distance between two indexed points.

`MetricSpace.metric` is a convenience wrapper for `dist`.

## Representations

| Type | Representation |
|:-----|:---------------|
| `MetricSpace.PointsSpace` | Vector of objects supporting `dist(p, q)` |
| `MetricSpace.MPointsSpace` | Matrix columns with Euclidean distances |
| `AbsPntSeq` | Point sequence with Euclidean distances between vertices |
| `MetricSpace.PermutMetric` | Reordered indices over an existing space |
| `MetricSpace.SpherePMetric` | Angular separation about a fixed base point |

`MetricSpace.original` recovers source indices from a permutation;
`MetricSpace.swap!` exchanges positions. `MetricSpace.dist_real` measures a
distance from an indexed point to an external point for supported representations.

```@example metrics
import BasicCompGeometry as CG

space = CG.MetricSpace.MPointsSpace([0.0 3.0 0.0; 0.0 4.0 4.0])
permuted = CG.MetricSpace.PermutMetric(space, [3, 1, 2])
(count=size(space), distance=CG.dist(space, 1, 2),
 original=CG.MetricSpace.original(permuted, 1))
```

## Greedy permutations

`greedy_permutation_naive` selects successive points farthest from the selected
centers, using Gonzalez's algorithm. It returns index, distance, and nearest
center arrays. For a request of `k` centers, use the first `k` returned indices:

```@example metrics
indices, distances, nearest_centers = CG.MetricSpace.greedy_permutation_naive(space, 2)
indices[1:2]
```

For custom spaces, subtype `CG.MetricSpace.AbsFMS` and implement `Base.size`
and `CG.dist` for the new type.
