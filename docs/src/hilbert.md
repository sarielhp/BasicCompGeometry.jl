# Hilbert spaces and virtual points

A [`HilbertSpace`](@ref) is a dataset whose points are known only through their
squared distances: the feature space of a kernel (possibly infinite-dimensional), a
precomputed distance function, or high-dimensional coordinates read once. Curves
through its points are sequences of virtual points ([`VPoint`](@ref)), which have no
coordinates; distances, convex combinations and closest points on them are computed
from the squared distances of the input points, each computed at most once and cached.

```@example hilbert
using BasicCompGeometry, LinearAlgebra
using BasicCompGeometry: local_frame, cache_stats

# Twelve frames of a toy "video", compared under a Gaussian kernel. The squared
# distance in feature space is 2 - 2 exp(-|x - y|^2 / s^2), written with expm1 so
# that close frames do not lose digits to cancellation.
frames = [[sin(0.3t + k) for k in 1:50] for t in 1:12]
gauss2(x, y) = -2expm1(-sum(abs2, x - y) / 25)
S = HilbertSpace(frames, gauss2)

P = vcurve(S, 1:6)     # a curve through frames 1..6
Q = vcurve(S, 7:12)
x = convex_comb(P[2], P[3], 0.25)   # a point on an edge of P
dist(x, Q[1]), geom_length(P)
```

## Constructing a space

| Call | Squared distance of items `x`, `y` |
|:-----|:-----------------------------------|
| `HilbertSpace(items, dist2)` | `dist2(x, y)`, given directly |
| `HilbertSpace(items; kernel = k)` | `k(x, x) + k(y, y) - 2k(x, y)`, clamped at 0 |
| `HilbertSpace(points)` | Euclidean, for coordinate vectors |
| `HilbertSpace(vpoints)` | the squared distances of virtual points of another space |

Give `dist2` directly when it can be computed accurately; through a kernel, the
squared distance of close points is a small difference of large numbers. The
function must be pure, thread-safe, and come from a Hilbert space (a
positive-definite kernel), and return finite non-negative values.

`vpoint(S, i)` is the input point `i`, `vcurve(S, idx)` the curve through the input
points `idx`, and `VPoint(S, a, b, t)` the point `(1 - t) x_a + t x_b` on the segment
of any two input points. `BasicCompGeometry.is_generator(x)` tells whether a virtual
point is an input point or derived from two of them. A space is also a finite metric
space (`size`, `dist(S, i, j)`).

## Geometry

`dist`, `dist_sq` and `convex_comb` work on virtual points of one space, and so do
`geom_length`, `prefix_lengths`, `push_distinct!` and `split_edges` on their curves.
Two building blocks let Euclidean code run on virtual points:

- [`local_frame`](@ref BasicCompGeometry.local_frame) places up to four virtual points, involving at most four input
  points, isometrically in ℝ³; Euclidean code then runs on the returned `Point{3}`s.
- `BasicCompGeometry.inner_diff(x, y, z, w)` is `⟨x - y, z - w⟩`; `x - y` is a formal
  difference accepted by `dot`.

```@example hilbert
A, B, C = local_frame(P[1], P[2], Q[1])
dist(A, C) ≈ dist(P[1], Q[1])
```

A convex combination of two virtual points must lie on a segment of two input points:
points of different segments throw an `ArgumentError`. To go further, re-base the
points: `HilbertSpace(points)` makes them the input points of a new space.

Functions that need coordinates (`dim`, `translate!`, `scale!`, bounding boxes, `at`,
`sample_uniformly`, `simplify`, `Matrix`, file output, the simplification hierarchy)
do not accept virtual curves.

## Precision

Computing from squared distances is less accurate than computing from coordinates when
points are close. For a point at distance δ from a segment whose points are at distance
at most D, the error is about ε D² / δ, and a true distance 0 comes out as about √ε D
(about `1e-8 D`). Use thresholds relative to `BasicCompGeometry.vscale` (the largest
distance among the input points involved), not absolute ones.

## The cache

Squared distances are cached in blocks of 256 × 256 pairs (512 KB), allocated when
first touched; a lookup is two loads and a comparison. Threads may share a space: a
value two threads need at the same time may be computed twice, with the same result.
`HilbertSpace(...; max_bytes = m)` limits the cache (further values are computed
without being stored), `empty!(S)` frees it, and `cache_stats(S)` counts the calls of
the distance function.

```@example hilbert
cache_stats(S)
```
