# BasicCompGeometry.jl

*A comprehensive library for basic computational geometry in Julia.*

`BasicCompGeometry.jl` provides efficient, high-dimensional primitives and algorithms for computational geometry, leveraging Julia's multiple dispatch and the `StaticArrays.jl` ecosystem.

A core concept in the library is the `AbsPntSeq`, which is treated as a **sequence of points**. This can represent a simple point set, a polygonal chain, or a classical closed polygon. `AbsPolygon` is provided as an alias for backward compatibility.

## Features

- **Multi-Dimensional Primitives**: Support for Points, Segments, Lines, Point Sequences (PntSeq), and Axis-Aligned Bounding Boxes (BBox) in any dimension.
- **Zero-Copy Matrix Integration**: Use `MatPntSeq` to treat columns of a matrix as points without copying memory.
- **Multiple Coordinate Types**: Works seamlessly with `Float64`, `Int64`, and other numeric types.
- **Geometric Predicates**: Fast checks for turns (left/right) and containment.
- **Curve Algorithms**: Hausdorff distance-based simplification and uniform resampling.
- **Spatial Decomposition**: Efficient Bounding Box Trees (BBT) and Well-Separated Pairs Decompositions (WSPD).
- **2D Transformations**: Efficient translation and rotation for planar geometry.
- **Visualization (Optional)**: High-quality PDF visualization for trees and polygons via `Cairo.jl` and `Colors.jl` (only loaded if these packages are present).

## Algorithms

The library implements a variety of classic and modern geometric algorithms:

- **Convex Hull**: 2D (Monotone Chain) and 3D (Gift-wrapping).
- **Diameter**: Exact $O(n^2)$ and $(1+\epsilon)$-approximation via WSPD.
- **Nearest Neighbor Search**: Exact, $c$-approximate (via BBT), and "silly" (fast descent) versions.
- **Spatial Decomposition**: WSPD, Bounding Box Tree (BBT), and Minimum Volume Bounding Box (MVBB).
- **Curve Processing**: Hausdorff distance-based simplification and uniform arc-length resampling.

- **Metric Space Algorithms**: Greedy permutation (incremental furthest-point sampling).
- **Advanced Optimization**: Longest Convex Subset calculation.
- **Geometric Predicates**: Fast `turn_sign`, `is_left_turn`, `is_right_turn`, and `is_collinear` checks.

## Data Structures

`BasicCompGeometry` provides highly optimized, type-safe data structures:

- **Geometric Primitives**: `Point{D, T}`, `Segment{D, T}`, `Line{D, T}`, and `BBox{D, T}`.
- **Point Sequences**: `PntSeq{D, T}` and the matrix-backed `MatPntSeq{D, T}`.
- **Tree Structures**: Hierarchical `BBT.Tree` and `BBT.Node`.
- **Metric Spaces**: `AbsFMS` interface and `PointsSpace` wrappers.
- **Utilities**: `VArray` for efficient permutation and original index tracking.

## Metric Spaces

The library provides a robust framework for working with **Finite Metric Spaces (FMS)**, allowing algorithms to operate on abstract indices while the underlying representation handles the distance calculations.

### The `AbsFMS` Interface
All metric spaces subtype `AbsFMS`. They represent a set of $n$ points identified by indices $1 \dots n$. The core interface requires:
- `size(space)`: Number of points.
- `dist(space, i, j)`: Distance between points at indices $i$ and $j$.

### Provided Metric Space Types
- **`PointsSpace{T}`**: A simple wrapper around a `Vector` of any objects that support the `dist(p1, p2)` function.
- **`MPointsSpace{T}`**: A matrix-backed space where each column is treated as a point in $\mathbb{R}^d$. It uses optimized Euclidean distance calculations.
- **`AbsPntSeq`**: Any point sequence (like `PntSeq` or `MatPntSeq`) automatically implements the `AbsFMS` interface.
- **`PermutMetric`**: A powerful decorator that creates a "view" of an existing metric space under a specific permutation or subset of indices.
- **`SpherePMetric`**: A specialized space where the distance between two points $i$ and $j$ is the **angle** they form at a fixed base point $b$.

### Metric Space Algorithms
- **Greedy Permutation**: Generates a furthest-point ordering of the metric space.
    - `greedy_permutation_naive`: Standard $O(n^2)$ implementation of Gonzalez's algorithm.
    - `greedy_permutation_vanity`: A variation that breaks distance ties using a secondary "vanity" score.

## Quick Start

```julia
using BasicCompGeometry

# Create points
p1 = point(0.0, 0.0)
p2 = point(3.0, 4.0)

# Calculate distance
d = dist(p1, p2) # 5.0

# Work with segments
seg = Segment(p1, p2)
len = geom_length(seg)

# Bounding boxes
bb = BBox(p1, p2)
is_inside(point(1.0, 1.0), bb) # true

# Point Sequences (Polygons)
ps = PntSeq([p1, p2, point(0.0, 4.0)])
cardin(ps) # 3
```

## API Reference

### Points
```@docs
Point
point
dist
dist_sq
convex_comb
is_left_turn
is_right_turn
rand_point
rand_gaussian
```

### Segments & Lines
```@docs
Segment
Line
geom_length
at
nn_point
dist_point_segment
dist_segment_segment
bisection_point
```

### Point Sequences
```@docs
AbsPntSeq
PntSeq
MatPntSeq
Points
cardin
prefix_lengths
simplify
sample_uniformly
rand_pnt_seq
write_plt
```

### Bounding Boxes
```@docs
BBox
width
height
middle
diam
max_dist
is_inside
expand!
expand_add!
```

### Nearest Neighbor Search (BBT)
```@docs
BBT.exact_naive_scan
BBT.approx_nn
BBT.silly_nn
BBT.hybrid_nn
```

### Minimum Volume Bounding Box (MVBB)
```@docs
MVBB.approx_mvbb
MVBB.approx_diam
MVBB.OBBox
MVBB.volume
```

### Transformations (2D)
```@docs
translation
rotation
apply_transform
```

### Algorithms
```@docs
hausdorff_simplify
hausdorff_dist_subseg
distance_infty
centroid
match_price
```
