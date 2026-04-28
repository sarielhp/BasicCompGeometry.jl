# BasicCompGeometry.jl

*A comprehensive library for basic computational geometry in Julia.*

`BasicCompGeometry.jl` provides efficient, high-dimensional primitives and algorithms for computational geometry, leveraging Julia's multiple dispatch and the `StaticArrays.jl` ecosystem.

A core concept in the library is the `AbsPntSeq`, which is treated as a **sequence of points**. This can represent a simple point set, a polygonal chain, or a classical closed polygon. `AbsPolygon` is provided as an alias for backward compatibility.

## Features

- **High-Dimensional Primitives**: Support for Points, Segments, Lines, Point Sequences (PntSeq), and Axis-Aligned Bounding Boxes (BBox) in any dimension.
- **Multiple Coordinate Types**: Works seamlessly with `Float64`, `Int64`, and other numeric types.
- **Geometric Predicates**: Fast checks for turns (left/right) and containment.
- **Curve Algorithms**: Hausdorff distance-based simplification and uniform resampling.
- **2D Transformations**: Efficient translation and rotation for planar geometry.

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
