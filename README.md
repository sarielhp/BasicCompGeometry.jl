# BasicCompGeometry.jl

[![Build Status](https://github.com/sarielhp/BasicCompGeometry.jl/workflows/Documentation/badge.svg)](https://github.com/sarielhp/BasicCompGeometry.jl/actions)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://sarielhp.github.io/BasicCompGeometry.jl/stable)

`BasicCompGeometry.jl` is a comprehensive Julia library for basic computational geometry. It provides efficient, high-dimensional primitives and algorithms, designed to be performant and easy to use through Julia's multiple dispatch system.

## Overview

This library provides a flat, idiomatic hierarchy for geometric types and algorithms that are useful across various computational geometry tasks.

## Features

- **Multi-Dimensional Primitives**: Support for Points, Segments, Lines, Polygons, and Axis-Aligned Bounding Boxes in any dimension (2D, 3D, and high-D).
- **Zero-Copy Matrix Integration**: Use `MatPolygon` to treat columns of a matrix as points without copying memory.
- **Coordinate Agnostic**: Works with `Float64`, `Int64`, and other numeric types.
- **Fast Predicates**: Optimized checks for left/right turns, point-in-box containment, etc.
- **Distance Metrics**: Generic `dist` function for point-point, point-segment, segment-segment, and box-box distances.
- **Curve Algorithms**: Hausdorff distance-based simplification and uniform resampling of polygonal curves.
- **Planar Geometry**: Homogeneous 2D transformations (translation and rotation).

## Quick Start

```julia
using BasicCompGeometry

# Create a 2D point
p1 = point(0.0, 0.0)
p2 = point(3.0, 4.0)

# Euclidean distance
println(dist(p1, p2)) # 5.0

# Bounding Box containment
bb = BBox(p1, p2)
is_inside(point(1.5, 2.0), bb) # true

# Hausdorff Simplification
poly = rand_polygon(2, Float64, 1000)
simplified, indices = hausdorff_simplify(poly, 0.01)

# Zero-copy matrix wrapper
M = rand(2, 500)
mp = MatPolygon(M)
d = exact_diameter(mp)
```

## AbsPolygon Interface

The library provides an `AbsPolygon{D, T}` abstract interface. All geometric algorithms (BBT, WSPD, Diameter, etc.) are implemented against this interface, allowing them to work seamlessly with different storage backends:
- `Polygon{D, T}`: Standard `Vector{Point{D, T}}` backed representation.
- `MatPolygon{D, T}`: Matrix-backed representation (`D x N` matrix) for zero-copy integration with existing datasets.

By using Julia's parametric type system, these abstractions incur **zero runtime overhead**.

## Documentation

For detailed information on all types and functions, please see the [Latest Documentation](https://sarielhp.github.io/BasicCompGeometry.jl/dev).

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/sarielhp/BasicCompGeometry.jl")
```

## Origins

The code in this module was originally part of the `FrechetDist` package. It has been modernized and reorganized to follow Julia's best practices, including:
- A flat module hierarchy.
- Type-generic implementations.
- Integration with the `StaticArrays.jl` ecosystem.
- Full compatibility with standard `Base` methods through multiple dispatch.
