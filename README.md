# BasicCompGeometry.jl

[![Build Status](https://github.com/sarielhp/BasicCompGeometry.jl/workflows/Documentation/badge.svg)](https://github.com/sarielhp/BasicCompGeometry.jl/actions)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://sarielhp.github.io/BasicCompGeometry.jl/stable)

`BasicCompGeometry.jl` is a comprehensive Julia library for basic computational geometry. It provides efficient, high-dimensional primitives and algorithms, designed to be performant and easy to use through Julia's multiple dispatch system.

## Overview

This library was originally developed as a core component of the `FrechetDist` package. It has been refactored into a standalone module to provide a flat, idiomatic hierarchy for geometric types and algorithms that are useful across various computational geometry tasks.

## Features

- **Multi-Dimensional Primitives**: Support for Points, Segments, Lines, Polygons, and Axis-Aligned Bounding Boxes in any dimension (2D, 3D, and high-D).
- **Coordinate Agnostic**: Works with `Float64`, `Int64`, and other numeric types.
- **Fast Predicates**: Optimized checks for left/right turns, point-in-box containment, etc.
- **Distance Metrics**: Generic `dist` function for point-point, point-segment, segment-segment, and box-box distances.
- **Curve Algorithms**: Hausdorff distance-based simplification and uniform resampling of polygonal curves.
- **Planar Geometry**: Homogeneous 2D transformations (translation and rotation).

## Quick Start

```julia
using BasicCompGeometry

# Create a 2D point
p1 = npoint(0.0, 0.0)
p2 = npoint(3.0, 4.0)

# Euclidean distance
println(dist(p1, p2)) # 5.0

# Bounding Box containment
bb = BBox(p1, p2)
is_inside(npoint(1.5, 2.0), bb) # true

# Hausdorff Simplification
poly = rand_polygon(2, Float64, 1000)
simplified, indices = hausdorff_simplify(poly, 0.01)
```

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
