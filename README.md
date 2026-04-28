# BasicCompGeometry.jl

[![Build Status](https://github.com/sarielhp/BasicCompGeometry.jl/workflows/Documentation/badge.svg)](https://github.com/sarielhp/BasicCompGeometry.jl/actions)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://sarielhp.github.io/BasicCompGeometry.jl/dev)

`BasicCompGeometry.jl` is a Julia library for basic computational
geometry stuff. It contains some code I wrote over the last few years
in a self contained cleaner form.  Its a bit of a random collection of
things so far, but hopefully would become more coherent over time. The
main guideline is to have self-contained code that does not rely on
heavy external packages (some of the examples use heavier packages for
visualization).

## Disclaimer

Much of the code in this package (like the convex hull code in 2 and 3
dimensions) was written using gemini-cli. Writing robust and efficient
code for Computational Geometry is notoriously difficult - see
[CGAL](https://www.cgal.org/) for a library that does it right!
(Unlike this library!) Some things are naturally robust and should
work fine, but even convex-hull computation in 2d using floating point
can become tricky. 



## Overview

This library provides a flat, idiomatic hierarchy for geometric types and algorithms that are useful across various computational geometry tasks.

## Features

- **Multi-Dimensional Primitives**: Support for Points, Segments, Lines, Point Sequences (PntSeq), and Axis-Aligned Bounding Boxes in any dimension (2D, 3D, and high-D).
- **Zero-Copy Matrix Integration**: Use `MatPntSeq` to treat columns of a matrix as points without copying memory.
- **Coordinate Agnostic**: Works with `Float64`, `Int64`, and other numeric types.
- **Fast Predicates**: Optimized checks for left/right turns, point-in-box containment, etc.
- **Distance Metrics**: Generic `dist` function for point-point, point-segment, segment-segment, and box-box distances.
- **Curve Algorithms**: Hausdorff distance-based simplification and uniform resampling of polygonal curves.
- **Planar Geometry**: Homogeneous 2D transformations (translation and rotation).

## Algorithms

The library implements a variety of classic and modern geometric algorithms:

- **Convex Hull**:
    - **2D**: Monotone Chain algorithm ($O(n \log n)$).
    - **3D**: Gift-wrapping based implementation.
- **Diameter**:
    - **Exact**: Brute-force $O(n^2)$ calculation.
    - **Approximate**: $(1+\epsilon)$-approximation using WSPD ($O(n \log n + n/\epsilon^d)$).
- **Spatial Decomposition**:
    - **WSPD**: Well-Separated Pairs Decomposition.
    - **BBT**: Bounding Box Tree construction and expansion.
- **Curve Processing**:
    - **Simplification**: Hausdorff distance-based simplification.
    - **Resampling**: Uniform arc-length resampling.
- **Metric Space Algorithms**:
    - **Greedy Permutation**: Incremental furthest-point sampling ($O(n^2)$ or optimized via BBT).
- **Advanced Optimization**:
    - **Longest Convex Subset**: Dynamic programming based calculation.
- **Geometric Predicates**:
    - Fast `turn_sign`, `is_left_turn`, `is_right_turn`, and `is_collinear` checks.
- **Distance Calculations**:
    - Point-Point, Point-Segment, Point-Line, Segment-Segment.
    - Min/Max distance between Bounding Boxes.

## Data Structures

`BasicCompGeometry` provides highly optimized, type-safe data structures:

- **Geometric Primitives**:
    - `Point{D, T}`: High-dimensional points (using `StaticArrays`).
    - `Segment{D, T}`: Directed line segments.
    - `Line{D, T}`: Infinite lines defined by a point and a direction.
    - `BBox{D, T}`: Axis-aligned bounding boxes.
- **Point Sequences (Polygons)**:
    - `PntSeq{D, T}`: Standard sequence of points.
    - `MatPntSeq{D, T}`: Zero-copy view of a Julia `Matrix` as a point sequence.
- **Tree Structures**:
    - `BBT.Tree`: Hierarchical Bounding Box Tree.
    - `BBT.Node`: Internal nodes and leaves of the BBT.
- **Metric Spaces**:
    - `AbsFMS`: Abstract interface for Finite Metric Spaces.
    - `PointsSpace`: Wrapper for treated point sets as metric spaces.
- **Utilities**:
    - `VArray`: Virtual array for efficient permutation and original index tracking.


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
- **`PermutMetric`**: A powerful decorator that creates a "view" of an existing metric space under a specific permutation or subset of indices. It includes a `swap!(space, i, j)` method to efficiently reorder the virtual indices.
- **`SpherePMetric`**: A specialized space where the distance between two points $i$ and $j$ is the **angle** they form at a fixed base point $b$ (i.e., the angular separation $\angle ibj$).

### Metric Space Algorithms
- **Greedy Permutation**: Generates a furthest-point ordering of the metric space. This is often used for $k$-center clustering or creating hierarchical approximations.
    - `greedy_permutation_naive`: Standard $O(n^2)$ implementation of Gonzalez's algorithm.
    - `greedy_permutation_vanity`: A variation that breaks distance ties using a secondary "vanity" score.

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
ps = rand_pnt_seq(2, Float64, 1000)
simplified, indices = hausdorff_simplify(ps, 0.01)

# Zero-copy matrix wrapper
M = rand(2, 500)
mp = MatPntSeq(M)
d = exact_diameter(mp)
```

## AbsPntSeq Interface

The library provides an `AbsPntSeq{D, T}` abstract interface representing a **sequence of points** (which can be viewed as a point set, a polygonal chain, or a classical polygon). All geometric algorithms (BBT, WSPD, Diameter, etc.) are implemented against this interface, allowing them to work seamlessly with different storage backends:
- `PntSeq{D, T}`: Standard `Vector{Point{D, T}}` backed representation.
- `MatPntSeq{D, T}`: Matrix-backed representation (`D x N` matrix) for zero-copy integration with existing datasets.

`AbsPolygon`, `Polygon`, and `MatPolygon` are provided as aliases for backward compatibility.

By using Julia's parametric type system, these abstractions incur **zero runtime overhead**.

## Documentation

For detailed information on all types and functions, please see the [Latest Documentation](https://sarielhp.github.io/BasicCompGeometry.jl/dev).

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/sarielhp/BasicCompGeometry.jl")
```

## Origins

The code in this module was originally part of the `FrechetDist`
package. It has been modernized and reorganized to follow Julia's best
practices, including:
- A flat module hierarchy.
- Type-generic implementations.
- Integration with the `StaticArrays.jl` ecosystem.
- Full compatibility with standard `Base` methods through multiple dispatch.
