# BasicCompGeometry.jl

[![Build Status](https://github.com/sarielhp/BasicCompGeometry.jl/workflows/Documentation/badge.svg)](https://github.com/sarielhp/BasicCompGeometry.jl/actions)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://sarielhp.github.io/BasicCompGeometry.jl/dev)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

`BasicCompGeometry.jl` is a collection of computational-geometry primitives,
spatial data structures, and algorithms for Julia. The core package keeps a
small dependency set; Cairo and LaTeX support load through package extensions.

## Status and scope

This package is suitable for research code, experiments, and figure generation.
The WSPD implementation is actively supported and tested for separation and
exact pair coverage. Fréchet distance is not implemented in this package.

Most predicates use ordinary floating-point arithmetic rather than exact or
adaptive predicates. Degenerate or ill-conditioned inputs can therefore produce
incorrect combinatorial results, especially in convex hull and arrangement code.
Use a robust geometry library when correctness on adversarial numerical inputs is
required.


## Overview

This library provides a flat, idiomatic hierarchy for geometric types and algorithms that are useful across various computational geometry tasks.

## Features

- **Multi-Dimensional Primitives**: Support for Points, Segments, Lines, Point Sequences (PntSeq), and Axis-Aligned Bounding Boxes in any dimension (2D, 3D, and high-D).
- **Zero-Copy Matrix Integration**: Use `MatPntSeq` to treat columns of a matrix as points without copying memory.
- **Coordinate Agnostic**: Works with `Float64`, `Int64`, and other numeric types.
- **Geometric Predicates**: Checks for left/right turns, collinearity, and containment.
- **Distance Metrics**: Generic `dist` function for point-point, point-segment, segment-segment, and box-box distances.
- **Curve Algorithms**: Hausdorff distance-based simplification and uniform resampling of polygonal curves.
- **Planar Geometry**: Homogeneous 2D transformations (translation and rotation).
- **Vector Figure Generation (`IpeDraw`)**: Programmatic generation of publication-ready Ipe 7 XML figures, cropped PDFs, and LaTeX fragments with native math formulas and geometric dispatches.

## Algorithms

The library implements a variety of classic and modern geometric algorithms:

- **Convex Hull**:
    - **2D**: Monotone Chain algorithm ($O(n \log n)$).
    - **3D**: Gift-wrapping based implementation.
- **Diameter**:
    - **Exact**: Brute-force $O(n^2)$ calculation (`exact_diameter`, `exact_diameter_subspace`).
    - **Approximate**: $(1+\epsilon)$-approximation using WSPD ($O(n \log n + n/\epsilon^d)$) (`approx_diameter`, `approx_diameter_subspace`).
- **Nearest Neighbor Search**:
    - **Exact**: Brute-force $O(n)$ linear scan (`exact_naive_scan`).
    - **Approximate**: $c$-approximation using BBT ($O(\log n)$ for many distributions) (`approx_nn`).
    - **Silly**: Simple tree descent for extremely fast, heuristic results (`silly_nn`).
    - **Hybrid**: Combines Silly and Approximate for better performance (`hybrid_nn`).
- **Spatial Decomposition**:
    - **WSPD**: Well-Separated Pairs Decomposition.
    - **BBT**: Bounding Box Tree construction and expansion.
    - **MVBB**: Approximate Minimum Volume Bounding Box ($O(n \log n + n \text{ poly}(1/\epsilon))$).
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

`BasicCompGeometry` provides typed geometric data structures:

- **Geometric Primitives**:
    - `Point{D, T}`: High-dimensional points (using `StaticArrays`).
    - `Segment{D, T}`: Directed line segments.
    - `Line{D, T}`: Infinite lines defined by a point and a direction.
    - `BBox{D, T}`: Axis-aligned bounding boxes.
- **Curves & Splines**:
    - `Ellipse{T}` & `EllipticArc{T}`: Exact bounding boxes, containment, area, and uniform sampling.
    - `CubicBezier{D, T}`: de Casteljau evaluation, exact BBox, subdivision, and adaptive flattening.
    - `CubicSpline{D, T}`: Composite spline with Catmull-Rom ($C^1$) and Natural Spline ($C^2$) interpolation.
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

# Nearest Neighbor Search (using Bounding Box Tree)
tree = BBT.Tree_init(mp)
BBT.Tree_fully_expand(tree)
q = point(0.5, 0.5)

# Exact, Approximate (c=1.5), and Heuristic ("Silly") NN
d1, p1, i1 = BBT.exact_naive_scan(tree, q)
d2, p2, i2 = BBT.approx_nn(tree, q, 1.5)
d3, p3, i3 = BBT.silly_nn(tree, q)
```

## AbsPntSeq Interface

The library provides an `AbsPntSeq{D, T}` abstract interface representing a **sequence of points** (which can be viewed as a point set, a polygonal chain, or a classical polygon). All geometric algorithms (BBT, WSPD, Diameter, etc.) are implemented against this interface, allowing them to work seamlessly with different storage backends:
- `PntSeq{D, T}`: Standard `Vector{Point{D, T}}` backed representation.
- `MatPntSeq{D, T}`: Matrix-backed representation (`D x N` matrix) for zero-copy integration with existing datasets.

`AbsPolygon`, `Polygon`, and `MatPolygon` are provided as aliases for backward compatibility.

The interfaces preserve concrete coordinate and storage types so Julia can specialize
the corresponding algorithms.

## Documentation

For detailed information on all types and functions, please see the [Latest Documentation](https://sarielhp.github.io/BasicCompGeometry.jl/dev).

## Installation

Until the package is registered, install it directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/sarielhp/BasicCompGeometry.jl")
```

## Examples

Ready-to-run example scripts are documented in
[`examples/README.md`](examples/README.md).

## Vector Figure Generation (`IpeDraw`)

The `IpeDraw` submodule provides programmatic generation of publication-ready vector figures using the [Ipe extensible drawing editor](http://ipe.otfried.org/) format (`.ipe`).

The compact API dispatches on geometry types, fits world coordinates to the page,
keeps label offsets and marks in page units, and renders publication PDFs through
Ipe so that geometry and LaTeX labels remain vector content:

```julia
using BasicCompGeometry
using BasicCompGeometry.IpeDraw

c1 = Circle(point(0.0, 0.0), 1.0)
c2 = Circle(point(1.0, 0.0), 1.0)
p, q = sort(intersections(c1, c2); by = p -> -p.y)

figure("output/disks.pdf"; fit=[c1, c2], margin=12) do fig
    disk = Style(fill_opacity=0.2, pen=:heavier)
    draw!(fig, c1; style=disk, stroke=:blue, fill=:lightblue)
    draw!(fig, c2; style=disk, stroke=:darkgreen, fill=:lightgreen)
    mark!(fig, [p, q])
    label!(fig, p, raw"p"; offset=(6, 4), anchor=:southwest)
    label!(fig, q, raw"\sqrt{q}"; offset=(6, -6), anchor=:northwest)
end
```

This writes `output/disks.pdf` and retains the editable
`output/disks.ipe` source. Pass `keep_source=false` when only the PDF is
wanted, `preview=true` to open the generated source in Ipe, or target a `.ipe`
path to skip PDF compilation. PDF targets require
Ipe's `ipetoipe` command on `PATH`.

Page-space helpers make common figure furniture independent of the fitted world
coordinates. Named themes keep repeated styles together, while scoped insets
and clipping preserve the surrounding viewport:

```julia
theme = Theme(region=Style(fill=:lightblue, stroke=:blue, fill_opacity=0.2))

inset(fig, PageBox(400, 320, 150, 140); fit=detail_bbox) do zoom
    clip_to(zoom, detail_bbox) do clipped
        draw!(clipped, circles; style=theme.region)
    end
end
legend!(fig, ["feasible" => theme.region]; position=:northwest)
scale_bar!(fig, 0.5; label=raw"1/2")
```

- **Native Geometric Dispatches**: Direct methods for `Point`, `Segment`, `BBox`, `Circle`, `CircleArc`, `Ellipse`, `EllipticArc`, `CubicBezier`, and `CubicSpline`.
- **Advanced Primitives**: Smooth splines (`draw_spline!`), approximating B-splines (`draw_bspline!`), polygons with holes (`draw_polygon_with_holes!`), and scoped groups with affine transforms (`ipe_group`).
- **Conceptual & Algorithmic Helpers**: Primitives for intervals, spans, dimension lines, and curved arrows (`draw_bar!`, `draw_span!`, `draw_dimension!`, `draw_curved_arrow!`).
- **LaTeX Math Formulas**: Native math support (accepts `LaTeXStrings` `L"..."` and standard LaTeX strings) with automatic XML escaping.
- **Hatch Patterns & Opacities**: Vector hatch patterns (`:hatch`, `:crosshatch`, `:vertical`, `:horizontal`, `:falling`, `:rising`) and 10%–90% opacity fills.
- **Multi-Layer & Multi-View**: Easily define layers and progressive presentation views.
- **Figure Composition**: Scoped inset viewports, editable clipping groups, named themes, legends, and scale bars.
- **Self-Contained Style**: Bundled 8-inch canvas with auto-crop (`crop="yes"`, `bbox="cropbox"`), extended pens, and rich academic color palettes.
- **Automated Compilation**: Generates editable `.ipe`, cropped vector `.pdf`, and companion `_fig.tex` LaTeX wrappers.

### Quick Example

```julia
using BasicCompGeometry
using BasicCompGeometry.IpeDraw
using LaTeXStrings

open_ipe("output/interval_demo"; caption="Interval estimation", label="fig:demo") do cv
    add_layer!(cv, "intervals")
    add_layer!(cv, "labels")

    set_layer!(cv, "intervals")
    draw_bar!(cv, 50.0, 350.0, 100.0; label_left=L"0", label_right=L"1")
    draw_span!(cv, 100.0, 250.0, 100.0; fill=:lightgreen, stroke=:darkgreen)
    draw_dimension!(cv, 100.0, 250.0, 130.0; label=L"\Delta \le \varepsilon", arrow=:both)

    set_layer!(cv, "labels")
    draw_point!(cv, 175.0, 100.0; stroke=:darkblue, fill=:darkblue, shape=:disk)
    draw_label!(cv, 175.0, 85.0, L"\mu"; halign=:center)
end
```
See [`examples/ipe_conceptual_figure.jl`](examples/ipe_conceptual_figure.jl) for a complete example.

## Visualization & Optional Dependencies

Starting with Julia 1.9, `BasicCompGeometry` uses **Package Extensions** to keep the core library lightweight. Visualization features (like `BBT.Tree_draw`) are only available when you explicitly load the following packages in your environment:

- **Cairo.jl**
- **Colors.jl**

```julia
using BasicCompGeometry
using Cairo, Colors # This triggers the BBTCairoExt extension

# Now Tree_draw is available
BBT.Tree_draw(tree, "output/tree.pdf")
```

## Origins

The code in this module was originally part of the `FrechetDist`
package, but this package does not currently provide Fréchet-distance algorithms.
Sariel Har-Peled wrote the original geometry code. Gemini CLI and OpenAI Codex
have assisted substantially with later implementation, refactoring, tests,
documentation, and visualization interfaces.

The package has since been reorganized around:
- A flat module hierarchy.
- Type-generic implementations.
- Integration with the `StaticArrays.jl` ecosystem.
- Full compatibility with standard `Base` methods through multiple dispatch.

## Maintenance and license

Bug reports with a small reproducing example are welcome. Maintenance is
best-effort, with correctness bugs prioritized. The package is distributed under
the [MIT License](LICENSE).
