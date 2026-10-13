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


## Highlights

- **AI-Assisted Figure Generation with `IpeDraw`**:
  A powerful workflow for researchers is prompting AI coding agents (such as Google Gemini, Claude, or ChatGPT) to write self-contained Julia scripts that compute geometric algorithms and produce publication-ready vector figures (`.ipe` and `.pdf`). With the native `IpeDraw` submodule, output figures retain editable vector layers, native LaTeX formulas, automatic page fitting, custom styles/themes, legends, scale bars, and inset zoom viewports.
- **Zero-Ceremony Standalone Scripts with [`QuickEnv.jl`](https://github.com/sarielhp/QuickEnv.jl)**:
  Instead of managing virtual environments or remembering `--project` flags, standalone figure scripts can simply start with `using QuickEnv` followed by `using BasicCompGeometry, BasicCompGeometry.IpeDraw`. QuickEnv automatically isolates dependencies without polluting the global environment, allowing scripts to be executed directly via `./script.jl` or `julia script.jl`.
- **Live Tested AI Prompt**:
  See below for a live test prompt and the resulting standalone script in [`examples/convex_hull_figure.jl`](examples/convex_hull_figure.jl), generated end-to-end by an AI agent (Gemini Flash).
- **Exact & Pluggable Predicate Kernels**:
  Easily opt into exact arithmetic decisions (`kernel=ExactKernel()`) for robust geometric decisions (arbitrary-precision integer/rational coordinates and algebraic radical isolation) without changing point storage types or incurring overhead on default floating-point calls.

### AI-Assisted Figure Generation: Live Test

To demonstrate how AI agents can leverage `BasicCompGeometry.jl` and `QuickEnv.jl` to generate figures on demand, here is an actual prompt tested with **Gemini Flash**:

> **Prompt:**  
> *"Using Julia, QuickEnv, and BasicCompGeometry (with IpeDraw), write a standalone executable script that generates a publication-quality Ipe/PDF figure illustrating the convex hull and diameter of a 2D point set:*  
> *1. Generate 35 deterministic 2D points inside a $[0, 10]^2$ box.*  
> *2. Compute their 2D convex hull with `convex_hull`.*  
> *3. Find the diameter pair $(p, q)$ among the hull vertices and the diameter segment.*  
> *4. Draw the bounding box, points (as small gray dots), the convex hull (as a translucent light-blue polygon with a navy border), hull vertices, and the diameter segment (as a dashed red line).*  
> *5. Add LaTeX labels: label vertex $p$ with $p$, vertex $q$ with $q$, and label the midpoint of the diameter segment with $D = \mathrm{diam}(P)$.*  
> *6. Add a legend in the top-left corner with LaTeX math labels ($\mathcal{CH}(P)$ and $D$), a scale bar of length $2.0$, and an inset zoom box in the top-right corner focusing on vertex $p$.*  
> *7. Save the output to `output/convex_hull_figure.pdf` and keep the editable `.ipe` file."*

The complete, working script produced from this prompt is available at [`examples/convex_hull_figure.jl`](examples/convex_hull_figure.jl). Running it directly:
```bash
./examples/convex_hull_figure.jl
```
runs completely standalone via `QuickEnv` and produces both the editable `output/convex_hull_figure.ipe` and the compiled vector `output/convex_hull_figure.pdf`.

## Overview

### Optional predicate kernels

Ordinary calls need no setup and use floating-point predicates. Opt into exact
arithmetic only where needed:

```julia
convex_hull(points)                       # floating point, as usual
convex_hull(points; kernel=ExactKernel()) # exact decisions, same point types
is_left_turn(p, q, r; kernel=ExactKernel())
```

The same kernel is passed through every predicate in a kernel-aware algorithm.
`FloatingKernel()` uses `Float64`; `FloatingKernel(Float32)` is also available.
`ExactKernel()` uses arbitrary-precision integers and rationals and compares
circular radicals algebraically. Finite floating-point inputs are interpreted
as their stored binary values. No global setting, context object, or new point
type is required. Existing integer predicate results now use floating-point
arithmetic by default; request the exact kernel when integer exactness matters.

Currently kernel-aware: 2D turn predicates, lexicographic point comparison,
distance-threshold comparisons, sphere containment, right-circle-boundary
comparisons, and the 2D convex hull. Other algorithms and coordinate
constructions have not been converted. An exact predicate cannot repair an
already rounded intersection. `BigFloat` alone is not an exact kernel.

Custom kernels subtype `AbstractKernel` and implement kernel-first predicate
methods such as `turn_sign(kernel, p, q, r)` and `compare_xy(kernel, p, q)`.
Missing methods do not fall back silently. A future certified filtered kernel
can use this same interface; no filtered kernel is claimed in this release.

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

The basic library was originally written by Sariel Har-Peled (initially developed as part of the `FrechetDist` package, with large portions based on computational geometry algorithms and code originally written in C++ 25 years ago).

The package has since been reorganized around:
- A flat module hierarchy.
- Type-generic implementations.
- Integration with the `StaticArrays.jl` ecosystem.
- Full compatibility with standard `Base` methods through multiple dispatch.

## AI Contribution & Code Ownership

More recent versions, modernizations, extensions, tests, documentation, and predicate kernels were written and generated with the assistance of AI coding tools (including Google Gemini CLI and OpenAI Codex). 

As the author and maintainer, I fully understand the codebase, curate its mathematical and algorithmic design, take full responsibility for all code, and actively support maintenance, bug fixes, and future development.

## Maintenance and license

Bug reports with a small reproducing example are welcome. Maintenance is
best-effort, with correctness bugs prioritized. The package is distributed under
the [MIT License](LICENSE).
