# Getting started

## Installation

Requires Julia 1.9 or newer. Install the development version described by these
docs directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/sarielhp/BasicCompGeometry.jl")
```

For a release available in Julia's General registry, use
`Pkg.add("BasicCompGeometry")`. These development docs describe the v0.5.0 API;
older releases have a larger default export list.

For a local checkout, run `Pkg.develop(path="/path/to/BasicCompGeometry")` in
your working environment. Optional figure dependencies are described in
[Vector figures](figures.md) and [Cairo and LaTeX](cairo.md).

## A first calculation

```@example start
using BasicCompGeometry

p = point(0.0, 0.0)
q = point(3.0, 4.0)
segment = Segment(p, q)
box = BBox(p, q)

(distance=dist(p, q),
 segment_length=geom_length(segment),
 inside=is_inside(point(1.0, 1.0), box))
```

Point sequences supply the input for many algorithms:

```@example start
points = PntSeq([p, q, point(0.0, 4.0), point(1.0, 2.0)])
hull = convex_hull(points; kernel=ExactKernel())
(hull_vertices=length(hull), diameter=exact_diameter(points))
```

`ExactKernel()` makes the hull's predicate decisions exact on the supplied
coordinates. Other algorithms, including diameter calculations, retain their
own arithmetic; see [Kernels and numerical scope](kernels.md).

Use [Geometry and point sequences](geometry.md) for storage and curve operations,
or [Algorithms and spatial structures](algorithms.md) for spatial searches.
