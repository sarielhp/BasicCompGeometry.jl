# BasicCompGeometry.jl

Computational-geometry primitives, spatial data structures, algorithms, and
vector figures for Julia. Points, segments, and bounding boxes work in arbitrary
dimensions; matrix-backed point sequences let existing datasets share storage
with geometry algorithms. The core has a small dependency set, with Cairo and
LaTeX support provided through optional package extensions.

Start with [installation and a first example](getting-started.md), then choose a guide:

| Task | Guide |
|:-----|:------|
| Select imports or migrate to v0.5.0 | [Imports and Extras](imports.md) |
| Represent points, boxes, and polygonal curves | [Geometry and point sequences](geometry.md) |
| Choose predicate arithmetic | [Kernels and numerical scope](kernels.md) |
| Compute hulls, diameters, or nearest neighbors | [Algorithms and spatial structures](algorithms.md) |
| Work with indexed distances and greedy permutations | [Finite metric spaces](metrics.md) |
| Generate editable Ipe figures and vector PDFs | [Vector figures](figures.md) |
| Draw with Cairo or add LaTeX annotations | [Cairo and LaTeX](cairo.md) |

The API reference in the sidebar lists the default API, the supported `Extras`
helpers, and the public operations in feature submodules. See the
[example catalog](https://github.com/sarielhp/BasicCompGeometry.jl/blob/main/examples/README.md)
for ready-to-run programs.

The package targets research code, experiments, and figure generation. Most
algorithms use floating-point decisions; kernel-aware operations can use exact
predicates. Read [the numerical limits](kernels.md) before relying on degenerate
or ill-conditioned inputs. Fréchet distance is not implemented here.

[About the package](about.md) covers its origins, code ownership, maintenance,
and license.
