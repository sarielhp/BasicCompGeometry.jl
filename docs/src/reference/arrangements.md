# Halfplanes and arrangements

Halfplane geometry, line bundles, and vertical decompositions. Helpers such as
`clip`, `complement`, `faces`, and `locate` require qualification, explicit
imports, or [Extras](../imports.md). See
[Kernels and numerical scope](../kernels.md) for floating-point limitations.

| Accessor | Result |
|:---------|:-------|
| `boundary_line(halfplane)` | Boundary line; equivalent to `boundary(halfplane)` |
| `trapezoids(decomposition)` | Vertical trapezoids in an `ArrVertDecomp` |
| `view_box(decomposition)` | Bounding box of the decomposition's view |
| `lines(decomposition)` | The decomposition's lines; available through `Extras` |

```@autodocs
Modules = [BasicCompGeometry]
Pages = ["Halfplanes.jl", "Bundles.jl", "Arrangements.jl"]
Filter = obj -> obj in Main.documented_api
```
