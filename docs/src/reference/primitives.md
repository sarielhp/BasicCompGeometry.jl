# Primitives and transformations

See [Geometry and point sequences](../geometry.md) for usage. Generic helpers
such as `distance`, `width`, `height`, and `at` require qualification, explicit
imports, or [Extras](../imports.md).

```@autodocs
Modules = [BasicCompGeometry]
Pages = ["Points.jl", "Segments.jl", "Planes.jl", "BBoxes.jl", "Transforms2D.jl"]
Filter = obj -> obj in Main.documented_api
```

```@docs
BasicCompGeometry.distance
```
