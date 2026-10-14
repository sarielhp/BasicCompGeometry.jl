# Curves and conics

Bézier curves, composite splines, hyperbolas, parabolas, and curve polygons.
Helpers such as `derivative` and specialized conic accessors belong to
[Extras](../imports.md).

`Curve2DF` is the `Float64` specialization of `Curve2D`.
`subdivide(bezier, t=0.5)` returns the two Bézier curves obtained by splitting
at parameter `t`; it is equivalent to `split(bezier, t)`.

```@docs
BasicCompGeometry.AbsCurve2D
```

```@autodocs
Modules = [BasicCompGeometry]
Pages = ["Bezier.jl", "Hyperbolas.jl", "Parabolas.jl", "CurvePolygon.jl"]
Filter = obj -> obj in Main.documented_api || obj === split
```
