# Circles and ellipses

Spheres, circles, circular arcs, ellipses, and elliptic arcs. Accessors such as
`center`, `r_major`, `r_minor`, and `rotation_angle` require qualification,
explicit imports, or [Extras](../imports.md).

| Accessor | Result |
|:---------|:-------|
| `r_major(ellipse_or_arc)` | Major semi-axis radius |
| `r_minor(ellipse_or_arc)` | Minor semi-axis radius |
| `alpha1(elliptic_arc)` | Starting ellipse parameter angle |
| `alpha2(elliptic_arc)` | Ending ellipse parameter angle |
| `end_point(elliptic_arc)` | Point at the ending angle; paired with `start_point` |

```@autodocs
Modules = [BasicCompGeometry]
Pages = ["Spheres.jl", "Ellipses.jl"]
Filter = obj -> obj in Main.documented_api
```
