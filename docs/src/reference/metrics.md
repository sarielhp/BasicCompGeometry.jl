# Metric spaces

See [Finite metric spaces](../metrics.md) for the interface and examples.
`MetricSpace` exports its space types and indexed-distance operations;
`VirtArray` is an auxiliary namespace available through [Extras](../imports.md).

```@docs
BasicCompGeometry.AbsFMS
```

```@autodocs
Modules = [BasicCompGeometry.MetricSpace]
Filter = obj -> obj in Main.documented_api || obj === size
```

## Hilbert spaces and virtual points

See [Hilbert spaces and virtual points](../hilbert.md) for the guide.

```@docs
BasicCompGeometry.AbsVirtualPoint
HilbertSpace
VPoint
vpoint
vcurve
BasicCompGeometry.local_frame
BasicCompGeometry.is_generator
BasicCompGeometry.is_vertex_curve
BasicCompGeometry.same_space
BasicCompGeometry.cache_stats
Base.empty!(::HilbertSpace)
BasicCompGeometry.vdist2
BasicCompGeometry.line_params
BasicCompGeometry.inner_diff
BasicCompGeometry.VDiff
BasicCompGeometry.vscale
BasicCompGeometry.dist_sq(::VPoint, ::VPoint)
BasicCompGeometry.convex_comb(::VPoint, ::VPoint, ::Real)
```

## Virtual arrays

`VirtArray.VArray` wraps a vector and an index permutation without reordering
the original storage. `VirtArray.orig_index(array, i)` returns the underlying
vector index at position `i`; `VirtArray.swap!(array, i, j)` exchanges virtual
positions.

```@autodocs
Modules = [BasicCompGeometry.VirtArray]
Private = false
```
