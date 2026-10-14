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

## Virtual arrays

`VirtArray.VArray` wraps a vector and an index permutation without reordering
the original storage. `VirtArray.orig_index(array, i)` returns the underlying
vector index at position `i`; `VirtArray.swap!(array, i, j)` exchanges virtual
positions.

```@autodocs
Modules = [BasicCompGeometry.VirtArray]
Private = false
```
