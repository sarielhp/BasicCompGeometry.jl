# Simplification hierarchies

`PntSeqHierarchy` builds palettes and hierarchies of simplifications of a
point sequence, controlled by a measure of subcurve error. The default measure
is the Hausdorff distance to the chord; supply a
[`HierarchyMeasure`](@ref BasicCompGeometry.PntSeqHierarchy.HierarchyMeasure)
to control a different error, such as a Frechet distance. The names are not
exported: use qualified calls such as
`BasicCompGeometry.PntSeqHierarchy.hierarchy_init(P)`.

```@autodocs
Modules = [BasicCompGeometry.PntSeqHierarchy]
Private = false
```
