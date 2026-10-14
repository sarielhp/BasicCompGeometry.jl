# Point sequences

See [Geometry and point sequences](../geometry.md) for representations,
zero-copy matrix storage, and curve operations. Helpers such as `Points`,
`prefix_lengths`, and `write_plt` belong to [Extras](../imports.md).

```@autodocs
Modules = [BasicCompGeometry]
Pages = ["PntSeqs.jl", "PntSeqHausdorff.jl"]
Filter = obj -> obj in Main.documented_api || obj === Matrix
```

## Matching and collection helpers

```@autodocs
Modules = [BasicCompGeometry]
Pages = ["BasicCompGeometry.jl"]
Filter = obj -> obj in (BasicCompGeometry.centroid, BasicCompGeometry.distance_infty, BasicCompGeometry.match_price, BasicCompGeometry.VecPoint2I, BasicCompGeometry.VecPntSeq2F, BasicCompGeometry.VecPolygon2F)
```
