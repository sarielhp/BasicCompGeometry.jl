# Data files

`ReadWrite` reads and writes matrix datasets in binary formats, and reads
planar and higher-dimensional curves from text files (`read_plt`, the inverse of
`write_plt`, and `read_geolife_plt` for GeoLife trajectories).
Use `import BasicCompGeometry as CG` and qualified calls such as
`CG.ReadWrite.read_fbin(path)`, or import the namespace through
[Extras](../imports.md).

```@autodocs
Modules = [BasicCompGeometry.ReadWrite]
Private = false
```
