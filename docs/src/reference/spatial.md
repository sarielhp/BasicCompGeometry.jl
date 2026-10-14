# Spatial structures

See [Algorithms and spatial structures](../algorithms.md) for examples. These
operations belong to the `BBT`, `WSPD`, and `MVBB` submodules.

## Bounding box trees

| Call | Purpose |
|:-----|:--------|
| `BBT.Tree_init(points)` | Create a tree over an `AbsPntSeq` and initially expand its root |
| `BBT.Tree_fully_expand(tree)` | Refine the entire tree |
| `BBT.Tree_refine_node(tree, node)` | Split a node when possible |
| `BBT.get_min_orig_index(tree, node)` | Minimum input index stored under a node |
| `BBT.get_max_orig_index(tree, node)` | Maximum input index stored under a node |

```@autodocs
Modules = [BasicCompGeometry.BBT]
Private = false
```

## Well-separated pairs

```@autodocs
Modules = [BasicCompGeometry.WSPD]
Private = false
```

The exported helpers without separate docstrings are:

| Helper | Purpose |
|:-------|:--------|
| `WSPD.get_reps(state, pair)` | Return one representative point from each node |
| `WSPD.reps_orig_indexes(state, pair)` | Return the representatives' input indices |
| `WSPD.min_orig_index(state, pair)` | Minimum input index in either node |
| `WSPD.max_orig_index(state, pair)` | Maximum input index in either node |
| `WSPD.get_orig_ranges(state, pair)` | Return the bounding input-index range for each node |

These bounding ranges need not enumerate exactly the members of each node.

## Oriented bounding boxes

`MVBB.volume(box)` returns the volume of an `OBBox`.
`MVBB.rotating_calipers_min_area(hull)` takes a 2D convex hull and returns a
rectangle description with `vertices` (four corners) and `area` fields.

```@autodocs
Modules = [BasicCompGeometry.MVBB]
Private = false
```
