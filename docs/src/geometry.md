# Geometry and point sequences

## Primitives

`Point{D,T}`, `Segment{D,T}`, `Line{D,T}`, and `BBox{D,T}` represent geometry
in dimension `D` with coordinate type `T`. Points use `StaticArrays`; aliases
such as `Point2F` and `Point3F` select common dimensions and coordinate types.
For 2D points, use `.x` and `.y` to access coordinates.

```@example geometry
using BasicCompGeometry

p = point(0.0, 0.0)
q = point(3.0, 4.0)
bb = BBox(p, q)
(x=q.x, y=q.y, distance=dist(p, q), inside=is_inside(point(1.5, 2.0), bb))
```

`dist` handles point, segment, and box distances. `nn_point` finds the nearest
point on a segment. `bbox` constructs bounds for supported objects, and
`union_bbox` combines boxes. General helpers such as `width`, `height`, `middle`,
and `area` belong to [Extras](imports.md).

The package also supplies planes, spheres and circles, circular and elliptic
arcs, conics, Bézier curves, splines, and curve polygons. Consult the reference
pages for [primitives](reference/primitives.md),
[circles and ellipses](reference/round-shapes.md), and
[curves and conics](reference/curves.md).

## Point sequences

`AbsPntSeq{D,T}` represents a point set, polygonal chain, or polygon. Algorithms
such as hull construction, diameter, BBT, and WSPD accept this interface.

| Representation | Storage |
|:---------------|:--------|
| `PntSeq{D,T}` | A `Vector{Point{D,T}}` |
| `MatPntSeq` | A zero-copy view of columns of a `D × N` matrix |

Sequences support indexing, iteration, `length`, and assignment. `PntSeq`
also supports `push!`, `pop!`, and `empty!`. `size(sequence)` returns the number
of points as an integer, following this package's finite metric space interface.
`dim(sequence)` returns the coordinate dimension.

```@example geometry
M = [0.0 3.0 0.0; 0.0 4.0 4.0]
points = MatPntSeq(M)
M[1, 2] = 6.0
(count=length(points), dimension=dim(points), second_point=points[2])
```

The matrix and sequence share storage; edits through either representation
affect the other. `Matrix(points)` creates a coordinate matrix copy.
`CG.Points(points)` returns the underlying point vector when using
`import BasicCompGeometry as CG`.

Compatibility aliases include `Polygon = PntSeq`, `CG.AbsPolygon = AbsPntSeq`,
and `CG.MatPolygon = MatPntSeq`. The type does not itself distinguish an open
chain from a closed polygon; check each operation's closure convention.

## Polygonal curve operations

`geom_length` sums consecutive edge lengths. `sample_uniformly` resamples an
open chain by arc length. `simplify` skips vertices close to the last kept
vertex; `hausdorff_simplify` uses the Hausdorff-based simplification algorithm.
Both simplifiers return the simplified sequence and the kept input indices.

`push_distinct!` appends a point unless it repeats the last vertex, `split_edges`
inserts the midpoint of every edge, and `reverse` reverses the vertex order. The last
two return new sequences and leave their input unchanged. `ReadWrite.read_plt` reads a point sequence from a
text file, and `ReadWrite.read_geolife_plt` reads a GeoLife trajectory.

```@example geometry
curve = PntSeq([point(0.0, 0.0), point(0.5, 0.0), point(1.0, 0.0)])
simplified, kept = hausdorff_simplify(curve, 0.01)
(length=geom_length(curve), kept=kept, samples=collect(sample_uniformly(curve, 5)))
```

`translate!` and `scale!` mutate coordinates. `translation`, `rotation`, and
`apply_transform` provide planar transformations. The
[point sequence reference](reference/sequences.md) also documents interpolation,
random generation, and compatibility helpers.
