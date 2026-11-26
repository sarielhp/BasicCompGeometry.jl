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
| `PntSeq{D,T}` | A `Vector{Point{D,T}}`, stored contiguously |
| `GPntSeq{D,T,P}` | A `Vector{P}` of points of type `P`; `PntSeq{D,T}` is `GPntSeq{D,T,Point{D,T}}` |
| `MatPntSeq` | The columns of a `D × N` matrix, read as static `Point{D,T}`; shares the matrix's memory (or a copy, with `copy = true`) |
| `MatVecPntSeq` | The columns of a `D × N` matrix, read as read-only views (`ColumnPoint`); D is read from the matrix, the type is `AbsPntSeq{AnyDim,T}`; sharing as for `MatPntSeq` |

Sequences support indexing, iteration, `length`, and assignment. `PntSeq` and `GPntSeq`
also support `push!`, `pop!`, and `empty!`. `size(sequence)` returns the number
of points as an integer, following this package's finite metric space interface.
`dim(sequence)` returns the coordinate dimension.

```@example geometry
M = [0.0 3.0 0.0; 0.0 4.0 4.0]
points = MatPntSeq(M)
M[1, 2] = 6.0
(count=length(points), dimension=dim(points), second_point=points[2])
```

`CG.Points(points)` returns the underlying point vector when using
`import BasicCompGeometry as CG`.

### Shared memory and copies

Constructors do not copy their input unless asked to. Which operations reach the
input:

| Operation | `MatPntSeq(M)` (default) | `MatPntSeq(M; copy = true)` | `PntSeq(v)` from a `Vector` of points |
|:----------|:-------------------------|:----------------------------|:--------------------------------------|
| Later change to the input (`M[i, j] = x`, `v[j] = p`) | Seen by the sequence | Not seen | Seen by the sequence |
| `P[j] = p`, `translate!(P, v)`, `scale!(P, s)` | Writes into `M` | Does not reach `M` | Writes into `v` |
| `push!`, `pop!`, `empty!` | Not defined (fixed size) | Not defined (fixed size) | Resize `v` itself |
| Writing through a point, `P[j][i] = x` | Impossible: points are immutable copies | Impossible | Impossible |
| `Matrix(P)` | A new matrix | A new matrix | A new matrix |

`MatVecPntSeq` follows the `MatPntSeq` columns, with one difference: its points are
views of the matrix, not copies. A point read before the matrix changes shows the
change, and writing through a point is an error (the views are read-only).

A view or an adjoint passed to `MatPntSeq` shares the memory of its parent array, so
the same rules apply to the parent. `copy = true` stores a new dense `Matrix{T}`;
this also makes the points of a non-contiguous view (such as a subset of rows)
contiguous. For a `PntSeq` that must not share, pass a copy: `PntSeq(copy(v))`.

```@example geometry
M = [0.0 1.0 2.0; 0.0 0.0 1.0]
shared = MatPntSeq(M)
independent = MatPntSeq(M; copy = true)
translate!(shared, Point(10.0, 10.0))   # writes into M
(M = M, independent_first = independent[1])
```

Use `copy = true` when the matrix will change while the sequence is in use, or when
the sequence will be translated, scaled or assigned to and the matrix must stay as it
was. Copying costs one pass over the coordinates.

`MatPntSeq` and `MatVecPntSeq` cannot grow or shrink. Resizing a vector that shares a matrix's memory
does not behave consistently in Julia: `push!` moves it to new memory, but on Julia
1.13 a `pop!` followed by a `push!` writes into the matrix again. To change the
number of points, build a `PntSeq` from `Matrix(P)` or from the points.

### Static or `Vector` points

A `Point{D,T}` is a static vector: its length is part of its type, so the compiler
unrolls every loop over the coordinates, and creating a point costs nothing. Code is
therefore compiled anew for each dimension, and both the compile time and the
unrolled code grow with D. A `GPntSeq{D,T,Vector{T}}` holds `Vector` points instead,
whose length is known only at run time: their code is compiled once for all
dimensions and runs plain vectorized loops, but every point created is a heap
allocation.

```@example geometry
G = GPntSeq{3,Float64,Vector{Float64}}([[0.0, 0.0, 0.0], [1.0, 2.0, 2.0]])
push!(G, [2.0, 2.0, 3.0])
(count = length(G), point_type = eltype(G))
```

Measured on the distance kernels of a Fréchet distance code (FrechetDist's
`example_higher/point_types.jl`; AMD Ryzen 7 5700G, Julia 1.13, nanoseconds per call):

| D | Point–point distance, static / `Vector` | Point–segment distance, static / `Vector` | Static code compile time |
|---:|---:|---:|---:|
| 2 | 2.0 / 2.0 | 3.8 / 9.2 | 0.14 s |
| 16 | 2.2 / 3.4 | 16.1 / 23.8 | 0.17 s |
| 64 | 8.5 / 9.2 | 80.6 / 63.6 | 0.23 s |
| 128 | 41.6 / 16.3 | 193.5 / 114.0 | 0.33 s |
| 1024 | 299.2 / 125.2 | 2558.3 / 849.9 | 4.52 s |
| 2048 | 712.8 / 236.0 | 7617.5 / 1697.2 | 13.49 s |

Use static points up to D = 16 and `Vector` points from D = 128; in between the two
are close. `Vector` code compiles in 0.2 s at every D, while a whole Fréchet distance
algorithm on static points takes over ten minutes to compile at D = 6000. Creating
points is the exception: a new static point is free, so code that creates many points
stays faster with static points up to about D = 512.

For data that arrives as a matrix, the same choice is between `MatPntSeq` (static
points, for low dimensions) and `MatVecPntSeq` (column views, for high dimensions).
Reading a `MatVecPntSeq` point copies nothing, the coordinates stay contiguous in the
matrix, nothing is stored per point, and its type does not depend on the dimension:
`MatVecPntSeq(M)` is type-stable, and code compiled for it serves every D. Its column
views run the kernels above about as fast as `Vector` points, with identical results:
the distances within 8% at D = 16 and 128, creating a point up to 16% slower at
D = 128, and every kernel faster at D = 1024. Functions written for any point type accept
it (among them `dist`, `dist_sq` and `convex_comb` on its points, and `geom_length`, `prefix_lengths` and `split_edges` on the sequence); those that build `Point{D,T}` from the dimension do not.

```@example geometry
X = rand(1000, 5)                      # 5 points in dimension 1000
H = MatVecPntSeq(X)
(type = typeof(H), dimension = dim(H), d12 = dist(H, 1, 2))
```

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
