# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.7.1] - 2026-11-27

### Added
- `HilbertSpace`: a dataset known only through the squared distances of its points (a kernel function, a distance function, or coordinates read once), possibly in an infinite-dimensional space such as a kernel feature space. Each squared distance is computed at most once and cached in blocks allocated on first use; threads can share a space. `max_bytes` limits the cache, `empty!` frees it, `cache_stats` counts the calls of the distance function. A space is an `AbsFMS`, so the metric-space algorithms take it.
- `VPoint`: a virtual point `(1 - t) x_a + t x_b` on the segment of two input points of a space (24 bytes, canonical, so `==` and `hash` are geometric), with `vpoint` and `vcurve`. `dist`, `dist_sq` and `convex_comb` work on virtual points (from at most four input points per call, without allocating), and `geom_length`, `prefix_lengths`, `push_distinct!` and `split_edges` on their curves. `local_frame` places up to four virtual points isometrically in ℝ³ so that Euclidean code can run on them; `is_generator`, `line_params`, `inner_diff` (and `x - y` as a formal difference accepted by `dot`) and `vscale` are building blocks for algorithms (through `Extras`). A convex combination of points of different segments throws; `HilbertSpace(vpoints)` re-bases points into a new space.
- `AbsVirtualPoint`: `GPntSeq{D,T,P}` accepts virtual point types as `P` (with `D = AnyDim`, `T = Float64`). The bound on `P` was `AbstractVector{T}`; for coordinate points nothing changes, and their compiled code is identical.
- Guide page "Hilbert spaces and virtual points".


## [0.7.0] - 2026-11-26

### Added
- `GPntSeq{D,T,P}`: a point sequence whose point type `P` is a parameter, for example `Vector{T}` points in high dimensions, which are compiled once for all D instead of being unrolled for each D. `PntSeq{D,T}` is now an alias for `GPntSeq{D,T,Point{D,T}}`, with the same constructors, memory layout and compiled code as before; only `nameof` of the type and serialized objects see the new name. `push!` on a sequence accepts any `AbstractVector` and converts it to the point type. For points whose length is not part of their type (`Vector{T}`), the constructor and `push!` check that every point has D coordinates (for `D = AnyDim`, as many as the first point, which also gives `dim`); static points need no check, and their code is unchanged. The `GPntSeq` docstring and the geometry guide say which point type is faster in which dimension.
- `MatVecPntSeq(M; copy = false)`: a matrix-backed sequence whose points are read-only views of the columns (`ColumnPoint`), for high dimensions. The dimension is read from the matrix; the type is `AbsPntSeq{AnyDim,T}` (`AnyDim` is a new marker for a dimension known only at run time), so the constructor is type-stable and its code serves every D. Memory is shared with the matrix as for `MatPntSeq`; the size is fixed. It supports indexing, iteration, `dim`, `dist(P, i, j)`, `P[j] = p`, `translate!`, `scale!` and `Matrix(P)`, as well as `dist`, `dist_sq` and `convex_comb` on its points, and `geom_length`, `prefix_lengths` and `split_edges` on the sequence; functions that build `Point{D,T}` from the dimension do not accept it.
- Points of run-time length (`Vector` points, `ColumnPoint`s): `dist`, `dist_sq` and `convex_comb` take any two `AbstractVector`s of equal length (plain loops; `dist` does not allocate), so `geom_length` and `prefix_lengths` work on them; `push_distinct!` takes any `GPntSeq` (checking the point length); `split_edges` of a `Vector`-point sequence or a `MatVecPntSeq` returns a `GPntSeq` with `Vector` points. The static methods are unchanged and still chosen for static points.
- `MatPntSeq(M; copy = true)` copies the matrix into a new `Matrix{T}`, so the sequence and `M` are independent. The default (`copy = false`) shares memory as before. The `MatPntSeq` docstring and the geometry guide ("Shared memory and copies") now list exactly which operations reach the input matrix or vector, for `MatPntSeq` with and without copying and for `PntSeq`.

### Fixed
- DataStructures compat is `0.19`: `approx_nn`, `approx_diam` and `approx_mvbb` call `popfirst!` on a `PriorityQueue`, which DataStructures 0.18 does not define, so they failed with a `MethodError` wherever 0.18 was resolved (for example next to a package that pins `^0.18`).

## [0.6.2] - 2026-11-15

### Added
- `triangulate_polygon(...; algorithm = :auto)` (also accepted by `polygon_diagonals` and `minkowski_sum_decomposition`): the small-polygon `:ear_cut` fast path for polygons with at most 64 vertices, `:sweep` above. Opt-in; the default stays `:sweep`.

### Changed
- `triangulate_ear_cut` (and `triangulate_polygon(...; algorithm = :ear_cut)`) has a fast path for polygons with at most 64 vertices, adapted from the small-polygon kernel of [Triage](https://github.com/sarielhp/triage): reflex vertices in a 64-bit mask, ears tested lazily while walking the ring. On random simple polygons it is about 3× faster for 12–64 vertices (1.8× for 6) and allocates a third as much. On this path validation is one pass (counterclockwise orientation and consecutive equal vertices); a repeated vertex that is not consecutive is treated as any other non-simple input (an `ArgumentError` or a valid triangulation) instead of being rejected up front. Larger polygons are unchanged.
- Point-in-triangle tests of all ear-cutting algorithms reject points outside the triangle's bounding box before the orientation tests (95–99% of tests on random polygons), compared in the kernel's own coordinates (a point outside the box is outside the closed triangle). `triangulate_ear_cut` is 30–40% faster on polygons with 64–512 vertices.
- `triangulate_ear_cut` (`:ear_cut`) indexes the reflex vertices of large non-monotone polygons in the unshrunk-cell Kd-tree of `triangulate_ear_cut_u_kd_tree` instead of a plain list (kept when there are at most 16): 3× faster on 1,000-vertex random star polygons, 11× at 10,000, and 35× on 10,000-vertex combs, where it was quadratic. The small-polygon and monotone paths are unchanged.

### Deprecated
- `triangulate_ear_cut_z_order` and `triangulate_polygon(...; algorithm = :ear_cut_z_order)` now call `triangulate_ear_cut`; `triangulate_ear_cut_kd_tree` and `:ear_cut_kd_tree` now call `triangulate_ear_cut_u_kd_tree`. Neither index was faster than the unshrunk-cell Kd-tree on any tested input. Both print a deprecation warning and will be removed in 0.7. The Morton Z-order implementation is removed; `BBT.FlatTree2D` stays.

### Fixed
- A test polygon for `:ear_cut` on non-monotone input with few reflex vertices was self-intersecting; it is replaced by a simple one.

## [0.6.1] - 2026-11-14

### Added
- `add_view!(canvas, layers; autocrop=true)` adds Ipe's special `VIEWBBOX` layer to the view, so that each view of a page (each page of the compiled PDF) is cropped to its own visible objects rather than to the whole page. A `"VIEWBBOX"` listed in `layers` is registered on the canvas automatically. Example `examples/ipe_view_autocrop.jl` and a section in the vector-figures guide.

## [0.6.0] - 2026-11-14

### Changed
- Requires Julia 1.10 or newer (the LTS); Julia 1.9 is no longer supported.

### Fixed
- `is_inside` and `tangent_intersections_cp` on `CurvePolygon2D`: on Julia 1.10 each curve went through dynamic dispatch (a 5-type curve union is not union-split there), allocating 160 bytes per curve and making `tangent_intersections_cp` type-unstable.
- `read_plt` and `read_geolife_plt` closed the file only after reading every line; a malformed line left it open until garbage collection (on Windows the file could not be deleted).
- Docstrings for every exported name. `convex_hull`, `exact_diameter`, `approx_diameter` and their subspace variants had theirs written as bare `raw"""..."""` strings, which Julia does not attach (now `@doc raw`); 28 other names (the `BBT`, `WSPD`, `MVBB` and `VirtArray` modules and their accessors, `Curve2DF`, ellipse and arrangement accessors, `subdivide`, `boundary_line`) had none.
- `build_rmq` has its docstring again (`?build_rmq`); since 0.5.17 it was attached to an internal helper.

## [0.5.18] - 2026-11-14

### Fixed
- `interpolate_natural_spline(pts; closed = true)` (and `draw_spline!(...; method = :natural, closed = true)`) returned an open spline for three or more points, ignoring `closed`. It now builds a periodic C² spline with a closing segment.
- `visibility_polygon` with a vantage point within rounding distance of a polygon vertex or edge (e.g. a vertex returned by `polygon_kernel`) returned a region larger than the polygon or threw "vantage point must lie inside". Under floating kernels such a point is now snapped onto the vertex or edge (tolerance `1e-10` times the bounding-box extent), and the edge containing a boundary vantage point is ignored as an occluder. Two further floating-point defects found while fixing it: a ray cast towards a vertex could miss that vertex on an edge nearly parallel to the ray (the side test used the rounded point `q + u`; the hit is now taken at the vertex itself), and collinear-vertex pruning used an absolute threshold that dropped genuine corners next to short edges, cutting chords across the result (now a scale-free angle test, also used by `polygon_kernel`). `ExactKernel` results are unchanged.
- `dist_segment_segment` in 2D: returns 0 for crossing segments and otherwise the smallest endpoint-to-segment distance, instead of solving the normal equations. Nearly parallel crossing segments gave errors of order eps/θ² (e.g. 3e-8 instead of 0 for length-200 segments), depended on the endpoint order, and below θ ≈ 1e-8 returned the endpoint distance.
- Polygon triangulation (`triangulate_polygon`, `triangulate_monotone`, `triangulate_ear_cut*`, `triangulate_sweep`, `triangulate_seidel`) never returns clockwise or overlapping triangles: the result is a valid triangulation (zero-area triangles allowed for weakly simple input) or an `ArgumentError` saying the polygon is not simple or is numerically degenerate for the kernel. `monotone_chain_indices` (used by the `:ear_cut` monotone fast path) and the monotone-piece triangulation shared by `:sweep` and `:seidel` split the upper and lower chains with a local turn test at the leftmost vertex. That test is exactly zero on a zero-width spike and can round to the wrong sign under `FloatingKernel`, and either way it swapped the chains. For example, `:ear_cut` on the strictly simple quadrilateral `(0.74, 0.21), (-0.35, 0.56), (-0.03, -0.56), (-0.07, -0.42)` returned a clockwise triangle whose areas summed to 1.25 times the polygon's. The chains now come from the polygon's signed area (pieces are CCW by construction). Fan triangles of the monotone stack triangulation and the final triangle of ear clipping, previously emitted unchecked, now throw when they come out clockwise (crossing chains, or inconsistent floating-point predicates). With this, all ear-cutting algorithms triangulate weakly simple polygons, and `:sweep`/`:seidel` either triangulate them or reject them.

## [0.5.17] - 2026-11-14

### Fixed
- `at`, `sample_uniformly`, and `dist(::AbsPntSeq, i, j)`: eliminated type instabilities (`Union{Point{D,T}, Point{D,Float64}}` / `Union{Float32, Float64}`) for `Float32` and integer point sequences by consistently using `float(T)` (`zero(float(T))` when `i == j`), and fixed `InexactError` in `sample_uniformly` on integer point sequences.
- `visibility_polygon` and `polygon_kernel`: integer coordinate inputs under `FloatingKernel()` now promote to `Float64` instead of throwing `InexactError`, and `Float32` / `Rational{Int64}` coordinate types are preserved with homogeneous, type-stable return types (including `_ray_edge_intersect`).
- `WSPD` and `approx_diameter`: `WSPD.init` now derives its distance scalar type from the `BBT.Tree` rather than `sep` and accepts any `sep::Real` (fixing `MethodError` on `Float32` point sequences and integer `sep` values), while `WSPD.separation` preserves the tree's distance scalar type `S` and `WSPD.PD` parameterizes the concrete point-sequence type `P`.
- `dist_subspace`, `nn_point`, and `bisection_point`: eliminated value-dependent union return types for integer, `Float32`, and `Rational` inputs so degenerate/nonpositive branches return the exact same promoted coordinate/scalar types as the general branches.
- `point_on(::AbsPntSeq; closed=true)`: standardized all return branches to `Point{D,float(T)}` (eliminating `Union{Point{D,T}, Point{D,Float64}}` type instability on integer and `Float32` sequences) and eliminated `vcat` / `zeros` allocations via a zero-allocation two-pass walk with zero-length edge guarding.

- `dist_segment_segment`: two zero-length segments no longer throw `SingularException`.
- `intersect_line_curve` on a zero-span `CircleArc` returns the tangent hit instead of `NaN`.

### Compatibility notes
- `WSPD.PD` and `PntSeqHierarchy.SimplificationHierarchy` gained a type parameter; code that spells out the old parameter list must be updated.
- Under `FloatingKernel()`, `visibility_polygon` and `polygon_kernel` always return `Float64` points for integer input; `at`, `nn_point`, and `point_on` always return float points.
- `vertical_ray_segment` rejects infinite rational coordinates (`1//0`) with `ArgumentError`.
- `read_plt` uses the first `dim` numbers of each line and ignores trailing fields that used to raise an error.
- `BBT.Tree_draw` fills with alpha exactly 0.1 (was 26/255) and draws nodes level by level; the pictures are visually unchanged.

### Changed
- `rand_point` and `rand_gaussian`: added `Val{D}` overloads for compile-time dimension inference and replaced `randn(T, D)` in `rand_gaussian` with stack-allocated `ntuple` construction (zero heap allocations).
- `dist_segment_segment`: replaced heap `Vector`/`Matrix` construction and `m \ c` linear solve with zero-allocation scalar Cramer's rule; parallel segments are detected by `Δ > 0` with no tolerance, so nearly parallel crossing segments no longer fall back to the endpoint distance.
- `apply_transform`: widened `Matrix{T}` to `AbstractMatrix{T}` and replaced augmented-vector heap allocations with inline `muladd` evaluation (zero heap allocations).
- `union_bbox(objects::AbstractVector)`: replaced vararg splatting and `map(bbox, objects)` with single-pass streaming reduction via `_union_bbox_iter` (10,000-box allocations reduced from ~13.2 MB to 96 bytes).
- `RMQ` and `ChainIndex`: eliminated `Core.Box` captures and missing comparator specialization in `_sparse_levels` and `build_rmq` (with `Base.@constprop :aggressive` and `Val`-dispatched backends), making `build_rmq` and `build_chain_index` concretely inferred and removing per-query heap allocations in `ChainWindow` queries.
- `BBT.get_min_orig_index` / `BBT.get_max_orig_index`: replaced dynamic `getfield`/`setfield!` reflection and generator allocations with direct typed field access and a tight loop, eliminating heap allocations during traversal.
- `triangulate_monotone`: eliminated inner-loop stack vector reallocations on alternating upper/lower chains by reusing the allocated buffer and pre-sizing capacity; eliminated redundant `orient_sign` predicate evaluations by emitting counterclockwise triangles directly from chain classifications (`is_upper`). Measured up to ~3x speedup under `FloatingKernel()` (4k allocations reduced to 38) and ~97x speedup under `ExactKernel()` (280k allocations reduced to 173).
- `_check_distinct_vertices`: eliminated intermediate $O(n)$ tuple array allocation by sorting vertices directly by coordinate projection with `MergeSort`, reducing memory allocation by ~57% and eliminating GC stalls on large polygons.
- `_sweep_partition_diagonals`: eliminated redundant $O(n)$ heap vector allocation for `prev_v` via inline cyclic index computation, and pre-sized `diags` with `sizehint!(diags, n - 3)`.
- `visibility_polygon` and `polygon_kernel`: eliminated `Core.Box` closure captures via coordinate-type function barriers and helper extraction, avoided intermediate edge-tuple allocation, and reused Sutherland–Hodgman scratch buffers across clipping iterations (40-vertex star allocations reduced from 656 KB to 4.9 KB and 385 KB to 5.0 KB; 200-vertex star allocations reduced from 15.0 MB to 22.8 KB and 7.0 MB to 20.6 KB).
- `is_convex`, `minkowski_sum`, and `minkowski_sum_decomposition`: eliminated `Core.Box` captures, concrete-typed `Point{2}` / `Segment{2}` wrapper vectors (`_pts_vec`), replaced intermediate rotated/edge vector allocations with direct cyclic indexing and zero-allocation `_TrianglePieces` views, and reused merge/deduplication scratch buffers across triangle pairs (`is_convex` 3.5 KB -> 128 B; `minkowski_sum` 10.8 KB -> 1.1 KB and 8.8 KB -> 960 B; 40-vertex star `minkowski_sum_decomposition` 5.1 MB -> 274 KB).
- `LongestConvexSubset`: parameterized `LCSContext` and `RotateOrder` on the concrete sequence type `P<:AbsPntSeq{2,T}`, preallocated `PointTurn` arrays in `pt_init`, and hoisted `lcs_init` out of the `target` loop in `compute_largest_convex_subset` (60-point allocations reduced from ~192.5 MB to ~175 KB).
- `MVBB.approx_diam` and `MVBB.approx_diam_subspace`: parameterized priority-queue keys on the concrete `BBT.Node{D,T,S,V}` type via `_init_node_pair_pq` and specialized `_push_child_pairs!` on the upper-bound function.
- `ConvexHull3D.convex_hull_3d`: replaced per-horizon-edge conflict vector concatenation and `unique!` with a reusable generation-counter stamp array (`seen_stamp`), reducing 400-point spherical shell allocations from ~4.64 MB to ~2.69 MB.
- `vertical_decomposition` and `trapezoid_chains`: concrete-typed `_ActiveTrap{R_type}` double-buffering, preallocated slab buffers with $O(k)$ two-pointer advancement, binary-searched zero-allocation `_vd_walls`, and indexed `x_left` neighbor linking in `_tc_link_neighbors!` (160-line decomposition allocations reduced from ~2.8 GiB to ~24 MiB and chain linking from ~1.06 s to ~9 ms).
- `vertices(...; unique=true)` (`_unique_points`): replaced $O(m^2)$ scan with sorted $x$-window deduplication preserving original order.
- `sample_uniformly`: replaced $O(mn)$ repeated `at` / `prefix_lengths` calls with a single-pass $O(m + n)$ scan (1,000-point resampling allocations reduced from ~8.1 MB to ~24 KB).
- `PntSeqHierarchy.SimplificationHierarchy`: parameterized `plt::Palette{T,R}` on `R<:SplitRule` to eliminate dynamic dispatch in `hierarchy_approx` / `hierarchy_build!`, and replaced `0.0` base-case returns in `_fill!` / `_p_extract!` with `zero(eltype(...))`.
- `vertical_ray_segment` and `first_positive_vertical_contact`: added exact `x`-span / ray-side rejection before `BigInt` / `Rational{BigInt}` widening, eliminating heap allocations on non-spanning segments.
- `CubicBezier` and `CubicSpline` (`bbox`, `flatten`, `interpolate_natural_spline`): replaced heap vector allocations and coordinate splatting with stack-allocated `ntuple(..., Val(D))` extrema and control-point construction, restoring concrete `BBox{D,T}` return-type inference and eliminating per-segment `BBox`/`PntSeq`/slice allocations.
- `CurvePolygon2D` (`is_inside`, `tangent_intersections_cp`, `intersect_line_curve`): replaced per-curve intersection vector allocations and `prefix_lengths` array construction with stack-allocated `_LineCurveHits` and streaming arc-length accumulation (zero heap allocations on `Segment2F`/`CircleArc2F` polygons).
- `bbox(::EllipticArc)`: replaced temporary point-vector and generator allocations with direct scalar extrema tracking.
- `IpeDraw`: parameterized `Style{NT<:NamedTuple}`, routed `draw!`/`mark!`/`label!` through concrete style function barriers with compile-time `NamedTuple` option filtering, stabilized `_apply_tf` to always return `Point{2,Float64}`, replaced per-label alignment `Dict` allocation and multi-pass `replace` escaping with constant lookup and single-pass `IOBuffer` streaming, and streamed `draw_polygon!`, `draw_spline!`, `draw_bspline!`, `draw_polygon_with_holes!`, and `draw_label!` directly into pre-sized `IOBuffer`s (`draw_bspline!` now throws `ArgumentError` for fewer than 3 control points).
- `ReadWrite.read_plt` and `ReadWrite.read_geolife_plt`: added a `Val{D}` function barrier with `Base.@constprop :aggressive` in `read_plt` for concrete `PntSeq{D,Float64}` inference (throwing `ArgumentError` for `dim < 1`), streamed files via `eachline` instead of `readlines`, and replaced regex/comma `split` allocations with zero-allocation `SubString` scanning.
- `ReadWrite.write_fbin` and `ReadWrite.write_ibin`: replaced lazy `Adjoint` scalar IO (`write(f, vecs')`) with bounded (≤ 1 MiB) chunked transposition into a contiguous buffer dispatched to bulk `unsafe_write`.
- `BBT.Tree_draw` (`BBTCairoExt`): pre-collected internal nodes by depth in a single $O(N)$ pass instead of $\Theta(d \cdot N)$ unpruned per-page DFS traversals, and replaced per-node `Colors.parse(Colors.Colorant, ...)` lookups with precomputed concrete `RGB`/`RGBA` constants (`tree.root === nothing` now throws `ArgumentError`).

## [0.5.16] - 2026-11-13

### Fixed
- Ear clipping (`:ear_cut` and the Z-order / Kd-tree variants) now throws `ArgumentError` when no ear remains instead of silently returning a partial triangulation (e.g. on numerically degenerate Float64 input; use `ExactKernel()` or `:sweep`).
- Polygon validation (`_validate_ccw_polygon`) rejects duplicate vertices for `:ear_cut`, `:seidel`, and the indexed ear cutters, which previously crashed, returned zero-area triangles, or truncated output.
- `triangulate_seidel`: trapezoid walks throw a clear `ArgumentError` on inconsistent floating-point predicates instead of `BoundsError` / internal errors or a possible infinite loop; triangles are emitted CCW from the chain structure instead of being silently re-oriented; its docstring is attached again.
- `BBT.Tree_init(...; leaf_size > 1)`: `approx_nn`, `silly_nn`, and min/max original-index queries scan every point of multi-point leaves (they returned wrong nearest neighbors). `leaf_size < 1` now throws.
- Flat Kd-trees (`Tree_init_flat`, `Tree_init_flat_u`): construction no longer hangs on ±Inf/NaN coordinates (non-finite input is rejected), original ids are validated (out-of-range ids wrote out of bounds), and huge integer coordinate spans no longer overflow.
- `BBT.Node` / `BBT.Tree` have compact `show` methods; the new `parent` field made printing any node dump the entire tree.

### Changed
- Ear clipping (`:ear_cut`, `:ear_cut_z_order`, `:ear_cut_kd_tree`, `:ear_cut_u_kd_tree`) uses a lazy re-testing ear queue: neighbors of a clipped ear are re-queued at the back and tested only when dequeued. Measured new/old time: u_kd and kd about 0.55-0.8x on star, coast and spiral inputs (neutral on altstar and combs), Z-order 0.13x on spirals. **The ear order changed, so `:ear_cut` and the indexed ear cutters return different diagonals than before** (they still agree with each other).
- `:ear_cut_kd_tree` / `:ear_cut_u_kd_tree`: Float64 queries pre-scan each Kd-tree leaf with a branch-free loop and run the exact scan only on a possible hit (0.70x on random stars, 0.78x on coastlines; up to ~10% slower on inputs with sparse leaves).
- `triangulate_seidel` is about 2.5x faster and allocates about 2.6x less memory (isbits `Int32` trapezoid/DAG records, flat piece storage, reused scratch buffers); output is unchanged for a given `rng`. Inputs above `typemax(Int32) ÷ 16` vertices now throw `ArgumentError` (use `:sweep`).
- `triangulate_polygon(...; algorithm = :ear_clipping)` is accepted again as an alias of `:ear_cut`.
- Corrected complexity claims for `:ear_cut_z_order` (degrades toward quadratic on long thin ears and far outliers) and the Kd-tree ear cutters (subquadratic in practice, `O(n r)` worst case).

## [0.5.15] - 2026-11-13

### Added
- Z-order ear clipping (`triangulate_ear_cut_z_order`, `:ear_cut_z_order`) with a Morton index over active reflex vertices.
- Kd-tree ear clipping (`triangulate_ear_cut_kd_tree`, `:ear_cut_kd_tree`) over the new flat array-backed `BBT.FlatTree2D` (`BBT.Tree_init_flat`), and the unshrunk-cell variant (`triangulate_ear_cut_u_kd_tree`, `:ear_cut_u_kd_tree`, `BBT.FlatUTree2D`, `BBT.Tree_init_flat_u`).
- `Kd_tree` alias for `BBT`; `parent` navigation for `BBT.Node` and flat-tree node ids; `leaf_size` and `compute_diam` keywords for `BBT.Tree_init`.
- `algorithm` keyword for `minkowski_sum_decomposition`; exported `triangulate_ear_cut`.
- Example `examples/polygon_triangulation.jl` comparing all triangulation algorithms.

### Changed
- **Breaking:** the default `triangulate_polygon` / `polygon_diagonals` / `minkowski_sum_decomposition` algorithm is now `:sweep` (previously ear clipping). Results use different diagonals, and inputs that are not strictly simple (e.g. a vertex touching another edge) are rejected with `ArgumentError`.
- **Breaking:** `:ear_clipping` was renamed `:ear_cut` (the old name is accepted again since the next release).

## [0.5.14] - 2026-11-13

### Added
- Deterministic worst-case `O(n log n)` time and `O(n)` space plane-sweep simple polygon triangulation (`triangulate_sweep` and `triangulate_polygon(...; algorithm = :sweep)` in `src/SweepTriangulation.jl`), sorting events with `MergeSort` under `compare_xy` symbolic shear `(x + εy, y)`, maintaining active lower-boundary edges in a preallocated AVL tree (`_SweepAVLTree`), and triangulating monotone pieces in `O(n)` total time via `_triangulate_from_partition`.
- Exported `triangulate_sweep` (`src/BasicCompGeometry.jl`, `test/runtests.jl`), documented it in `docs/src/algorithms.md` and `docs/src/reference/algorithms.md`, and added invariant, degeneracy, and `O(n log n)` scaling tests in `test/test_monotone_triangulation.jl`.

### Fixed
- Ear clipping (`src/MonotoneTriangulation.jl`, now `_ear_clip_run!`) tracks reflex vertices in a doubly-linked list (`_EarClipState`) and candidates in an ear queue so it runs in `O(n r)` time (`O(n)` for `r == 0`) rather than `Θ(n^2)`.
- `convex_hull_3d` (`src/ConvexHull3D.jl`) now maintains reverse conflict adjacency `point_conflicts` incrementally instead of scanning all faces per point, matching its expected `O(n log n)` runtime.
- `discrete_hull_columns` (`src/DiscreteHull.jl`) now uses linear-time `_lattice_hull_sorted` on unit-step columns to achieve its documented `O(r)` runtime without `O(r log r)` sorting.
- Added complexity and scaling regression tests in `test/test_convex_hull_3d.jl` and `test/test_monotone_triangulation.jl`.

## [0.5.13] - 2026-11-13

### Changed
- Replaced the previous `O(n^3)` placeholder in `triangulate_seidel` with Raimund Seidel's (1991) expected `O(n log* n)` randomized incremental polygon triangulation algorithm (`src/SeidelTriangulation.jl`), using symbolic shear `(x + εy, y)` via `AbstractKernel` predicates (`compare_xy` / `orient_sign`), `O(log* n)` phases interleaved with `O(n)` polygon boundary tracing, linear-time monotone piece extraction, and linear-time monotone piece triangulation.
- Updated `triangulate_seidel` tests (`test/test_monotone_triangulation.jl`) to check trapezoidal-map and DAG invariants, degenerate inputs (vertical/horizontal edges and collinear vertices), `log* n` phase bounds, and linear predicate-count scaling.

## [0.5.12] - 2026-11-01

### Added
- Discrete (lattice) convex hulls, `ConvexHull(C ∩ Z²)`, by the output-sensitive algorithm of Har-Peled (CGTA 1998): `discrete_hull` makes O(|hull| log diam(C)) chord queries, e.g. O(r^(2/3) log r) for a disk of radius r.
  - Bodies are accessed only through `chord_range` (exact `(⌈t_lo⌉, ⌊t_hi⌋)` of a line chord). Exact bodies: `LatticeDisk` (rational center and radius) and `ConvexPolygonBody` (rational vertices); `FloatChordBody` wraps any floating-point chord oracle.
  - `find_edge_direction` (continued-fraction edge search in a lattice cone), `discrete_hull_start` (Lemma 3.10 start vertex), `discrete_hull_columns` (naive O(r) baseline), and `lattice_convex_hull` (exact integer monotone chain).
  - Example `examples/discrete_hull_disk.jl`: a disk of radius 10^7 in under 0.1 s, verified against the naive hull.
- Discrete hulls documented in the manual and README.

### Fixed
- Documentation build: the RMQ, chain-index and incidence-predicate exports now appear in the algorithm reference, and an unresolvable `ArgumentError` cross-reference in the `polygon_kernel` docstring was removed.

## [0.5.11] - 2026-10-30

### Changed
- The `--help` output of `tools/release`, `tools/bump`, `tools/publish` and `tools/push` describes what each tool does.

## [0.5.10] - 2026-10-26

### Added
- `build_rmq` / `range_argmin`: static linear-space, linear-preprocessing, O(1)-query range-minimum index (rightmost minimum on ties; `:sparse` reference backend).
- `build_chain_index` / `chain_window` (+ `nedges`, qualified `vertex`, `edge`, `bounds`, `leftmost`, `point_at`): shared cyclic-chain geometry with distinct closing occurrence, O(1) windows and exact interpolation.
- `vertical_ray_segment` (`NoHit`, `PointHit`, `OverlapHit`) and `first_positive_vertical_contact`: exact closed vertical-ray contacts.
- `open_sector_contains` and `incoming_sector_contains`: exact angular-sector predicates, including slit tips and the symbolic incoming direction.
- All constructions are exact (widened to `BigInt`/`Rational{BigInt}`); floating-point coordinates are rejected.

## [0.5.9] - 2026-10-25

### Added
- Minkowski sum of 2D polygons and shapes:
  - `minkowski_sum` computes the exact Minkowski sum $P \oplus Q = \{ p + q \mid p \in P, q \in Q \}$ in optimal $O(n + m)$ time for convex polygons by merging directed edges cyclically ordered by polar angle.
  - Overloaded for `Point{2}`, `Segment{2}`, `AbsPntSeq{2}`, and `AbstractVector{<:Point{2}}` (producing translations, segments, stadium swept polygons, and parallelograms).
  - Optional `hull=true` computes the convex hull $\mathrm{CH}(P \oplus Q) = \mathrm{CH}(P) \oplus \mathrm{CH}(Q)$ in $O(n \log n + m \log m)$ time for arbitrary input polygons.
  - `minkowski_sum_decomposition` decomposes non-convex polygons into convex pieces (triangles) and computes convex Minkowski sums for each pair.
  - `is_convex` predicate verifies polygon convexity (supporting collinear edge vertices, orientation-agnostic, and ruling out self-intersecting star shapes).
  - Supports both `FloatingKernel()` and `ExactKernel()` (preserving exact rational coordinates).
  - Added `Base.eltype` implementation for `AbsPntSeq` preserving element type across collections.

## [0.5.8] - 2026-10-25

### Added
- Polygon kernel and star-shapedness:
  - `polygon_kernel` computes the convex kernel $\mathrm{Ker}(P) = \{ q \in P \mid \forall p \in P, \overline{qp} \subset P \}$ of points that see the entire simple polygon $P$.
  - `is_star_shaped` tests whether a simple polygon has a non-empty kernel.
  - Supports both `FloatingKernel()` and `ExactKernel()` (preserving exact rational coordinates).
  - Fixed comprehension type inference in `Arrangements.jl` for compatibility across Julia 1.9, 1.10, and 1.11+.

## [0.5.7] - 2026-10-25

### Added
- Visibility polygon from a vantage point:
  - `visibility_polygon` computes the star-shaped visibility polygon $\mathrm{Vis}(q, P)$ consisting of all points inside counterclockwise simple polygon $P$ visible from vantage point $q \in P$.
  - Handles vantage points strictly in the interior, along boundary edges, and at convex or reflex vertices.
  - Supports both `FloatingKernel()` (Float64) and `ExactKernel()` (preserving exact rational coordinates without floating-point approximations).
  - Extended `is_inside(q, poly)` predicate for testing point containment in simple polygons.

## [0.5.6] - 2026-10-24

### Added
- Seidel's randomized incremental polygon triangulation algorithm:
  - `triangulate_seidel` implements Raimund Seidel's randomized incremental polygon triangulation algorithm in expected $O(n \log^* n)$ time, decomposing the polygon into monotone sub-polygons via reflex cusp resolution and triangulating each piece in linear time.
  - `triangulate_polygon` now supports an `algorithm` keyword argument accepting `:ear_clipping` (default) or `:seidel`, as well as a reproducible `rng` keyword.

## [0.5.5] - 2026-10-24

### Added
- General simple polygon triangulation:
  - `triangulate_polygon` triangulates an arbitrary simple counterclockwise polygon with $n \ge 3$ vertices, returning $n - 3$ internal diagonals and $n - 2$ counterclockwise triangles.
  - Automatically detects $x$-monotone or $y$-monotone polygons via `is_x_monotone` and `is_y_monotone` and triangulates them in linear time $O(n)$; general simple polygons are triangulated via ear clipping with reflex acceleration in $O(n r)$ time where $r$ is the count of reflex vertices.
  - `PolygonTriangulation` result container supporting property access (`res.diagonals`, `res.triangles`) and destructuring `(diagonals, triangles) = res`.
  - `is_x_monotone` and `is_y_monotone` predicates for checking strict polygon monotonicity.
  - `polygon_diagonals` returning internal diagonals as sorted index pairs `(i, j)` with $i < j$.
  - `MonotoneTriangulation` is preserved as a backward-compatible alias for `PolygonTriangulation`.

## [0.5.4] - 2026-10-24

### Added
- Exact geometric primitives and predicates for rational and exact kernels:
  - `orient` returning the exact orientation determinant and `orient_sign` returning the normalized sign preserving zero.
  - `on_segment` testing closed segment membership including endpoints and vertical segments.
  - `proper_crossing` testing strict interior segment-segment crossing.
  - `midpoint` and `affine_interpolate` for exact constructed points without floating-point fallback.
  - `segment_at_x` for exact nonvertical segment evaluation at a given x-coordinate (enforcing nonvertical precondition).
  - `compare_angle` and `angle_less` for counterclockwise angular order of direction vectors starting from the positive x-axis in $[0, 2\pi)$ with same-direction equivalence and nonzero vector precondition enforcement.
  - `segment_key` returning canonical undirected segments with endpoints ordered lexicographically; `Segment` now supports iteration, tuple unpacking `(p, q) = seg`, and integer indexing.
  - `id_pair` returning canonical sorted index pairs.
  - `twice_area` computing signed double area using the shoelace formula and `polygon_segments` extracting cyclic directed boundary segments.
- Linear-time monotone polygon triangulation:
  - `triangulate_monotone` decomposes a simple counterclockwise x-monotone polygon into $k - 3$ distinct internal diagonals and $k - 2$ counterclockwise triangles in deterministic $O(k)$ time and $O(k)$ space without comparison sorting or coordinate dictionaries.
  - `MonotoneTriangulation` container type supporting property access (`res.diagonals`, `res.triangles`) and destructuring `(diagonals, triangles) = res`.
  - `monotone_chains`, `monotone_chain_indices`, `merge_chains`, and `stack_diagonals` helpers.

## [0.5.2] - 2026-10-20

### Changed
- No behaviour change: `IpeDraw.jl` is split into topic files under `src/IpeDraw/`, and
  over-complex functions in the Cairo extensions, arrangements, `BBT`, halfplanes,
  longest convex subset, `ConvexHull3D`, `MVBB` and `Bezier` were split into helpers
  (verified bit-identical on seeded random inputs against 0.5.1).
- `BBox{D,T}` now has an explicit keyword constructor instead of one generated by
  `@kwdef`; `BBox{D,T}(; f_init, mini, maxi)` behaves as before.
- Release tooling: `tools/release` (tests, bump, push), `tools/publish` (CI wait and
  registration) and `tools/bump` replace `scripts/bump_version`.

## [0.5.1] - 2026-10-19

### Added
- `PntSeqHierarchy` submodule: palettes of recursive curve subdivision (midpoint or
  farthest-vertex split rules), palette extraction by error threshold or level, and
  `SimplificationHierarchy`, a width-indexed set of simplifications driven by a
  pluggable `HierarchyMeasure`. The default measure uses the Hausdorff distance;
  other measures (such as a Frechet distance) supply their own error, simplification
  and certification functions. Not exported by default.
- `push_distinct!` and `split_edges` for `PntSeq` (both in `Extras`, which now exports
  91 names), and `Base.reverse(::PntSeq)`.
- `ReadWrite.read_plt` (the inverse of `write_plt`: comma or whitespace separated
  coordinates, any dimension) and `ReadWrite.read_geolife_plt` for GeoLife trajectories.
- `push_distinct!` and `split_edges` for `PntSeq` (both in `Extras`, now 91 names),
  `Base.reverse(::PntSeq)`, and `ReadWrite.read_plt` / `ReadWrite.read_geolife_plt`
  for text files of planar curves.

## [0.5.0] - 2026-10-16

### Added
- Coverage for ordinary constructors, accessors, curve intersections, binary-file
  subsets and errors, bichromatic WSPD, hull diagnostics, and drawing helpers.
- Reproducible examples environment and bounded smoke tests for all 22 example
  programs, including all three prophet modes, with output validation.
- Opt-in `BasicCompGeometry.Extras` submodule exporting 89 additional supported
  public API names, with their original bindings.
- CI publication of unit, example, and external-tool integration coverage to
  Codecov, alongside downloadable coverage reports.

### Changed
- Removed the unused `MetricSpace.greedy_permutation_vanity` routine and its
  private helper from the active API; preserved their source under `archive/`.
- Fixed `angle(::EllipticArc)` to extend `Base.angle` and made Hausdorff
  simplification accept one-point curves.
- Fixed stale polygon-drawing calls, bounding-box-tree color rendering, late Cairo
  loading in the halfplane example, and the sphere example's table backend.
- **Breaking:** Reduced the default export list to 130 public API names. The
  89 names listed below are available through qualification, explicit imports,
  or `Extras`. Their definitions, signatures, and behavior are unchanged.
- Core geometry types, useful coordinate aliases, common operations, and main
  algorithms remain exported. Specialized helpers, auxiliary namespaces,
  less central compatibility aliases, and Cairo/LaTeX utilities move to `Extras`.
- Updated internal modules, extensions, examples, tests, and documentation to
  import additional helpers explicitly. Feature submodules keep their own exports.
- Shortened the README and reorganized documentation into focused guides and
  API reference pages, with executable examples and strict exported-docstring
  and reference checks.
- Consolidated the development-registration helper under one lowercase filename
  for case-insensitive filesystems.
- Removed the compiled C++ test executable and a generated example PDF from
  tracked release files; their source programs remain available.

### Migration

Existing calls such as `width(bb)` after `using BasicCompGeometry` must either
qualify the function, explicitly import it, or opt into `Extras`:

```julia
import BasicCompGeometry as CG
bb = CG.BBox(CG.point(0.0, 0.0), CG.point(2.0, 3.0))
CG.width(bb) # 2.0

# Alternatively, keep unqualified calls with an explicit import:
using BasicCompGeometry: width
width(bb) # 2.0
```

To import both the default API and all additional public helpers:

```julia
import BasicCompGeometry as CG
using .CG, .CG.Extras

bb = BBox(point(0.0, 0.0), point(2.0, 3.0))
width(bb) # 2.0
```

`Extras` exports only the explicitly listed public API names below. It shares
their original bindings and does not export internal implementation utilities.

| Group | Names available through `Extras` or explicit imports |
|---|---|
| Generic geometry helpers | `area`, `boundary`, `clip`, `complement`, `depth`, `direction`, `distance`, `faces`, `height`, `lines`, `locate`, `vertices`, `width` |
| Geometry accessors | `a_coeff`, `alpha1`, `alpha2`, `axis_direction`, `b_coeff`, `c_coeff`, `center`, `end_point`, `p_coeff`, `r_major`, `r_minor`, `rotation_angle`, `start_point`, `vertex` |
| Curve queries | `at`, `derivative` |
| Compatibility and container aliases | `AbsPolygon`, `MatPolygon`, `Polygon_random_sphere`, `VecPntSeq2F`, `VecPoint2I`, `VecPolygon2F`, `VertTrapezoid`, `defaultview`, `rand_polygon`, `random_halfplane` |
| Specialized geometry helpers | `after_tangent_to_polygon`, `before_tangent_to_polygon`, `bisection_point`, `compare_right_arc_to_point`, `compare_right_arcs`, `convex_coef`, `dist_point_segment`, `dist_segment_segment`, `dist_subspace`, `distance_infty`, `hausdorff_dist_subseg`, `match_price`, `max_dist_subspace`, `pnt_seq_random_sphere`, `prefix_lengths` |
| Bounding-box helpers | `bound!`, `d_max`, `d_min`, `expand!`, `expand_add!`, `init!`, `middle` |
| Cairo and LaTeX utilities | `Canvas`, `add_latex_macros!`, `add_latex_packages!`, `append_latex_preamble!`, `cairo_draw_halfplane`, `cairo_draw_halfplanes`, `cairo_draw_latex`, `cairo_draw_latex_page`, `cairo_draw_points`, `cairo_draw_polygon`, `cairo_draw_setup`, `cairo_set_line_width`, `description`, `get_file_path`, `get_latex_preamble`, `latex_to_pdf`, `open_canvas`, `pdf_merge`, `read_latex_snippet`, `reset_latex_preamble!`, `set_latex_preamble!` |
| File I/O | `read_halfplanes`, `write_halfplanes`, `write_plt` |
| Auxiliary namespaces | `ReadWrite`, `VirtArray` |
| Sequence helpers | `Points`, `cardin` |

Canonical spellings such as `default_view` and `rand_halfplane` remain exported.
Use `length(ps)` for a point sequence's cardinality, or explicitly import
`cardin`. `at` and `point_on` are distinct operations; moving `at` does not imply
that they use identical parameter conventions. Visualization still requires
loading Cairo and, for LaTeX utilities, LaTeXStrings.

`BasicCompGeometry.depth` still supports both halfplane queries and bounding-box
trees. `BasicCompGeometry.BBT.depth` remains exported within the `BBT` submodule.
When combining geometry packages, core type names such as `Point`, `Line`,
`Polygon`, `Circle`, and `Sphere` may still overlap; use module qualification or
explicit imports to select the intended type.

## [0.4.1] - 2026-10-13

### Changed
- Prepared the first General registry release with macOS and Windows CI alongside Linux.
- Added downloadable coverage reports to CI.
- Resolve documentation dependencies against the checked-out package instead of a committed manifest.
- Release tooling waits for CI and requests registration; TagBot handles tags and releases.

## [0.4.0] - 2026-10-13

### Added
- Pluggable geometric predicate kernels (`AbstractKernel`, `FloatingKernel`, `ExactKernel`).
- Exact predicates for arbitrary-precision integer and rational coordinates, and exact binary float values.
- Exact algebraic radical isolation for circle boundary comparisons (`compare_right_arcs`, `compare_right_arc_to_point`) without floating-point square root approximations.
- Predicates `compare_xy`, `compare_squared_distance`, `compare_distance`.
- Kernel-aware 2D convex hull via Monotone Chain.

### Changed
- `turn_sign` and related orientation predicates (`is_left_turn`, `is_right_turn`, etc.) now accept an optional `kernel` parameter, defaulting to `FloatingKernel{Float64}()`.
- `turn_sign` now returns `Float64` by default for integer and `Float32` inputs; use `kernel=ExactKernel()` when exact arithmetic is required.
- Non-finite coordinates in predicates now raise `DomainError`.

## [0.3.0] - 2026-10-05

### Added
- Geometry-aware Ipe figure API with automatic fitting, reusable styles, page-space labels, and direct PDF targets that retain editable Ipe sources.
- Circle-circle intersections and bounding-box unions.
- Release CI, package-quality checks, and adversarial WSPD tests.

### Changed
- Prepared version 0.3.0 with a permanent package UUID.
- Made Cairo, Colors, and LaTeXStrings optional extension dependencies.
- Reduced the core dependency set and completed compatibility bounds.
- Replaced componentwise `min` and `max` methods for points with `point_min` and
  `point_max`, avoiding method piracy against `StaticArrays.SVector`.
- Removed the root manifest so each supported Julia release resolves compatible
  dependency versions.
- Clarified supported functionality, maintenance expectations, and numerical limitations.

## [0.2.3] - 2026-07-30

### Added
- `Curve2D{T}` type alias (`Union{Segment{2, T}, CircleArc{T}}`) and `CurvePolygon2D{T}` struct representing polygons formed by curves between consecutive vertices.
- Total length computation `geom_length(cp::CurvePolygon2D)`.
- Parametric position evaluation `point_on(cp::CurvePolygon2D, t)` and `point_on(pnt_seq::AbsPntSeq, t; closed=false)`.
- `direction(alpha)` returning 2D unit vector at angle $\alpha$.
- `tangent_line(poly, alpha)` returning tangent line to convex polygon with outer normal `direction(alpha)`.
- `before_tangent_to_polygon` and `after_tangent_to_polygon` computing parametric intersection locations $u, v$ on `CP` for tangent lines to inner polygon $Q$.
- Extended multi-page visualizations in `examples/polar_bisectory.jl`:
  - Page 6: CurvePolygon2D `CP` outer boundary representation.
  - Page 7: Tangent line to inner polygon $Q$ at angle $\alpha$, tangency point $q$, and intersection points $u, v$ on `CP`.
  - Pages 8 & 9: Periodic plots of $f(\alpha)$ and $g(\alpha)$ with vertical interval segments between curves.
  - Page 10: Staircase algorithm $S(x_1)$ with ray-shooting and self-intersection termination.
- Animated MP4 movie (`output/staircase.mp4`) and 200 frame PNG generation (`output/frames/%06d.png`).

## [0.2.2] - 2026-07-23

### Added
- `Sphere{D, T}`, `Circle{T}`, and `CircleArc{T}` primitives with type aliases (`Circle2F`, `Sphere2F`, `Sphere3F`, `CircleArc2F`).
- Circular inversion (`invert`) for `Point`, `Line`, `Sphere`, and `Segment` in 2D.
- 2D `Line` $\leftrightarrow$ `Plane` conversion constructors (`Line(::Plane)` and `Plane(::Line)`).
- `polar_line(p)` and `polar(Line, p)` methods for explicit 2D Point $\iff$ Line duality; 2D `polar(p)` now defaults to `Line{2, T}`.
- In-place bounding box merging method `bound!(bb1::BBox, bb2::BBox)`.
- Multi-page visualization outputs in `examples/polar_bisectory.jl`: `bisectory_random.pdf`, `bisectory_square.pdf`, and `bisectory_20.pdf`.
- Comprehensive unit test suite in `test/test_spheres.jl`.

## [0.2.1] - 2026-07-23

### Added
- Hyperplane (`Plane` / `Plane2F`) geometry type with $D$-dimensional normal/point representation and point/plane/line conversion constructors.
- Point-plane and line-plane duality / polarity transformations (`polar`).
- Cairo extension integration (`CairoExt`) for rendering point sequences (`cairo_draw_points`) and polygons (`cairo_draw_polygon`).
- Example `polarity.jl` demonstrating point-plane and line-plane duality round-trip operations.
- Example `polar_bisectory.jl` demonstrating 2D convex hull generation in the disk $r \le 0.5$, origin drawing, supporting line polarity, vertex polar lines, and 500-sample perimeter bisector line polar point visualization across 4 output PDF pages.

### Fixed
- Cleared Cairo path buffer state (`Cairo.new_path`) prior to rendering points and polygon edges in `CairoExt.jl` to prevent stray line artifacts between point sets and polygons.
- Properly scoped Cairo line width settings in user/device coordinate transformations when rendering PDF/PNG outputs.

## [0.2.0] - 2026-04-29

### Added
- Subspace diameter algorithms: `approx_diameter_subspace` and `exact_diameter_subspace`.
- `MVBB` (Minimum Volume Bounding Box) module for 3D point sets.
- Rotating Calipers algorithm for 2D minimum area rectangle (`MVBB.rotating_calipers_min_area`).
- `approx_diam` in `MVBB` using BBT-based spatial decomposition.
- `OBBox` type for arbitrarily oriented bounding boxes.

### Changed
- Exported subspace diameter functions in the main `BasicCompGeometry` module.
- Improved `BBT` (Bounding Box Tree) with lazy expansion and optimized node refinement.

## [0.1.0] - 2026-04-28

### Added
- Core geometric primitives: `Point`, `Segment`, `Line`, and `BBox` (Bounding Box) in arbitrary dimensions.
- `PntSeq` (Point Sequence) and `MatPntSeq` (Matrix-backed sequence) types.
- 2D and 3D Convex Hull algorithms.
- Well-Separated Pairs Decomposition (WSPD).
- Exact and approximate ($1+\epsilon$) Euclidean diameter algorithms.
- Hausdorff distance-based curve simplification (`hausdorff_simplify`).
- Uniform arc-length resampling of point sequences.
- Bounding Box Tree (BBT) for efficient spatial decomposition.
- Comprehensive documentation using `Documenter.jl`.
- GitHub Actions for automated documentation deployment.

### Changed
- Refactored `Polygon` to `PntSeq` to better represent general sequences of points.
- Renamed internal files to match the new `PntSeq` naming convention.
- Modernized docstrings and improved LaTeX compatibility.
