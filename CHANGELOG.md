# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
