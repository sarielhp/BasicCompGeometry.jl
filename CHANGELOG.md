# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
