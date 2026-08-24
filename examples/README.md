# Examples

This directory contains runnable example scripts demonstrating the
capabilities of `BasicCompGeometry.jl`. Each example is self-contained
and can be run with:

```bash
julia --project=examples examples/<script>.jl [args...]
```

Entries are sorted by modification date (newest first).

---

## `misc/monte_carlo_area.jl` — Monte Carlo Area Estimation of Curved Voronoi Cell

Constructs the Voronoi cell of a query point with respect to the 4 boundary lines
of the unit square $[0, 1]^2$ as a `CurvePolygon2D` formed by 4 `ParabolaArc` curves.
Estimates the cell area using Monte Carlo sampling across $N = 10^4, 10^5, 10^6, 10^7$
points using the ray-casting `is_inside` point-in-curved-polygon test and compares
convergence against the boundary discretization reference area.

**Usage:**
```bash
julia --project=examples examples/misc/monte_carlo_area.jl
```

---

## `polarity.jl` — Point-Plane Polarity

Demonstrates the involutive point-plane polarity with respect to the
unit sphere: `polar(polar(x)) = x`. In 2D, a plane is a line, so
`polar(Line)` works the same as `polar(Plane)`. Verifies the round-trip
property numerically.

---

## `prophet.jl` — Prophet Voronoi Diagram (Periodic 3×3 Grid)

Computes Voronoi diagrams on the flat torus (periodic 3×3 replication)
and identifies the "prophet" site with the largest cell area. Three
modes:

- **Mode 1** (prefix): Renders Voronoi diagrams for prefixes
  `{p_1, ..., p_k}` from k=4 to N, tracking the prophet cell.
- **Mode 2** (powers of 2): Compares N = 2^3 … 2^12 with calibrated
  styling (ink area ≤ 25%).
- **Mode 3** (default, 9 pages): Final diagram, prefix diagram, text
  pages with LaTeX, area-threshold highlighting, and vertex analysis.

**Usage:**
```
julia --project=examples examples/prophet.jl          # Mode 3 (N=1000)
julia --project=examples examples/prophet.jl 1 [N]    # Mode 1
julia --project=examples examples/prophet.jl 2        # Mode 2
julia --project=examples examples/prophet.jl 3 [N]    # Mode 3
```

**Output:** `output/prophet.pdf`

---

## `vor_flower.jl` — Voronoi Flower

Generates a 10-page interleaved PDF visualization exploring concentric
scaled regular polygon sites and their Voronoi diagrams, with summary
pages alternating between figures and text.

**Output:** `output/flower.pdf`

---

## `bisector_cell.jl` — Bisector Cell of Tangent Lines

Generates 10 random lines tangent to the unit square, computes their
bisector curves (parabolas) with the origin, and constructs the cell in
their arrangement that contains the origin. Outputs a single-page PDF.

**Output:** `output/bisector_cell.pdf`

---

## `sq_vs_point.jl` — Voronoi Cell of a Point vs. the Unit Square

Computes the Voronoi cell of a query point with respect to the four
boundary lines of the unit square. The bisector of a point and a line
is a parabola. Renders a 9-page PDF document:
- **Page 1**: Interior Voronoi cell for $p = (0.30, 0.40)$ (4 parabolic arcs).
- **Page 2**: 3D surface plot of $f(p) = \text{area}(\text{Vor}(p))$ across $[0.01, 0.99]^2$.
- **Page 3**: Superlevel set region where $f(p) > 1/6$.
- **Page 4**: Corner Voronoi cell for $p = (0.05, 0.05)$ with magnified inset showing the 2-arc parabolic lens.
- **Page 5**: Mathematical analysis and explanation of Page 4.
- **Page 6**: Edge-adjacent Voronoi cell for $p = (0.50, 0.10)$ (4 parabolic arcs).
- **Page 7**: Mathematical analysis and explanation of Page 6.
- **Page 8**: Edge-adjacent Voronoi cell for $p = (0.50, 0.02)$ (4 parabolic arcs).
- **Page 9**: Mathematical analysis and explanation of Page 8.

Also generates an animated GIF showing the continuous deformation of the Voronoi cell as the site moves vertically from $(0.50, 0.50)$ down to $(0.50, 0.001)$ (60 frames at 6 fps, with the final frame held for 10 frames).

**Output:** `output/sq_vs_point.pdf`, `output/sq_vs_point.gif`

---

## `sphere_diameter.jl` — Diameter on the Sphere

Computes the exact and approximate diameter of random point sets on
the sphere in various dimensions. Demonstrates BBT, WSPD, and diameter
approximation algorithms. Outputs a table of results.

---

## `polar_bisectory.jl` — Polar Duality and Bisector Arcs

A comprehensive 10-page PDF (plus animation frames) demonstrating:
- Convex hull computation
- Polar duality with respect to the unit circle
- Bisector inversion arcs (CircleArc from inverted hull edges)
- `CurvePolygon2D` construction from circular arcs
- Tangent line intersection functions (`before_tangent_to_polygon`,
  `after_tangent_to_polygon`)
- Periodic angle plots of intersection parameters
- Staircase algorithm with horizontal ray-shooting

**Outputs:**
- `output/bisectory_random.pdf`
- `output/bisectory_square.pdf`
- `output/bisectory_20.pdf`
- `output/frames/` (200 PNG frames)
- `output/staircase.mp4`

---

## `points_and_bins.jl` — Points, Grid Cells, and Voronoi Cells

A 15-page PDF demonstration:
1. Empty 10×10 grid
2. N=300 random points with a singleton cell highlighted
3–12. Progressive zoom into the singleton cell
13. Curved polygon (bisector cell against grid cell boundary)
14. Zoomed-out view with neighbors
15. Voronoi cell (red) and curved polygon (green) overlaid

**Output:** `output/points_and_bins.pdf`

---

## `longest_convex_subset.jl` — Longest Convex Subset

Finds the largest subset of random points in the plane that form a
convex polygon, using the `LongestConvexSubset` module. Renders the
result to a PDF.

**Output:** `output/longest_convex_subset.pdf`

---

## `ann.jl` — Nearest Neighbor Search Comparison

Compares four nearest neighbor search strategies on random point sets:
exact linear scan, approximate (`approx_nn`), heuristic ("silly"), and
hybrid. Reports approximation quality (ratios) and running time.

**Usage:** `julia --project=examples examples/ann.jl <n> <d> <q> [c]`
- `n`: number of points to store in the BBT
- `d`: dimension
- `q`: number of query points
- `c`: approximation parameter for `approx_nn` (default: 1.0)

---

## `dev/` — Development / Verification Scripts

Scripts in `examples/dev/` are used for internal verification and
benchmarking of subspace diameter approximation algorithms:

| Script | Purpose |
|---|---|
| `verify_subspace_diam.jl` | Compares WSPD-based and FST-based `approx_diam_subspace` against exact naive (100 tests, n=100) |
| `verify_multi_scale.jl` | Multi-scale verification across varying dimensions and epsilons |
| `compare_subspace_diam.jl` | Compares implementations at scale (100 tests, n=10000) |
| `compare_subspace_diam_batch.jl` | Batch comparison for a single (n, epsilon) pair |
| `benchmark_subspace_diam.jl` | High-scale benchmarking across multiple n values |