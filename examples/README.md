# Examples

This directory contains runnable example scripts demonstrating the
capabilities of `BasicCompGeometry.jl`. Each example is self-contained
and can be run with Julia 1.11 or newer (some entry points use `@main`).
Prepare the environment against this checkout once:

```bash
julia --startup-file=no --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
```

Then run an example with:

```bash
julia --project=examples examples/<script>.jl [args...]
```

For reproducible runs in this environment, set `QUICKENV_DISABLE_AUTO=true`.
Otherwise scripts importing QuickEnv may select its automatic shared environment.
Figure examples require Ipe, a LaTeX installation, `pdftocairo`, and `qpdf`;
animated examples also require `ffmpeg`. The Julia library itself supports Julia 1.9.

## Automated smoke tests

`test/example_smoke.yml` records what every program does, its small-run arguments,
expected output files, and PDF page counts. `test/example_smoke.rb` runs each
program in a fresh Julia process, with a 240-second timeout and its own output
directory under `output/tmp/example-smoke/`. An inventory check fails if a new
example has no smoke case. All three `prophet.jl` modes have separate cases.

```bash
ruby test/example_smoke.rb
# Collect library and extension coverage as well:
BCG_EXAMPLE_COVERAGE=true ruby test/example_smoke.rb
# Run one case while fixing it:
ruby test/example_smoke.rb --name test_example_sq_vs_point
```

The runner enables `BCG_EXAMPLE_SMOKE=true` to reduce benchmark sizes, sample
counts, heatmap grids, and animation frames while retaining the main computation
and rendering paths. It checks numerical approximation bounds and area/polarity
invariants, validates PDFs and page counts, parses Ipe/SVG XML, and checks PNG
signatures and GIF frame counts. It does not assess visual layout or exhaustively
test every command-line option; the polar-bisector movie modes are not included.
CI merges these runs and the external-tool integration tests into Codecov.

To run a small example directly:

```bash
QUICKENV_DISABLE_AUTO=true BCG_EXAMPLE_SMOKE=true \
  julia --project=examples examples/sq_vs_point.jl
```

`BCG_EXAMPLE_OUTPUT_DIR` can redirect generated files into a chosen output
subdirectory. Ordinary runs keep the original input sizes and animation lengths.

Entries are sorted by modification date (newest first).

Starting with v0.5.0, the default API contains 130 names; 89 additional public
names are available through qualification, explicit imports, or
`using BasicCompGeometry.Extras`. Examples explicitly import the shape accessors,
specialized helpers, and Cairo/LaTeX utilities they use. The `ann.jl` example uses
`depth` from its existing `BasicCompGeometry.BBT` import. Figure scripts still
load Cairo, Colors, or LaTeXStrings as required by the relevant extensions.

For examples that need the default API and all additional public helpers, use:

```julia
import BasicCompGeometry as CG
using .CG, .CG.Extras
```

---

## `convex_hull_figure.jl` — AI-Generated Standalone Convex Hull & Diameter Figure

Demonstrates a standalone figure script using `QuickEnv`, `BasicCompGeometry`, and `IpeDraw`.
Generates 35 deterministic 2D points, computes their 2D convex hull, finds the diameter pair and segment,
and renders a publication-ready Ipe/PDF figure with LaTeX formulas, a named theme, legend, scale bar,
and a zoomed-in inset viewport. Can be run completely standalone without `--project` flags via `QuickEnv`.

**Usage:**
```bash
./examples/convex_hull_figure.jl
# or:
julia examples/convex_hull_figure.jl
```

**Output:** `output/convex_hull_figure.ipe`, `output/convex_hull_figure.pdf`

---

## `disks_figure.jl` — Compact Geometry-Aware Figure

Uses the compact drawing API, a named theme, legend, scale bar, and clipped
inset to draw two intersecting disks with automatic world-to-page fitting. It
writes an editable Ipe source and a publication PDF. Pass `--preview` to open
the result in Ipe after generation.

**Usage:** `julia --project=. examples/disks_figure.jl`

**Output:** `output/disks_figure.ipe`, `output/disks_figure.pdf`

---

## `draw_halfplanes_arrangement.jl` — Arrangement Faces & Depth Heatmap Visualization

Reads a saved set of halfplanes (or generates random ones), computes their arrangement
and default view bounding box, reconstructs the convex polygonal faces inside the view
from chains of adjacent vertical trapezoids, and renders a 2-page publication-quality PDF:
- **Page 1**: All arrangement faces filled with distinct soft colors and stroking the boundary lines.
- **Page 2**: Arrangement depth heatmap (colored by face depth using the Viridis palette) with inward directional tick marks ("whiskers") along every halfplane boundary line.

**Usage:**
```bash
./examples/draw_halfplanes_arrangement.jl [--input file.txt] [--output file.pdf] [--n 12]
# Examples:
./examples/draw_halfplanes_arrangement.jl --input output/pruned_halfplanes_k5.txt
./examples/draw_halfplanes_arrangement.jl --n 10 --output output/random_arrangement.pdf
```

---

## `prune_halfplanes_depth.jl` — Minimal Halfplane Subset Preserving Arrangement Center Depth $\ge k$

Receives a target depth parameter $k$ and generates $6k$ random halfplanes (or reads an
initial set from a text file). Computes the arrangement's default view and vertical decomposition,
extracts the centroid of every vertical trapezoid, and greedily prunes halfplanes as long
as the depth at every center remains $\ge k$. Stops when no further halfplane can be legally
removed, outputs the final set size, and saves the minimal halfplane set to a text file.
Can read saved text files back as the initial halfplane set.

**Usage:**
```bash
./examples/prune_halfplanes_depth.jl [k] [--input file.txt] [--output file.txt]
# Examples:
./examples/prune_halfplanes_depth.jl 5
./examples/prune_halfplanes_depth.jl 4 --input output/pruned_halfplanes_k5.txt
```

---

## `rand_halfplanes_depth.jl` — Random Halfplanes Arrangement & Vertical Decomposition Depth

Generates $n$ (default 100) random halfplanes, computes their default view rectangle,
constructs the vertical decomposition of their boundary lines inside the view,
computes the centroid of each vertical trapezoid, evaluates its depth (the number
of halfplanes containing the center), and determines the total depth of the arrangement
(the minimum depth across all trapezoid centers).

**Usage:**
```bash
./examples/rand_halfplanes_depth.jl [n] [seed] [--plot]
# Or:
julia --project=. examples/rand_halfplanes_depth.jl 100 42
```

---

## `ipe_conceptual_figure.jl` — Conceptual Ipe Vector Figure Generation

Demonstrates programmatic generation of publication-ready, vector conceptual figures
with native LaTeX math formulas using `BasicCompGeometry.IpeDraw` and `LaTeXStrings`.
Automatically produces:
- An editable Ipe 7 XML document (`output/inverse_estimation_demo.ipe`)
- A cropped vector PDF (`output/inverse_estimation_demo.pdf`)
- A companion LaTeX figure fragment (`output/inverse_estimation_demo_fig.tex`) ready for `\input`.

**Usage:**
```bash
./examples/ipe_conceptual_figure.jl
```


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
is a parabola. Renders a 10-page PDF document:
- **Page 1**: Interior Voronoi cell for $p = (0.30, 0.40)$ (4 parabolic arcs).
- **Page 2**: 3D surface plot of $f(p) = \text{area}(\text{Vor}(p))$ across $[0.01, 0.99]^2$.
- **Page 3**: Superlevel set region where $f(p) > 1/6$.
- **Page 4**: Corner Voronoi cell for $p = (0.05, 0.05)$ with magnified inset showing the 2-arc parabolic lens.
- **Page 5**: Mathematical analysis and explanation of Page 4.
- **Page 6**: Edge-adjacent Voronoi cell for $p = (0.50, 0.10)$ (4 parabolic arcs).
- **Page 7**: Mathematical analysis and explanation of Page 6.
- **Page 8**: Edge-adjacent Voronoi cell for $p = (0.50, 0.02)$ (4 parabolic arcs).
- **Page 9**: Mathematical analysis and explanation of Page 8.
- **Page 10**: High-resolution heat map of $f(p) = \text{area}(\text{Vor}(p))$ across the unit square $[0, 1]^2$ with a Viridis color scale, colorbar, and key area reference values.

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
