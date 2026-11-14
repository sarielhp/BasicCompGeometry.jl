# AGENTS.md

Instructions, architecture guide, and workflow conventions for AI agents and developers working on `BasicCompGeometry.jl`.

---

## 1. Project Overview

**Boundary**: All work is confined to this project directory. Do not edit files outside `/home/sariel/prog/26/BasicCompGeometry/` without explicit permission.

`BasicCompGeometry.jl` is a high-performance Julia library providing fundamental computational geometry data structures and algorithms in arbitrary dimensions:

- **Geometric Primitives**: [`Point`](src/Points.jl), [`Segment`](src/Segments.jl), [`Line`](src/Segments.jl), [`Plane`](src/Planes.jl), [`BBox`](src/BBoxes.jl), [`PntSeq`](src/PntSeqs.jl), [`Sphere`](src/Spheres.jl), [`CircleArc`](src/Spheres.jl), [`Hyperbola`](src/Hyperbolas.jl), [`Parabola`](src/Parabolas.jl), [`ParabolaArc`](src/Parabolas.jl), [`CurvePolygon2D`](src/CurvePolygon.jl).
- **Algorithms**:
  - 2D & 3D Convex Hulls ([`src/ConvexHull.jl`](src/ConvexHull.jl), [`src/ConvexHull3D.jl`](src/ConvexHull3D.jl))
  - Minimum Volume Bounding Box ([`src/MVBB.jl`](src/MVBB.jl))
  - Bounding Box Trees / Spatial Indexing ([`src/BBT.jl`](src/BBT.jl))
  - Well-Separated Pair Decomposition ([`src/WSPD.jl`](src/WSPD.jl))
  - Exact & Approximate Diameter ([`src/Diameter.jl`](src/Diameter.jl))
  - Metric Spaces & Embeddings ([`src/MetricSpace.jl`](src/MetricSpace.jl))
  - Longest Convex Subset ([`src/LongestConvexSubset.jl`](src/LongestConvexSubset.jl))
  - Curve simplification hierarchies and palettes ([`src/PntSeqHierarchy.jl`](src/PntSeqHierarchy.jl))
  - Output-sensitive discrete (lattice) convex hulls of convex bodies ([`src/DiscreteHull.jl`](src/DiscreteHull.jl))

---

## 2. Package Extensions

The library uses Julia's package extension mechanism (`[weakdeps]` and `[extensions]`) to keep the core dependency footprint lightweight while providing rich visualization and documentation capabilities.

### A. `CairoExt` (Activated via `using Cairo`)
Defined in [`ext/CairoExt.jl`](ext/CairoExt.jl). Provides core 2D vector canvas transformations, primitives, and unified multi-page/animated output:
- `Canvas(path, cw, ch; fps=20, title=nothing)`: Unified canvas supporting `.pdf` (multi-page PDF), `.svg` (multi-page SVG sequence), `.png` (multi-frame PNG sequence), `.gif` (animated GIF via `ffmpeg`), and `.html` (interactive slide deck presentation directory with SVG slides and `index.html`).
- `open_canvas(f, path, cw, ch; fps=20, title=nothing)`: Opens canvas, passes to `f(canvas)`, and automatically finalizes on exit.
- `description(cr_or_canvas, text)`: Attaches a text caption to the active slide/page (rendered beneath the figure in HTML output).
- `get_file_path(canvas_or_path)`: Returns the absolute file path to the generated document (e.g. `.../index.html` for HTML canvas).
- `cairo_draw_setup(cr_or_canvas, bb, cw, ch, margin=20)`: Auto-scales and flips the y-axis (y-upwards mathematical orientation) to fit a bounding box onto a canvas.
- `cairo_draw_points(cr_or_canvas, points, radius=2)`: Draws point sequences as filled circles.
- `cairo_draw_polygon(cr_or_canvas, poly, close=true)`: Draws and strokes polygon boundaries.

### B. `BBTCairoExt` (Activated via `using Cairo, Colors`)
Defined in [`ext/BBTCairoExt.jl`](ext/BBTCairoExt.jl).
- `BBT.Tree_draw(tree, filename)`: Produces a multi-page PDF rendering level-by-level bounding boxes of a 2D Bounding Box Tree.

### C. `LaTeXCairoExt` (Activated via `using Cairo, LaTeXStrings`)
Defined in [`ext/LaTeXCairoExt.jl`](ext/LaTeXCairoExt.jl). Provides first-class LaTeX math typesetting, in-canvas annotations, snippet extraction, and PDF merging.

---

## 3. Adding LaTeX Annotations to Cairo Images

The `LaTeXCairoExt` extension makes it seamless to annotate geometric figures with publication-quality math formulas and labels:

### 1. Register Session Macros & Packages
Define common notation once at the top of a script:
```julia
using BasicCompGeometry, Cairo, LaTeXStrings

add_latex_packages!("amsmath", "amssymb", "amsfonts", "microtype")
add_latex_macros!(raw"""
\newcommand{\R}{\mathbb{R}}
\newcommand{\eps}{\varepsilon}
\DeclareMathOperator{\Vor}{Vor}
""")
```

### 2. In-Canvas Math Annotations (`cairo_draw_latex`)
Render mathematical labels directly onto the Cairo canvas:
```julia
# Draw math label at (x, y) with customizable alignment
cairo_draw_latex(cr, x_pos, y_pos, L"\Vor(p_i) \cap [0, 1]^2";
                 fontsize = 13.0,
                 halign = :left,      # :left, :center, :right
                 valign = :baseline)  # :baseline, :top, :center, :bottom
```
- **Caching**: Compiled LaTeX snippets are automatically cached in memory, ensuring fast redrawing across many frames or sites.

### 3. Extracting Macros from Existing `.tex` Files (`read_latex_snippet` / `append_latex_preamble!`)
Extract macros delimited by user-defined marker comments from existing LaTeX source files:
```julia
# Extract snippet between custom markers and append directly to session preamble:
append_latex_preamble!("macros.tex"; beg_marker="%%% MACROS BEG", end_marker="%%% MACROS END")

# Or read a snippet as a string:
snip = read_latex_snippet("figure.tex"; beg_marker="%%% PREAMBLE START", end_marker="%%% PREAMBLE END")
```

### 4. Generating Text Pages & Merging Multi-Page Documents
Combine Cairo diagram pages with LaTeX documentation pages:
```julia
# 1. Compile LaTeX text pages to PDF matching Cairo dimensions (bp)
latex_to_pdf(["\\section*{Page 1}\nExplanation...", "\\section*{Page 2}\nFormulas..."],
             "text_pages.pdf"; paperwidth=800, paperheight=800, margin=50)

# 2. Interleave diagram pages and text pages with qpdf
pdf_merge("final_output.pdf",
          ("diagrams.pdf", 1), ("text_pages.pdf", 1),
          ("diagrams.pdf", 2), ("text_pages.pdf", 2))
```

---

## 4. Development & Testing Guidelines

- **Test Re-run Policy**: Tests only need to be rerun when changes are made to the core library (`src/`) or test suite (`test/`). When only examples (`examples/`) are updated or added, do not rerun the full test suite.
- **Running Test Suite**:
  ```bash
  julia --project=. test/runtests.jl
  ```
- **Gate** (`tools/gate`, the standard Julia gate from `~/prog/standards/julia/templates`): `tools/lint` (strict `julia-audit` against `.julia-audit-baseline`, syntax, JuliaFormatter check), the temp-file policy, `tools/test` (`Pkg.test()`), and `tools/build` (precompile). The baseline lists known violations; new ones fail, and fixed ones must be removed from it. Format with JuliaFormatter (`.JuliaFormatter.toml`) before committing.
- **Golden outputs** ([`test/test_golden.jl`](test/test_golden.jl), part of `Pkg.test`): Ipe, binary and PLT outputs are compared byte for byte with [`test/golden/`](test/golden/). An intended output change means regenerating with `BCG_UPDATE_GOLDEN=1` and reviewing the diff; never edit a fixture to make a test pass.
- **Full Gate** (`tools/gate --full`): the gate plus property/reference tests on seeded random inputs ([`test/properties/properties.jl`](test/properties/properties.jl)), allocation baselines of the hot entry points ([`test/alloc/alloc_check.jl`](test/alloc/alloc_check.jl), baseline per Julia minor version, `--update` only after reviewing why allocation changed), JET `report_package` and ExplicitImports, and a JET runtime-dispatch check of the hot entry points ([`test/jet/jet_opt.jl`](test/jet/jet_opt.jl), its own environment in `test/jet/`). `tools/release` and review-cycle run it. When adding a performance-critical public function, add an entry point to that file; a failure means an abstract field or container type, or a type-unstable local, makes the code dispatch at run time.
- **Running Extension Tests**:
  ```bash
  julia --project=examples test/test_latex_cairo.jl
  ```
- **Running Examples**:
  ```bash
  julia --project=examples examples/prophet.jl [mode] [args...]
  ```

---

## 5. Coding Conventions

- **2D Point Coordinate Access**: For points in two dimensions, use `.x` and `.y` (e.g., `p.x`, `p.y`) to access coordinates instead of indexing (`p[1]`, `p[2]`), unless it complicates the code (e.g., in generic dimension-agnostic algorithms).
- **Line Widths & Radii across Vector (PDF) vs. Raster (PNG / GIF)**:
  - When rendering with `scale = cairo_draw_setup(cr_or_canvas, bb, cw, ch, margin)`:
    - **Stroke Line Widths**: Always specify stroke widths in **absolute device point/pixel sizes** (e.g. `Cairo.set_line_width(canvas, 1.5)` or `3.0`), rather than tiny user-coordinate bounding-box units (e.g. `0.005`). On raster surfaces (`.png`, `.gif` via `CairoImageSurface`), sub-pixel user-coordinate values get rounded to zero by Cairo's rasterizer, making borders and curves completely invisible even though vector PDF renderers might draw them.
    - **Point Radii**: For point markers drawn with `Cairo.arc(canvas, x, y, r, ...)`, scale the radius using `scale` returned from `cairo_draw_setup` (e.g. `r = 4.0 / scale` for a 4px dot), or use `cairo_draw_points(canvas, points, radius_px)`.

---

## 6. Temporary Files

Temporary files (test outputs, intermediate data, scratch plots, etc.) should be created under `output/tmp/`. The project root must stay clean — do not write files directly into `/` or any top-level directory besides `output/` and its subdirectories.

- **Creating temp files**: Use `mkpath(joinpath(@__DIR__, "..", "output", "tmp"))` then write inside it.
- **Cleanup**: If feasible, remove temp files at the end of the script.

---

## 7. Git Workflow

- **Auto-track new files**: After creating any new source file (example, script, extension, test, etc.), immediately stage and commit it. Do not leave new files untracked.
- **Output directory**: The `output/` directory is gitignored. All examples must write their output (PDFs, PNGs, SVGs, etc.) to `output/` or a subdirectory thereof.
- **Commit discipline**: Keep commits focused. Stage only the files relevant to the change.
- **Example documentation**: [`examples/README.md`](examples/README.md) describes every example program in detail. When adding or modifying an example, update this file accordingly.

- **External Dependencies**:
  - `lualatex` (or `pdflatex`): for compiling LaTeX text and snippets.
  - `pdftocairo`: for rasterizing in-canvas LaTeX snippet overlays.
  - `qpdf`: for merging and interleaving multi-page PDFs.

---

## 8. Keyword Triggers

- **"release"** / **"bump"**: Releases are frequent and do not register the package. Run `tools/release [--minor|--major]` from a clean `main` that matches GitHub: it runs the full gate (`tools/gate --full`), bumps the version (patch by default, via `tools/bump`), commits, and pushes through `tools/push`. It does not register. Use `--dry-run` to preview the next version.
- **"publish"**: Registration is rare. Run `tools/publish` only after major changes or important breaking changes, and only when explicitly asked; ordinary releases (`tools/release`) are not registered. When asked, run `tools/publish [release notes]` after a release: it waits for release CI on the pushed HEAD, then requests JuliaRegistrator registration (notes default to "See CHANGELOG.md."). TagBot creates the tag and GitHub release after registration merges. JuliaRegistrator must be enabled for this repository. Use `--dry-run` to preview.
- **Pushing**: Use `tools/push` to push the current branch to the configured `origin` and verify its remote commit. It never forces updates or pushes tags. Use `--dry-run` to preview. SSH remotes still require the explicit authorization specified by the global instructions.

## Autonomous review cycle (`review-cycle`)

The engine and review-audit live in ~/prog/26/ruby/review_cycle and run from PATH. The repository contract is ~/prog/standards/workflow/README.md; project choices belong in tools/review_cycle.json. Read project/language guidelines, agent instructions, and tooling settings before remediation.

| Command | Effect |
|---|---|
| `rtk proxy review-cycle --doctor` | Check host and configured sandbox prerequisites. |
| `rtk proxy review-cycle --setup` | Requested onboarding: configure local bws profiles/tools and project review files. |
| `rtk proxy review-cycle --no-push` | Requested audit/remediation with local commits, gate/build, and differential verification. |
| `rtk proxy review-cycle --loop --no-push` | Bounded loop across lenses; default maximum seven cycles. |
| `rtk proxy review-cycle --push` | Authorized upstream publication with remote SHA verification. |
| `rtk proxy review-cycle --status` | Inspect rotation and history. |
| `rtk proxy review-cycle --dry-run --no-push` | Audit/report without remediation or commits; writes reports/logs and runs pre-flight checks. |

Declare executable tools/gate and tools/build. Adopt all language helpers, entry points, and settings. The canonical gate owns analyzers, tests, and one build; tools/check delegates to it. Build validates repository artifacts and never installs or changes Git. Compiler/dependency caches may be external. Custom gates need a checked separate build unless gate_includes_build:true declares they include it.

The pipeline audits a lens, dispatches sandbox remediation, inspects its diff, runs the post-merge gate/build and differential audit, then writes summary/archive/state artifacts. Invalid configuration, missing commands, failed verification, and failed Git operations stop the run. Remediation or verification failures roll back to the cycle base with pre-existing untracked-file preservation. Ignored output/cache changes and external effects remain outside rollback.

The differential audit sees only the code patch (reports and CHANGELOG excluded). Every finding located in a changed file, or without a parseable location, blocks completion; findings about untouched files are logged as advisory. --skip-diff-audit records local needs_review, exits nonzero, and prevents publication. Engine evidence overrides agent-authored verification claims in summaries. Preserve regression coverage and add appropriate defect tests. Missing test edits remain advisory; --strict-tests rejects detected test weakening and --strict-diff rejects oversized diffs. Neither proves complete coverage.

Commits remain local by default; inspect existing explicit push:true settings and use --no-push when needed. Authorized publication requires a remote branch upstream and one push URL, pushes the exact intended commit, and verifies its destination SHA. Publication failures exit nonzero and retain verified local commits; inspect the remote before retrying. Batched pushes follow the same policy.

Review commits do not contain lifecycle structured verification or create release tags. Run the full project gate before authorized publication or use tools/bump/tools/release. Install only through the explicitly authorized project installation workflow.

Command approval permits execution, not unrelated work. Reuse session authorization. Keep sandboxing and manual approval review enabled; never enlarge writable roots, modify approval rules, install dependencies, or publish merely because these instructions are present.
