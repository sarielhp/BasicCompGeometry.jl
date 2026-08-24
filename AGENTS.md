# AGENTS.md

Instructions, architecture guide, and workflow conventions for AI agents and developers working on `BasicCompGeometry.jl`.

---

## 1. Project Overview

**Boundary**: All work is confined to this project directory. Do not edit files outside `/home/sariel/prog/26/BasicCompGeometry/` without explicit permission.

`BasicCompGeometry.jl` is a high-performance Julia library providing fundamental computational geometry data structures and algorithms in arbitrary dimensions:

- **Geometric Primitives**: [`Point`](src/Points.jl), [`Segment`](src/Segments.jl), [`Line`](src/Segments.jl), [`Plane`](src/Planes.jl), [`BBox`](src/BBoxes.jl), [`PntSeq`](src/PntSeqs.jl), [`Sphere`](src/Spheres.jl), [`CircleArc`](src/Spheres.jl), [`CurvePolygon2D`](src/CurvePolygon.jl).
- **Algorithms**:
  - 2D & 3D Convex Hulls ([`src/ConvexHull.jl`](src/ConvexHull.jl), [`src/ConvexHull3D.jl`](src/ConvexHull3D.jl))
  - Minimum Volume Bounding Box ([`src/MVBB.jl`](src/MVBB.jl))
  - Bounding Box Trees / Spatial Indexing ([`src/BBT.jl`](src/BBT.jl))
  - Well-Separated Pair Decomposition ([`src/WSPD.jl`](src/WSPD.jl))
  - Exact & Approximate Diameter ([`src/Diameter.jl`](src/Diameter.jl))
  - Metric Spaces & Embeddings ([`src/MetricSpace.jl`](src/MetricSpace.jl))
  - Longest Convex Subset ([`src/LongestConvexSubset.jl`](src/LongestConvexSubset.jl))

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

- **"bump"**: Run `ruby scripts/bump_version` to bump the patch version, commit, push, and tag on GitHub.

