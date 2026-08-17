# AGENTS.md

Instructions, architecture guide, and workflow conventions for AI agents and developers working on `BasicCompGeometry.jl`.

---

## 1. Project Overview

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
Defined in [`ext/CairoExt.jl`](ext/CairoExt.jl). Provides core 2D vector canvas transformations and primitives:
- `cairo_draw_setup(cr, bb, cw, ch, margin=20)`: Auto-scales and flips the y-axis (y-upwards mathematical orientation) to fit a bounding box onto a canvas.
- `cairo_draw_points(cr, points, radius=2)`: Draws point sequences as filled circles.
- `cairo_draw_polygon(cr, poly, close=true)`: Draws and strokes polygon boundaries.

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
- **External Dependencies**:
  - `lualatex` (or `pdflatex`): for compiling LaTeX text and snippets.
  - `pdftocairo`: for rasterizing in-canvas LaTeX snippet overlays.
  - `qpdf`: for merging and interleaving multi-page PDFs.
