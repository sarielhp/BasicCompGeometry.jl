# Cairo and LaTeX

## Loading extensions

Visualization methods load when their optional packages are present in the
current environment. Install them with `Pkg.add` as needed, then load:

| Capability | Packages to load |
|:-----------|:-----------------|
| Cairo canvases and 2D drawing | `Cairo` |
| `BBT.Tree_draw` | `Cairo`, `Colors` |
| LaTeX annotations and PDF utilities | `Cairo`, `LaTeXStrings` |

Cairo and LaTeX helper names belong to `Extras`. A typical drawing script starts
with:

```julia
import BasicCompGeometry as CG
using .CG, .CG.Extras
using Cairo
```

## Canvases

`open_canvas` opens a canvas and finalizes it when the block exits. Targets can
be multi-page PDF, SVG or PNG sequences, GIF animation, or an HTML slide deck.
GIF generation requires `ffmpeg`. The HTML target creates a presentation
directory; `get_file_path` returns the path to its `index.html`.

```julia
points = PntSeq([point(0.0, 0.0), point(3.0, 0.0), point(1.0, 2.0)])
bb = bbox(points)

open_canvas("output/triangle.pdf", 600, 400) do canvas
    cairo_draw_setup(canvas, bb, 600, 400, 20)
    Cairo.set_source_rgb(canvas, 0.1, 0.2, 0.6)
    cairo_draw_polygon(canvas, points; line_width=1.5)
    cairo_draw_points(canvas, points, 3)
    description(canvas, "A triangle and its vertices")
end
```

`cairo_draw_setup` scales and translates the world box onto the canvas and flips
the y-axis to mathematical orientation. Prefer absolute device widths such as
`1.5` with `cairo_set_line_width` or `Cairo.set_line_width(canvas, 1.5)`.
`cairo_draw_points` takes radii in device units. For a direct `Cairo.arc` call,
divide a desired device radius by the scale returned by `cairo_draw_setup`.

For tree visualizations, add `using Colors` and call
`BBT.Tree_draw(tree, "output/tree.pdf")` on a 2D tree.

## Math labels

Load `LaTeXStrings` as well as `Cairo` to activate LaTeX drawing. Register
packages and macros once, then draw labels inside a canvas block:

```julia
using LaTeXStrings

add_latex_packages!("amsmath", "amssymb", "amsfonts", "microtype")
add_latex_macros!(raw"\newcommand{\R}{\mathbb{R}}")

open_canvas("output/label.pdf", 600, 400) do canvas
    cairo_draw_setup(canvas, bb, 600, 400, 20)
    cairo_draw_polygon(canvas, points; line_width=1.5)
    cairo_draw_latex(canvas, 1.0, 1.0, L"p \in \R^2";
                     fontsize=13.0, halign=:center, valign=:baseline)
end
```

LaTeX rendering requires `lualatex` (or a configured `pdflatex`) and
`pdftocairo`. Compiled snippets are cached in memory. Text is rasterized before
being placed on Cairo output; use [IpeDraw](figures.md) when labels must remain
vector content in a PDF.

## Existing preambles and PDF composition

Extract marked macros from an existing `.tex` file:

```julia
append_latex_preamble!("macros.tex";
    beg_marker="%%% MACROS BEG", end_marker="%%% MACROS END")
```

`read_latex_snippet` returns the extracted text without adding it to the
session preamble. `set_latex_preamble!` replaces that preamble, while
`reset_latex_preamble!` restores defaults.

`latex_to_pdf` compiles full text pages. `pdf_merge` uses `qpdf` to combine files
or page ranges, for example:

```julia
latex_to_pdf([raw"\section*{Explanation} A geometric construction."],
             "output/text.pdf"; paperwidth=600, paperheight=400, margin=40)
pdf_merge("output/combined.pdf", ("output/triangle.pdf", 1), ("output/text.pdf", 1))
```

See the [Cairo and LaTeX reference](reference/cairo.md) for all helper signatures.
