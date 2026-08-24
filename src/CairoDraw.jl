"""
    cairo_draw_setup(cr, bb, canvas_width, canvas_height, margin)

Set up Cairo context `cr` with a scale and translate transform so that
the bounding box `bb` is mapped to the canvas of size `canvas_width` x `canvas_height`,
with `margin` pixels of padding on each side.
The y-axis is flipped so that the coordinate system matches the standard
mathematical orientation (y increases upward).
"""
function cairo_draw_setup end

"""
    cairo_draw_points(cr, points, radius)

Draw points as filled circles of the given `radius` using the Cairo context `cr`.
`points` is any iterable of `Point{2, Float64}`.
"""
function cairo_draw_points end

"""
    cairo_set_line_width(cr_or_canvas, width)

Set the stroke line width on a Cairo context or `Canvas`.
If `width < 1.0`, it is treated as a mathematical world-coordinate unit and automatically
scaled by the active transformation matrix scale factor (`scale = sqrt(|det(CTM)|)`) and rounded
up to at least 1.0 device pixel across all surfaces. If `width >= 1.0`, it is treated as an
explicit device pixel/point width.
"""
function cairo_set_line_width end

"""
    cairo_draw_polygon(cr, poly; line_width=nothing, close=true)

Draw the edges of a polygon (or point sequence) `poly` using the Cairo context `cr`.
If `line_width` is provided, `cairo_set_line_width` is called first.
If `close` is true, the polygon is closed by connecting the last vertex back to the first.
"""
function cairo_draw_polygon end

"""
    cairo_draw_latex(cr, x, y, latex_str; fontsize=12, halign=:left, valign=:baseline, dpi=300, scale=1.0)

Draw a LaTeX math or text string `latex_str` onto the Cairo context `cr` at `(x, y)`.
Requires `Cairo` and `LaTeXStrings`.
"""
function cairo_draw_latex end

"""
    cairo_draw_latex_page(canvas_or_cr, latex_str; compiler="lualatex", dpi=300)

Compile and render a full LaTeX page `latex_str` directly onto the current canvas page.
Requires `Cairo` and `LaTeXStrings`.
"""
function cairo_draw_latex_page end

"""
    latex_to_pdf(pages, output_path; paperwidth=800, paperheight=800, margin=50, packages=["amsmath", "amssymb", "amsfonts", "microtype"], compiler="lualatex", passes=2)

Compile one or more LaTeX page strings into a multi-page PDF document.
Page dimensions are specified in big points (bp).
Requires `Cairo` and `LaTeXStrings`.
"""
function latex_to_pdf end

"""
    pdf_merge(output_pdf, sources...)

Merge multiple PDF files or page ranges into `output_pdf` using `qpdf`.
Requires `Cairo` and `LaTeXStrings`.
"""
function pdf_merge end

"""
    read_latex_snippet(file_or_content; beg_marker=nothing, end_marker=nothing)

Read text from a `.tex` file or string. If `beg_marker` is specified,
extracts the portion between `beg_marker` and `end_marker`. If omitted, returns the whole text.
Requires `Cairo` and `LaTeXStrings`.
"""
function read_latex_snippet end

"""
    append_latex_preamble!(source; beg_marker=nothing, end_marker=nothing)

Append LaTeX packages or macros from a string or `.tex` file to the global session preamble.
If `source` is a file path, extracts delimited snippet if markers are present.
Requires `Cairo` and `LaTeXStrings`.
"""
function append_latex_preamble! end

"""
    set_latex_preamble!(source; beg_marker=nothing, end_marker=nothing)

Set the global session LaTeX preamble from a string or `.tex` file.
Requires `Cairo` and `LaTeXStrings`.
"""
function set_latex_preamble! end

"""
    add_latex_macros!(source; beg_marker=nothing, end_marker=nothing)

Alias for `append_latex_preamble!`.
Requires `Cairo` and `LaTeXStrings`.
"""
function add_latex_macros! end

"""
    add_latex_packages!(pkgs...)

Add one or more LaTeX package names to the global session package list.
Requires `Cairo` and `LaTeXStrings`.
"""
function add_latex_packages! end

"""
    reset_latex_preamble!()

Reset session LaTeX packages and custom macros to library defaults.
Requires `Cairo` and `LaTeXStrings`.
"""
function reset_latex_preamble! end

"""
    get_latex_preamble(; extra_packages=[], extra_preamble="")

Retrieve the combined LaTeX preamble string for the active session.
Requires `Cairo` and `LaTeXStrings`.
"""
function get_latex_preamble end

"""
    Canvas(filename, canvas_width, canvas_height; fps=20, title=nothing)

Create a unified 2D vector/raster canvas supporting `.pdf`, `.svg`, `.png`, `.gif`, and `.html`.
Allows drawing and multi-page / animated transitions (`Cairo.show_page`) with identical syntax.
For `.html`, generates an interactive presentation directory containing the SVGs and `index.html`.
Requires `Cairo`.
"""
function Canvas end

"""
    open_canvas(f::Function, filename, canvas_width, canvas_height; kwargs...)

Execute function `f(canvas)` and automatically finalize and flush with `Cairo.finish`.
Requires `Cairo`.
"""
function open_canvas end

"""
    description(canvas_or_cr, text)

Attach a descriptive text note to the current page/slide.
For HTML canvas output, this description is displayed beneath the figure in the generated interactive presentation.
For other formats, this call is a no-op unless embedded comments are supported.
"""
function description end

"""
    get_file_path(canvas_or_path)

Retrieve the full absolute path of the generated output file.
For HTML canvas output, returns the absolute path to `index.html`.
"""
function get_file_path end

export cairo_draw_setup, cairo_draw_points, cairo_draw_polygon, cairo_set_line_width, cairo_draw_latex, cairo_draw_latex_page, latex_to_pdf, pdf_merge
export read_latex_snippet, append_latex_preamble!, set_latex_preamble!, add_latex_macros!, add_latex_packages!, reset_latex_preamble!, get_latex_preamble
export Canvas, open_canvas, description, get_file_path