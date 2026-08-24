module CairoExt

using BasicCompGeometry
import Cairo
import Cairo: show_page, finish
using Printf

"""
    cairo_draw_setup(cr, bb, canvas_width, canvas_height, margin=20)

Set up Cairo context `cr` with a scale and translate transform so that
the bounding box `bb` is mapped to the canvas of size `canvas_width` x `canvas_height`,
with `margin` pixels of padding on each side.
The y-axis is flipped so that the coordinate system matches the standard
mathematical orientation (y increases upward).
"""
function BasicCompGeometry.cairo_draw_setup(
    cr, bb::BBox{2,T}, cw::Real, ch::Real, margin::Real=20
) where {T}
    bl = bottom_left(bb)
    tr = top_right(bb)
    bb_width = tr[1] - bl[1]
    bb_height = tr[2] - bl[2]
    bb_width == 0 && (bb_width = 1)
    bb_height == 0 && (bb_height = 1)

    scale = min((cw - 2margin) / bb_width, (ch - 2margin) / bb_height)

    C = Cairo
    C.reset_transform(cr)
    C.translate(cr, margin, ch - margin)
    C.scale(cr, scale, -scale)
    C.translate(cr, -bl[1], -bl[2])
    return scale
end

"""
    cairo_draw_points(cr, points, radius=2)

Draw points as filled circles of the given `radius` using the Cairo context `cr`.
`points` is any iterable of `Point{2, Float64}`.
"""
function BasicCompGeometry.cairo_draw_points(cr, points, radius::Real=2)
    C = Cairo
    m = C.get_matrix(cr)
    scale = sqrt(abs(m.xx * m.yy - m.xy * m.yx))
    effective_r = (scale > 0 && radius >= 0.5) ? (Float64(radius) / scale) : Float64(radius)
    for p in points
        C.new_path(cr)
        C.arc(cr, p[1], p[2], effective_r, 0.0, 2pi)
        C.fill(cr)
    end
    return
end

"""
    cairo_draw_polygon(cr, poly; line_width=nothing, close=true)

Draw the edges of a polygon (or point sequence) `poly` using the Cairo context `cr`.
If `line_width` is provided, `cairo_set_line_width` is called first.
If `close` is true, the polygon is closed by connecting the last vertex back to the first.
"""
function BasicCompGeometry.cairo_draw_polygon(cr, poly; line_width=nothing, close::Bool=true)
    C = Cairo
    n = length(poly)
    n == 0 && return

    !isnothing(line_width) && BasicCompGeometry.cairo_set_line_width(cr, line_width)
    C.new_path(cr)
    C.move_to(cr, poly[1][1], poly[1][2])
    for i = 2:n
        C.line_to(cr, poly[i][1], poly[i][2])
    end
    if close && n > 2
        C.close_path(cr)
    end
    C.stroke(cr)
    return
end

# --- Unified Multi-Page / Animated Canvas ---

"""
    Canvas(path::String, cw::Real, ch::Real; fps::Int=20, title=nothing)

Create a unified 2D vector / raster / web canvas supporting:
- `.pdf`: Multi-page PDF document
- `.svg`: Multi-page SVG sequence (e.g. `slide_001.svg`, `slide_002.svg`, ...)
- `.png`: Multi-frame PNG sequence (e.g. `frame_0001.png`, `frame_0002.png`, ...)
- `.gif`: High-quality animated GIF (via `ffmpeg` on finish)
- `.html`: Interactive presentation directory with SVG slides and `index.html`

Advance pages or animation frames with `Cairo.show_page(canvas)`.
Finalize and flush to disk with `Cairo.finish(canvas)`.
"""
mutable struct Canvas
    path::String
    cw::Float64
    ch::Float64
    fmt::Symbol             # :pdf, :svg, :png, :gif, :html
    page::Int
    fps::Int
    tempdir::String
    htmldir::String
    index_html_path::String
    title::String
    descriptions::Dict{Int, String}
    needs_new_surface::Bool
    surface::Cairo.CairoSurfaceBase
    cr::Cairo.CairoContext
end

BasicCompGeometry.description(c::Canvas, text::AbstractString) = (c.descriptions[c.page] = String(text); text)
BasicCompGeometry.description(cr::Cairo.CairoContext, text::AbstractString) = text
BasicCompGeometry.description(any, text::AbstractString) = text

BasicCompGeometry.get_file_path(c::Canvas) = (c.fmt == :html ? abspath(c.index_html_path) : abspath(c.path))
BasicCompGeometry.get_file_path(path::AbstractString) = abspath(path)

function _svg_page_filename(pattern::String, page::Int)
    if occursin("%", pattern)
        return Printf.format(Printf.Format(pattern), page)
    end
    base, ext = splitext(pattern)
    return "$(base)_$(lpad(page, 3, '0'))$(ext)"
end

function _png_page_filename(pattern::String, page::Int)
    if occursin("%", pattern)
        return Printf.format(Printf.Format(pattern), page)
    end
    base, ext = splitext(pattern)
    return "$(base)_$(lpad(page, 4, '0'))$(ext)"
end

function _ensure_surface!(c::Canvas)
    !c.needs_new_surface && return
    c.needs_new_surface = false
    if c.fmt == :svg
        next_file = _svg_page_filename(c.path, c.page)
        surf = Cairo.CairoSVGSurface(next_file, c.cw, c.ch)
        c.surface = surf
        c.cr = Cairo.CairoContext(surf)
    elseif c.fmt == :html
        next_file = joinpath(c.htmldir, @sprintf("page_%03d.svg", c.page))
        surf = Cairo.CairoSVGSurface(next_file, c.cw, c.ch)
        c.surface = surf
        c.cr = Cairo.CairoContext(surf)
    elseif c.fmt in (:png, :gif)
        surf = Cairo.CairoImageSurface(Int(round(c.cw)), Int(round(c.ch)), Cairo.FORMAT_ARGB32)
        c.surface = surf
        c.cr = Cairo.CairoContext(surf)
    end
    return
end

function BasicCompGeometry.Canvas(path::String, cw::Real, ch::Real; fps::Int=20, title::Union{String, Nothing}=nothing)
    ext = lowercase(splitext(path)[2])
    w = Float64(cw)
    h = Float64(ch)

    if ext == ".pdf"
        surf = Cairo.CairoPDFSurface(path, w, h)
        cr = Cairo.CairoContext(surf)
        return Canvas(path, w, h, :pdf, 1, fps, "", "", "", "", Dict{Int, String}(), false, surf, cr)
    elseif ext == ".svg"
        fmt_path = _svg_page_filename(path, 1)
        surf = Cairo.CairoSVGSurface(fmt_path, w, h)
        cr = Cairo.CairoContext(surf)
        return Canvas(path, w, h, :svg, 1, fps, "", "", "", "", Dict{Int, String}(), false, surf, cr)
    elseif ext == ".png"
        surf = Cairo.CairoImageSurface(Int(round(w)), Int(round(h)), Cairo.FORMAT_ARGB32)
        cr = Cairo.CairoContext(surf)
        return Canvas(path, w, h, :png, 1, fps, "", "", "", "", Dict{Int, String}(), false, surf, cr)
    elseif ext == ".gif"
        tmp = mktempdir()
        surf = Cairo.CairoImageSurface(Int(round(w)), Int(round(h)), Cairo.FORMAT_ARGB32)
        cr = Cairo.CairoContext(surf)
        return Canvas(path, w, h, :gif, 1, fps, tmp, "", "", "", Dict{Int, String}(), false, surf, cr)
    elseif ext in (".html", ".htm")
        title_str = isnothing(title) ? "" : title
        if basename(lowercase(path)) in ("index.html", "index.htm")
            htmldir = dirname(abspath(path))
            index_html_path = abspath(path)
            if isempty(title_str)
                title_str = basename(htmldir)
            end
        else
            base = splitext(path)[1]
            htmldir = abspath(base)
            index_html_path = joinpath(htmldir, "index.html")
            if isempty(title_str)
                title_str = basename(base)
            end
        end

        mkpath(htmldir)
        for f in readdir(htmldir)
            if endswith(lowercase(f), ".svg") || f == "index.html" || f == "page"
                rm(joinpath(htmldir, f), force=true)
            end
        end

        first_svg = joinpath(htmldir, "page_001.svg")
        surf = Cairo.CairoSVGSurface(first_svg, w, h)
        cr = Cairo.CairoContext(surf)
        return Canvas(path, w, h, :html, 1, fps, "", htmldir, index_html_path, title_str, Dict{Int, String}(), false, surf, cr)
    else
        error("Unsupported Canvas format: '$ext'. Expected .pdf, .svg, .png, .gif, or .html")
    end
end

function _embed_svg_description(svg_path::String, text::String)
    isempty(strip(text)) && return
    !isfile(svg_path) && return
    content = read(svg_path, String)
    occursin("<desc>", content) && return

    xml_clean_text = replace(text, "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")
    comment_clean_text = replace(text, "--" => "- -")

    injection = "\n<!-- $(comment_clean_text) -->\n<desc>$(xml_clean_text)</desc>\n"
    m = match(r"(<svg[^>]*>)", content)
    if m !== nothing
        pos = m.offset + length(m.match)
        new_content = content[1:pos] * injection * content[(pos + 1):end]
        write(svg_path, new_content)
    end
end

function Cairo.show_page(c::Canvas)
    if c.fmt == :pdf
        Cairo.show_page(c.cr)
        c.page += 1
    elseif c.fmt in (:svg, :html)
        Cairo.finish(c.surface)
        if c.fmt == :svg
            cur_file = _svg_page_filename(c.path, c.page)
            if haskey(c.descriptions, c.page)
                _embed_svg_description(cur_file, c.descriptions[c.page])
            end
        elseif c.fmt == :html
            cur_file = joinpath(c.htmldir, @sprintf("page_%03d.svg", c.page))
            if haskey(c.descriptions, c.page)
                _embed_svg_description(cur_file, c.descriptions[c.page])
            end
        end
        c.page += 1
        c.needs_new_surface = true
    elseif c.fmt == :png
        cur_file = _png_page_filename(c.path, c.page)
        Cairo.write_to_png(c.surface, cur_file)
        Cairo.finish(c.surface)
        c.page += 1
        c.needs_new_surface = true
    elseif c.fmt == :gif
        frame_file = joinpath(c.tempdir, @sprintf("frame_%05d.png", c.page))
        Cairo.write_to_png(c.surface, frame_file)
        Cairo.finish(c.surface)
        c.page += 1
        c.needs_new_surface = true
    end
    return c
end

function _write_html_presentation(c::Canvas)
    svg_files = filter(f -> occursin(r"^page_\d+\.svg$", f), readdir(c.htmldir))
    sort!(svg_files, by = f -> parse(Int, match(r"page_(\d+)\.svg", f).captures[1]))

    if isempty(svg_files)
        @warn "No SVG files found in $(c.htmldir) to generate HTML presentation"
        return
    end

    slides_json = "[" * join(["\"$f\"" for f in svg_files], ", ") * "]"
    desc_array = [get(c.descriptions, i, "") for i in 1:length(svg_files)]
    desc_json = "[" * join(["\"" * escape_string(d) * "\"" for d in desc_array], ", ") * "]"

    page_title = isempty(c.title) ? "Presentation" : c.title

    html_content = """
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>$(page_title)</title>
  <style>
    :root {
      --bg: #0f172a;
      --card-bg: #1e293b;
      --text: #f8fafc;
      --text-muted: #94a3b8;
      --accent: #38bdf8;
      --border: #334155;
    }
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body {
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      background: var(--bg);
      color: var(--text);
      min-height: 100vh;
      display: flex;
      flex-direction: column;
      align-items: center;
      justify-content: space-between;
    }
    header {
      width: 100%;
      padding: 12px 24px;
      background: var(--card-bg);
      border-bottom: 1px solid var(--border);
      display: flex;
      align-items: center;
      justify-content: space-between;
    }
    h1 { font-size: 1.1rem; font-weight: 600; }
    .controls { display: flex; gap: 8px; align-items: center; }
    button, select {
      background: var(--bg);
      color: var(--text);
      border: 1px solid var(--border);
      padding: 6px 14px;
      border-radius: 6px;
      cursor: pointer;
      font-size: 0.9rem;
      transition: all 0.15s ease;
    }
    button:hover:not(:disabled) {
      background: var(--accent);
      color: #0f172a;
      border-color: var(--accent);
    }
    button:disabled {
      opacity: 0.4;
      cursor: not-allowed;
    }
    .slide-container {
      flex: 1;
      display: flex;
      flex-direction: column;
      align-items: center;
      justify-content: center;
      width: 100%;
      padding: 20px;
      gap: 14px;
    }
    .slide-wrapper {
      background: white;
      border-radius: 8px;
      box-shadow: 0 10px 25px -5px rgba(0, 0, 0, 0.5);
      max-width: 90vw;
      max-height: 72vh;
      display: flex;
      align-items: center;
      justify-content: center;
      overflow: hidden;
    }
    .slide-wrapper img {
      width: 100%;
      height: 100%;
      max-height: 72vh;
      object-fit: contain;
      display: block;
    }
    .slide-description {
      background: var(--card-bg);
      border: 1px solid var(--border);
      border-left: 4px solid var(--accent);
      border-radius: 6px;
      padding: 12px 18px;
      max-width: 90vw;
      width: 100%;
      box-shadow: 0 4px 12px rgba(0, 0, 0, 0.2);
      font-size: 0.95rem;
      line-height: 1.5;
      color: #e2e8f0;
      white-space: pre-wrap;
    }
    footer {
      width: 100%;
      padding: 10px 24px;
      background: var(--card-bg);
      border-top: 1px solid var(--border);
      display: flex;
      justify-content: space-between;
      align-items: center;
      font-size: 0.85rem;
      color: var(--text-muted);
    }
    .kbd {
      background: var(--bg);
      border: 1px solid var(--border);
      border-radius: 4px;
      padding: 2px 6px;
      font-family: monospace;
      font-size: 0.8rem;
    }
  </style>
</head>
<body>
  <header>
    <h1>$(page_title)</h1>
    <div class="controls">
      <button id="prevBtn" onclick="prevSlide()">← Prev</button>
      <select id="slideSelect" onchange="jumpSlide(parseInt(this.value))"></select>
      <button id="nextBtn" onclick="nextSlide()">Next →</button>
      <button id="fullscreenBtn" onclick="toggleFullscreen()">⛶ Fullscreen</button>
    </div>
  </header>

  <main class="slide-container" id="mainContainer">
    <div class="slide-wrapper" id="slideWrapper">
      <img id="slideImg" src="" alt="Slide">
    </div>
    <div class="slide-description" id="slideDesc" style="display: none;"></div>
  </main>

  <footer>
    <span id="pageStatus">Slide 1 of $(length(svg_files))</span>
    <span>Use <span class="kbd">←</span> <span class="kbd">→</span> or <span class="kbd">Space</span> to navigate</span>
  </footer>

  <script>
    const slides = $slides_json;
    const descriptions = $desc_json;
    let currentIdx = 0;

    const img = document.getElementById('slideImg');
    const desc = document.getElementById('slideDesc');
    const select = document.getElementById('slideSelect');
    const prevBtn = document.getElementById('prevBtn');
    const nextBtn = document.getElementById('nextBtn');
    const pageStatus = document.getElementById('pageStatus');

    slides.forEach((s, idx) => {
      const opt = document.createElement('option');
      opt.value = idx;
      opt.textContent = `Slide \${idx + 1} (\${s})`;
      select.appendChild(opt);
    });

    function showSlide(idx) {
      if (idx < 0 || idx >= slides.length) return;
      currentIdx = idx;
      img.src = slides[currentIdx];
      select.value = currentIdx;
      prevBtn.disabled = currentIdx === 0;
      nextBtn.disabled = currentIdx === slides.length - 1;
      pageStatus.textContent = `Slide \${currentIdx + 1} of \${slides.length}`;

      const d = descriptions[currentIdx];
      if (d && d.trim().length > 0) {
        desc.textContent = d;
        desc.style.display = 'block';
      } else {
        desc.textContent = '';
        desc.style.display = 'none';
      }
    }

    function prevSlide() { showSlide(currentIdx - 1); }
    function nextSlide() { showSlide(currentIdx + 1); }
    function jumpSlide(idx) { showSlide(idx); }

    function toggleFullscreen() {
      if (!document.fullscreenElement) {
        document.documentElement.requestFullscreen().catch(err => {});
      } else {
        document.exitFullscreen();
      }
    }

    window.addEventListener('keydown', (e) => {
      if (e.key === 'ArrowRight' || e.key === ' ' || e.key === 'PageDown') {
        nextSlide();
      } else if (e.key === 'ArrowLeft' || e.key === 'PageUp') {
        prevSlide();
      } else if (e.key === 'Home') {
        showSlide(0);
      } else if (e.key === 'End') {
        showSlide(slides.length - 1);
      } else if (e.key.toLowerCase() === 'f') {
        toggleFullscreen();
      }
    });

    showSlide(0);
  </script>
</body>
</html>
"""
    write(c.index_html_path, html_content)
end

function Cairo.finish(c::Canvas)
    if c.fmt == :pdf
        Cairo.finish(c.surface)
    elseif c.fmt == :svg
        if !c.needs_new_surface
            Cairo.finish(c.surface)
        end
        total_pages = c.needs_new_surface ? c.page - 1 : c.page
        if !c.needs_new_surface && haskey(c.descriptions, total_pages)
            cur_file = _svg_page_filename(c.path, total_pages)
            _embed_svg_description(cur_file, c.descriptions[total_pages])
        end
        if total_pages == 1 && !occursin("%", c.path)
            orig_page1 = _svg_page_filename(c.path, 1)
            if isfile(orig_page1) && orig_page1 != c.path
                mv(orig_page1, c.path, force=true)
            end
            if haskey(c.descriptions, 1)
                _embed_svg_description(c.path, c.descriptions[1])
            end
        end
    elseif c.fmt == :html
        if !c.needs_new_surface
            Cairo.finish(c.surface)
        end
        total_pages = c.needs_new_surface ? c.page - 1 : c.page
        for (p, desc) in c.descriptions
            f = joinpath(c.htmldir, @sprintf("page_%03d.svg", p))
            _embed_svg_description(f, desc)
        end
        _write_html_presentation(c)
    elseif c.fmt == :png
        if !c.needs_new_surface
            cur_file = _png_page_filename(c.path, c.page)
            Cairo.write_to_png(c.surface, cur_file)
            Cairo.finish(c.surface)
        end
        total_pages = c.needs_new_surface ? c.page - 1 : c.page
        if total_pages == 1 && !occursin("%", c.path)
            orig_page1 = _png_page_filename(c.path, 1)
            if isfile(orig_page1) && orig_page1 != c.path
                mv(orig_page1, c.path, force=true)
            end
        end
    elseif c.fmt == :gif
        if !c.needs_new_surface
            frame_file = joinpath(c.tempdir, @sprintf("frame_%05d.png", c.page))
            Cairo.write_to_png(c.surface, frame_file)
            Cairo.finish(c.surface)
        end
        try
            if Sys.which("ffmpeg") !== nothing
                cmd = `ffmpeg -y -loglevel error -framerate $(c.fps) -i $(joinpath(c.tempdir, "frame_%05d.png")) -filter_complex "split[a][b];[a]palettegen[p];[b][p]paletteuse" $(c.path)`
                run(cmd)
            else
                @warn "ffmpeg not found in PATH; GIF was not generated at $(c.path)"
            end
        finally
            rm(c.tempdir, recursive=true, force=true)
        end
    end
    return
end

function BasicCompGeometry.open_canvas(f::Function, path::String, cw::Real, ch::Real; kwargs...)
    c = BasicCompGeometry.Canvas(path, cw, ch; kwargs...)
    try
        f(c)
    finally
        Cairo.finish(c)
    end
    return c
end

function Base.getproperty(c::Canvas, s::Symbol)
    if s === :cr || s === :surface
        _ensure_surface!(c)
    end
    return getfield(c, s)
end

Base.cconvert(::Type{Ptr{Cvoid}}, c::Canvas) = (_ensure_surface!(c); c.cr.ptr)
Base.unsafe_convert(::Type{Ptr{Cvoid}}, c::Canvas) = (_ensure_surface!(c); c.cr.ptr)

BasicCompGeometry.cairo_draw_setup(c::Canvas, bb::BBox{2,T}, cw::Real, ch::Real, margin::Real=20) where {T} =
    (_ensure_surface!(c); BasicCompGeometry.cairo_draw_setup(c.cr, bb, cw, ch, margin))

BasicCompGeometry.cairo_draw_points(c::Canvas, points, radius::Real=2) =
    (_ensure_surface!(c); BasicCompGeometry.cairo_draw_points(c.cr, points, radius))

BasicCompGeometry.cairo_draw_polygon(c::Canvas, poly; line_width=nothing, close::Bool=true) =
    (_ensure_surface!(c); BasicCompGeometry.cairo_draw_polygon(c.cr, poly; line_width=line_width, close=close))

# Cairo drawing primitives forwarding
for fn in (:save, :restore, :reset_transform, :new_path, :close_path, :stroke, :fill, :fill_preserve, :paint)
    @eval Cairo.$fn(c::Canvas) = (_ensure_surface!(c); Cairo.$fn(c.cr))
end

raw"""
    cairo_set_line_width(cr_or_canvas, a::Real)

Set the stroke line width on a Cairo context or `Canvas`, with automatic scaling and raster protection.

# Background & Problem Description
In computational geometry visualization, drawings are typically defined in mathematical
world coordinates (e.g., within a unit square [0, 1]^2 or arbitrary bounding box).
The canvas transformation is configured via `cairo_draw_setup(canvas, bb, cw, ch, margin)`,
which scales the world coordinates by a factor:
    `scale = min((cw - 2margin) / bb_width, (ch - 2margin) / bb_height)`

When setting line widths, code often passes values in one of two conventions:
1. **Mathematical / Bounding-box units**: e.g., `0.008` or `0.015` relative to a [0, 1] domain.
2. **Device / Pixel units**: e.g., `1.5`, `2.0`, or `3.0` pixels.

**The Discrepancy between Vector and Raster Surfaces:**
- **Vector Surfaces (`.pdf`, `.svg`)**: Vector viewers render sub-unit paths with vector precision,
  producing visible hairlines even for very small user-unit numbers.
- **Raster Surfaces (`.png`, `.gif` via `CairoImageSurface`)**: Cairo's underlying rasterizer (Pixman)
  interprets line widths directly in device pixels during scanline conversion. If a user passes a small
  mathematical fraction (e.g. `0.008`) without scaling it by `scale`, the device stroke width is `< 0.5 px`,
  which rounds down to **0 pixels**. As a result, polygon boundaries, grid lines, and curve strokes
  silently vanish and become completely invisible on PNGs and animated GIFs.

# Automatic Fix Implementation
To ensure consistent rendering across all supported formats without requiring manual conversions:
1. When targeting raster surfaces (`:png`, `:gif`) or when a fractional width `a < 1.0` is passed:
   - The active scaling factor `scale = sqrt(|det(CTM)|)` is extracted from the current transformation matrix.
   - The mathematical width is converted to device pixels: `w_dev = a * scale`.
   - The device width is rounded up to the nearest integer and clamped to at least `1.0 px`:
     `effective_w = max(1.0, ceil(a * scale))`.
2. When an explicit pixel width `a >= 1.0` is passed, the value is passed directly to Cairo as device pixels.
"""
function BasicCompGeometry.cairo_set_line_width(cr_or_canvas, a::Real)
    cr = cr_or_canvas isa Canvas ? (_ensure_surface!(cr_or_canvas); cr_or_canvas.cr) : cr_or_canvas
    if a < 1.0
        m = Cairo.get_matrix(cr)
        scale = sqrt(abs(m.xx * m.yy - m.xy * m.yx))
        if scale > 0
            effective_w = max(1.0, ceil(Float64(a) * scale))
            Cairo.set_line_width(cr, effective_w)
            return
        end
    end
    Cairo.set_line_width(cr, Float64(a))
end

Cairo.set_line_width(c::Canvas, a::Real) = BasicCompGeometry.cairo_set_line_width(c, a)

for fn in (:move_to, :line_to, :translate, :scale)
    @eval Cairo.$fn(c::Canvas, a, b) = (_ensure_surface!(c); Cairo.$fn(c.cr, a, b))
end

for fn in (:set_source_rgb, :rotate)
    @eval Cairo.$fn(c::Canvas, a, b, d...) = (_ensure_surface!(c); Cairo.$fn(c.cr, a, b, d...))
end

for fn in (:set_source_rgba, :rectangle)
    @eval Cairo.$fn(c::Canvas, a, b, d, e) = (_ensure_surface!(c); Cairo.$fn(c.cr, a, b, d, e))
end

Cairo.arc(c::Canvas, x, y, r, a1, a2) = (_ensure_surface!(c); Cairo.arc(c.cr, x, y, r, a1, a2))
Cairo.set_source_surface(c::Canvas, s, x, y) = (_ensure_surface!(c); Cairo.set_source_surface(c.cr, s, x, y))

# Cairo text primitives forwarding
Cairo.select_font_face(c::Canvas, family::AbstractString, slant, weight) =
    (_ensure_surface!(c); Cairo.select_font_face(c.cr, family, slant, weight))
Cairo.set_font_size(c::Canvas, size) = (_ensure_surface!(c); Cairo.set_font_size(c.cr, size))
Cairo.show_text(c::Canvas, text::AbstractString) = (_ensure_surface!(c); Cairo.show_text(c.cr, text))
Cairo.text_extents(c::Canvas, text::AbstractString) = (_ensure_surface!(c); Cairo.text_extents(c.cr, text))

end # module