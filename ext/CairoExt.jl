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
    for p in points
        C.new_path(cr)
        C.arc(cr, p[1], p[2], radius, 0.0, 2pi)
        C.fill(cr)
    end
    return
end

"""
    cairo_draw_polygon(cr, poly, close=true)

Draw the edges of a polygon (or point sequence) `poly` using the Cairo context `cr`.
If `close` is true, the polygon is closed by connecting the last vertex back to the first.
"""
function BasicCompGeometry.cairo_draw_polygon(cr, poly, close::Bool=true)
    C = Cairo
    n = length(poly)
    n == 0 && return

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
    Canvas(path::String, cw::Real, ch::Real; fps::Int=20)

Create a unified 2D vector / raster canvas supporting:
- `.pdf`: Multi-page PDF document
- `.svg`: Multi-page SVG sequence (e.g. `slide_001.svg`, `slide_002.svg`, ...)
- `.png`: Multi-frame PNG sequence (e.g. `frame_0001.png`, `frame_0002.png`, ...)
- `.gif`: High-quality animated GIF (via `ffmpeg` on finish)

Advance pages or animation frames with `Cairo.show_page(canvas)`.
Finalize and flush to disk with `Cairo.finish(canvas)`.
"""
mutable struct Canvas
    path::String
    cw::Float64
    ch::Float64
    fmt::Symbol             # :pdf, :svg, :png, :gif
    page::Int
    fps::Int
    tempdir::String
    surface::Cairo.CairoSurfaceBase
    cr::Cairo.CairoContext
end

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

function BasicCompGeometry.Canvas(path::String, cw::Real, ch::Real; fps::Int=20)
    ext = lowercase(splitext(path)[2])
    w = Float64(cw)
    h = Float64(ch)

    if ext == ".pdf"
        surf = Cairo.CairoPDFSurface(path, w, h)
        cr = Cairo.CairoContext(surf)
        return Canvas(path, w, h, :pdf, 1, fps, "", surf, cr)
    elseif ext == ".svg"
        fmt_path = _svg_page_filename(path, 1)
        surf = Cairo.CairoSVGSurface(fmt_path, w, h)
        cr = Cairo.CairoContext(surf)
        return Canvas(path, w, h, :svg, 1, fps, "", surf, cr)
    elseif ext == ".png"
        surf = Cairo.CairoImageSurface(Int(round(w)), Int(round(h)), Cairo.FORMAT_ARGB32)
        cr = Cairo.CairoContext(surf)
        return Canvas(path, w, h, :png, 1, fps, "", surf, cr)
    elseif ext == ".gif"
        tmp = mktempdir()
        surf = Cairo.CairoImageSurface(Int(round(w)), Int(round(h)), Cairo.FORMAT_ARGB32)
        cr = Cairo.CairoContext(surf)
        return Canvas(path, w, h, :gif, 1, fps, tmp, surf, cr)
    else
        error("Unsupported Canvas format: '$ext'. Expected .pdf, .svg, .png, or .gif")
    end
end

function Cairo.show_page(c::Canvas)
    if c.fmt == :pdf
        Cairo.show_page(c.cr)
        c.page += 1
    elseif c.fmt == :svg
        Cairo.finish(c.surface)
        c.page += 1
        next_file = _svg_page_filename(c.path, c.page)
        c.surface = Cairo.CairoSVGSurface(next_file, c.cw, c.ch)
        c.cr = Cairo.CairoContext(c.surface)
    elseif c.fmt == :png
        cur_file = _png_page_filename(c.path, c.page)
        Cairo.write_to_png(c.surface, cur_file)
        Cairo.finish(c.surface)
        c.page += 1
        c.surface = Cairo.CairoImageSurface(Int(round(c.cw)), Int(round(c.ch)), Cairo.FORMAT_ARGB32)
        c.cr = Cairo.CairoContext(c.surface)
    elseif c.fmt == :gif
        frame_file = joinpath(c.tempdir, @sprintf("frame_%05d.png", c.page))
        Cairo.write_to_png(c.surface, frame_file)
        Cairo.finish(c.surface)
        c.page += 1
        c.surface = Cairo.CairoImageSurface(Int(round(c.cw)), Int(round(c.ch)), Cairo.FORMAT_ARGB32)
        c.cr = Cairo.CairoContext(c.surface)
    end
    return c
end

function Cairo.finish(c::Canvas)
    if c.fmt == :pdf
        Cairo.finish(c.surface)
    elseif c.fmt == :svg
        Cairo.finish(c.surface)
        if c.page == 1 && !occursin("%", c.path)
            orig_page1 = _svg_page_filename(c.path, 1)
            if isfile(orig_page1) && orig_page1 != c.path
                mv(orig_page1, c.path, force=true)
            end
        end
    elseif c.fmt == :png
        cur_file = _png_page_filename(c.path, c.page)
        Cairo.write_to_png(c.surface, cur_file)
        Cairo.finish(c.surface)
        if c.page == 1 && !occursin("%", c.path)
            orig_page1 = _png_page_filename(c.path, 1)
            if isfile(orig_page1) && orig_page1 != c.path
                mv(orig_page1, c.path, force=true)
            end
        end
    elseif c.fmt == :gif
        frame_file = joinpath(c.tempdir, @sprintf("frame_%05d.png", c.page))
        Cairo.write_to_png(c.surface, frame_file)
        Cairo.finish(c.surface)
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

Base.cconvert(::Type{Ptr{Cvoid}}, c::Canvas) = c.cr.ptr
Base.unsafe_convert(::Type{Ptr{Cvoid}}, c::Canvas) = c.cr.ptr

BasicCompGeometry.cairo_draw_setup(c::Canvas, bb::BBox{2,T}, cw::Real, ch::Real, margin::Real=20) where {T} =
    BasicCompGeometry.cairo_draw_setup(c.cr, bb, cw, ch, margin)

BasicCompGeometry.cairo_draw_points(c::Canvas, points, radius::Real=2) =
    BasicCompGeometry.cairo_draw_points(c.cr, points, radius)

BasicCompGeometry.cairo_draw_polygon(c::Canvas, poly, close::Bool=true) =
    BasicCompGeometry.cairo_draw_polygon(c.cr, poly, close)

# Cairo drawing primitives forwarding
for fn in (:save, :restore, :new_path, :close_path, :stroke, :fill, :fill_preserve)
    @eval Cairo.$fn(c::Canvas) = Cairo.$fn(c.cr)
end

for fn in (:set_line_width, :paint)
    @eval Cairo.$fn(c::Canvas, a) = Cairo.$fn(c.cr, a)
end

for fn in (:move_to, :line_to, :translate, :scale)
    @eval Cairo.$fn(c::Canvas, a, b) = Cairo.$fn(c.cr, a, b)
end

for fn in (:set_source_rgb, :rotate)
    @eval Cairo.$fn(c::Canvas, a, b, d...) = Cairo.$fn(c.cr, a, b, d...)
end

for fn in (:set_source_rgba, :rectangle)
    @eval Cairo.$fn(c::Canvas, a, b, d, e) = Cairo.$fn(c.cr, a, b, d, e)
end

Cairo.arc(c::Canvas, x, y, r, a1, a2) = Cairo.arc(c.cr, x, y, r, a1, a2)
Cairo.set_source_surface(c::Canvas, s, x, y) = Cairo.set_source_surface(c.cr, s, x, y)

# Cairo text primitives forwarding
Cairo.select_font_face(c::Canvas, family::AbstractString, slant, weight) =
    Cairo.select_font_face(c.cr, family, slant, weight)
Cairo.set_font_size(c::Canvas, size) = Cairo.set_font_size(c.cr, size)
Cairo.show_text(c::Canvas, text::AbstractString) = Cairo.show_text(c.cr, text)
Cairo.text_extents(c::Canvas, text::AbstractString) = Cairo.text_extents(c.cr, text)

end # module