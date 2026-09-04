"""
    IpeDraw

A clean, robust Julia interface for programmatically generating Ipe 7 vector figures.
Integrates natively with `BasicCompGeometry` geometric primitives (Point, Segment, BBox, Polygon)
and provides high-level helpers for algorithmic and conceptual diagrams.
"""
module IpeDraw

using ..BasicCompGeometry
using Printf

export IpeCanvas, open_ipe, edit_ipe, add_preamble!
export draw_point!, draw_points!, draw_segment!, draw_box!, draw_polygon!
export draw_circle!, draw_arc!, draw_ellipse!, draw_elliptic_arc!
export draw_bezier!, draw_spline!, draw_bspline!, draw_polygon_with_holes!, ipe_group
export draw_bar!, draw_span!, draw_dimension!, draw_arrow!, draw_curved_arrow!
export draw_label!, set_layer!, add_layer!, add_view!, setup_transform!
export save_ipe, compile_pdf, save_figure_tex, export_figure

const DEFAULT_STYLE_FILE = normpath(joinpath(@__DIR__, "..", "assets", "default.ipe"))

"""
    IpeCanvas

A canvas holding geometry, layers, stylesheets, and metadata to be emitted as an Ipe 7 XML document.
"""
mutable struct IpeCanvas
    width::Float64
    height::Float64
    paper::String
    bbox::String
    preamble::String
    stylesheet_template::Union{String, Nothing}
    active_layer::String
    layers::Vector{String}
    views::Vector{String}
    elements::Vector{String}
    transform_fn::Union{Nothing, Function}

    function IpeCanvas(;
        width::Real = 576.0,
        height::Real = 504.0,
        paper::String = "576 504",
        bbox::String = "cropbox",
        preamble::String = "\\usepackage{amsmath,amssymb}\\def\\ipeMode{TRUE}\\def\\Sample{\\mathsf{R}}",
        template::Union{String, Nothing} = nothing,
        layer::String = "alpha"
    )
        new(
            Float64(width),
            Float64(height),
            paper,
            bbox,
            preamble,
            template,
            layer,
            [layer],
            String[],
            String[],
            nothing
        )
    end
end

# -----------------------------------------------------------------------------
# Attribute Formatting Helpers
# -----------------------------------------------------------------------------

function _fmt_attr(name::String, val::Union{Symbol, String, Nothing})
    val === nothing && return ""
    s = string(val)
    s == "none" && return ""
    if name == "arrow" || name == "rarrow"
        s = s == "pointed" ? "pointed/normal" : (s == "normal" ? "normal/normal" : s)
    end
    return " $name=\"$s\""
end

function _pt_str(p::Point{2})
    @sprintf("%.3f %.3f", p[1], p[2])
end

function _apply_tf(canvas::IpeCanvas, p::Point{2})
    canvas.transform_fn === nothing ? p : canvas.transform_fn(p)
end

_apply_tf(canvas::IpeCanvas, x::Real, y::Real) = _apply_tf(canvas, Point(Float64(x), Float64(y)))

# -----------------------------------------------------------------------------
# Layer Management
# -----------------------------------------------------------------------------

"""
    add_layer!(canvas, name)

Register a new layer on the canvas.
"""
function add_layer!(canvas::IpeCanvas, name::String)
    if !(name in canvas.layers)
        push!(canvas.layers, name)
    end
    return canvas
end

"""
    set_layer!(canvas, name)

Set the active drawing layer, automatically adding it if not present.
"""
function set_layer!(canvas::IpeCanvas, name::String)
    add_layer!(canvas, name)
    canvas.active_layer = name
    return canvas
end

"""
    add_view!(canvas, layers; active=nothing)

Add a view containing the specified list of layers.
"""
function add_view!(canvas::IpeCanvas, layers::Vector{String}; active=nothing)
    act = active === nothing ? layers[1] : string(active)
    push!(canvas.views, "<view layers=\"$(join(layers, " "))\" active=\"$act\"/>")
    return canvas
end

"""
    add_preamble!(canvas, snippet)

Append LaTeX packages or macros to the canvas LaTeX preamble.
"""
function add_preamble!(canvas::IpeCanvas, snippet::AbstractString)
    canvas.preamble = isempty(canvas.preamble) ? String(snippet) : canvas.preamble * "\n" * String(snippet)
    return canvas
end


# -----------------------------------------------------------------------------
# World-to-Canvas Affine Transformation
# -----------------------------------------------------------------------------

"""
    setup_transform!(canvas, world_bb::BBox{2}; margin=30.0, flip_y=false)

Automatically map a world bounding box onto the canvas area with given padding margin.
"""
function setup_transform!(
    canvas::IpeCanvas, world_bb::BBox{2, T};
    margin::Real = 30.0,
    flip_y::Bool = false
) where {T}
    bl = bottom_left(world_bb)
    tr = top_right(world_bb)
    w_w = tr[1] - bl[1]
    w_h = tr[2] - bl[2]
    w_w <= 0 && (w_w = 1.0)
    w_h <= 0 && (w_h = 1.0)

    avail_w = canvas.width - 2 * margin
    avail_h = canvas.height - 2 * margin
    scale = min(avail_w / w_w, avail_h / w_h)

    canvas.transform_fn = function(p::Point{2})
        x_c = margin + (p[1] - bl[1]) * scale
        y_c = flip_y ? (canvas.height - margin - (p[2] - bl[2]) * scale) : (margin + (p[2] - bl[2]) * scale)
        return Point(x_c, y_c)
    end
    return canvas
end

# -----------------------------------------------------------------------------
# Basic Geometry Drawing Primitives
# -----------------------------------------------------------------------------

"""
    draw_point!(canvas, p; stroke=:black, fill=:black, size=:normal, shape=:disk)

Draw a point as an Ipe symbol mark (e.g. `mark/disk(sx)`).
"""
function draw_point!(
    canvas::IpeCanvas, p::Point{2};
    stroke::Symbol = :black,
    fill::Symbol = :black,
    size::Symbol = :normal,
    shape::Symbol = :disk
)
    pt = _apply_tf(canvas, p)
    mark_name = shape == :disk ? "mark/disk(sx)" : (shape == :circle ? "mark/circle(sx)" : "mark/box(sx)")
    xml = "<use layer=\"$(canvas.active_layer)\" name=\"$mark_name\" pos=\"$(_pt_str(pt))\" size=\"$size\"$(_fmt_attr("stroke", stroke))$(_fmt_attr("fill", fill))/>"
    push!(canvas.elements, xml)
    return canvas
end

draw_point!(canvas::IpeCanvas, x::Real, y::Real; kwargs...) = draw_point!(canvas, Point(Float64(x), Float64(y)); kwargs...)

"""
    draw_points!(canvas, pts; stroke=:black, fill=:black, size=:normal, shape=:disk)

Draw multiple points.
"""
function draw_points!(canvas::IpeCanvas, pts; kwargs...)
    for p in pts
        draw_point!(canvas, p; kwargs...)
    end
    return canvas
end

"""
    draw_segment!(canvas, p1, p2; stroke=:black, pen=:normal, dash=:solid, arrow=:none, rarrow=:none)

Draw a straight line segment between two points.
"""
function draw_segment!(
    canvas::IpeCanvas, p1::Point{2}, p2::Point{2};
    stroke::Symbol = :black,
    pen::Symbol = :normal,
    dash::Symbol = :solid,
    arrow::Symbol = :none,
    rarrow::Symbol = :none
)
    q1 = _apply_tf(canvas, p1)
    q2 = _apply_tf(canvas, p2)
    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("pen", pen))$(_fmt_attr("dash", dash))$(_fmt_attr("arrow", arrow))$(_fmt_attr("rarrow", rarrow))"
    xml = "<path layer=\"$(canvas.active_layer)\"$attrs>\n$(_pt_str(q1)) m\n$(_pt_str(q2)) l\n</path>"
    push!(canvas.elements, xml)
    return canvas
end

draw_segment!(canvas::IpeCanvas, s::Segment{2}; kwargs...) = draw_segment!(canvas, s.p, s.q; kwargs...)
draw_segment!(canvas::IpeCanvas, x1::Real, y1::Real, x2::Real, y2::Real; kwargs...) = 
    draw_segment!(canvas, Point(Float64(x1), Float64(y1)), Point(Float64(x2), Float64(y2)); kwargs...)

"""
    draw_box!(canvas, bb; stroke=:black, fill=:none, pen=:normal, dash=:solid)

Draw a rectangular box (axis-aligned bounding box).
"""
function draw_box!(
    canvas::IpeCanvas, x1::Real, y1::Real, x2::Real, y2::Real;
    stroke::Symbol = :black,
    fill::Symbol = :none,
    pen::Symbol = :normal,
    dash::Symbol = :solid
)
    q1 = _apply_tf(canvas, x1, y1)
    q2 = _apply_tf(canvas, x2, y2)
    min_x, max_x = min(q1[1], q2[1]), max(q1[1], q2[1])
    min_y, max_y = min(q1[2], q2[2]), max(q1[2], q2[2])
    
    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("fill", fill))$(_fmt_attr("pen", pen))$(_fmt_attr("dash", dash))"
    xml = "<path layer=\"$(canvas.active_layer)\"$attrs>\n" *
          @sprintf("%.3f %.3f m\n%.3f %.3f l\n%.3f %.3f l\n%.3f %.3f l\nh\n</path>",
                  min_x, min_y, max_x, min_y, max_x, max_y, min_x, max_y)
    push!(canvas.elements, xml)
    return canvas
end

function draw_box!(canvas::IpeCanvas, bb::BBox{2}; kwargs...)
    bl = bottom_left(bb)
    tr = top_right(bb)
    draw_box!(canvas, bl[1], bl[2], tr[1], tr[2]; kwargs...)
end

"""
    draw_polygon!(canvas, poly; close=true, stroke=:black, fill=:none, pen=:normal, dash=:solid)

Draw a polygonal chain or closed polygon from an iterable of points.
"""
function draw_polygon!(
    canvas::IpeCanvas, poly;
    close::Bool = true,
    stroke::Symbol = :black,
    fill::Symbol = :none,
    pen::Symbol = :normal,
    dash::Symbol = :solid
)
    pts = [_apply_tf(canvas, p) for p in poly]
    length(pts) < 2 && return canvas

    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("fill", fill))$(_fmt_attr("pen", pen))$(_fmt_attr("dash", dash))"
    buf = IOBuffer()
    println(buf, "<path layer=\"$(canvas.active_layer)\"$attrs>")
    println(buf, "$(_pt_str(pts[1])) m")
    for i in 2:length(pts)
        println(buf, "$(_pt_str(pts[i])) l")
    end
    if close
        println(buf, "h")
    end
    print(buf, "</path>")
    push!(canvas.elements, String(take!(buf)))
    return canvas
end

"""
    draw_circle!(canvas, center, radius; stroke=:black, fill=:none, pen=:normal)

Draw a circle or disk.
"""
function draw_circle!(
    canvas::IpeCanvas, center::Point{2}, radius::Real;
    stroke::Symbol = :black,
    fill::Symbol = :none,
    pen::Symbol = :normal
)
    c = _apply_tf(canvas, center)
    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("fill", fill))$(_fmt_attr("pen", pen))"
    xml = "<path layer=\"$(canvas.active_layer)\"$attrs>\n" *
          @sprintf("%.3f 0 0 %.3f %.3f %.3f e\n</path>", radius, radius, c[1], c[2])
    push!(canvas.elements, xml)
    return canvas
end

draw_circle!(canvas::IpeCanvas, c::Sphere{2}; kwargs...) = draw_circle!(canvas, c.center, c.radius; kwargs...)
draw_circle!(canvas::IpeCanvas, x::Real, y::Real, radius::Real; kwargs...) = 
    draw_circle!(canvas, Point(Float64(x), Float64(y)), radius; kwargs...)

"""
    draw_arc!(canvas, center, radius, a1, a2; stroke=:black, pen=:normal, arrow=:none, rarrow=:none)
    draw_arc!(canvas, arc::CircleArc; ...)

Draw a circular arc from angle `a1` to `a2` (in radians).
"""
function draw_arc!(
    canvas::IpeCanvas, center::Point{2}, radius::Real, a1::Real, a2::Real;
    stroke::Symbol = :black,
    pen::Symbol = :normal,
    arrow::Symbol = :none,
    rarrow::Symbol = :none
)
    c = _apply_tf(canvas, center)
    p1 = Point(c[1] + radius * cos(a1), c[2] + radius * sin(a1))
    p2 = Point(c[1] + radius * cos(a2), c[2] + radius * sin(a2))
    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("pen", pen))$(_fmt_attr("arrow", arrow))$(_fmt_attr("rarrow", rarrow))"
    xml = "<path layer=\"$(canvas.active_layer)\"$attrs>\n" *
          @sprintf("%.3f 0 0 %.3f %.3f %.3f %.3f %.3f %.3f %.3f arc\n</path>",
                  radius, radius, c[1], c[2], p1[1], p1[2], p2[1], p2[2])
    push!(canvas.elements, xml)
    return canvas
end

draw_arc!(canvas::IpeCanvas, arc::CircleArc; kwargs...) = 
    draw_arc!(canvas, arc.center, arc.radius, arc.theta1, arc.theta2; kwargs...)

"""
    draw_ellipse!(canvas, center::Point{2}, r_major::Real, r_minor::Real; angle=0.0, ...)
    draw_ellipse!(canvas, ellipse::Ellipse; ...)

Draw an ellipse with given center, semi-axes `r_major` and `r_minor`, and orientation `angle` in radians.
"""
function draw_ellipse!(
    canvas::IpeCanvas, center::Point{2}, r_major::Real, r_minor::Real;
    angle::Real = 0.0,
    stroke::Symbol = :black,
    fill::Symbol = :none,
    pen::Symbol = :normal,
    dash::Symbol = :none,
    opacity::Union{Symbol, Nothing} = nothing,
    tiling::Symbol = :none
)
    c = _apply_tf(canvas, center)
    ca, sa = cos(Float64(angle)), sin(Float64(angle))
    a, b = Float64(r_major), Float64(r_minor)
    m11, m21 = a * ca, a * sa
    m12, m22 = -b * sa, b * ca
    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("fill", fill))$(_fmt_attr("pen", pen))$(_fmt_attr("dash", dash))$(_fmt_attr("opacity", opacity))$(_fmt_attr("tiling", tiling))"
    xml = "<path layer=\"$(canvas.active_layer)\"$attrs>\n" *
          @sprintf("%.4f %.4f %.4f %.4f %.3f %.3f e\n</path>", m11, m21, m12, m22, c[1], c[2])
    push!(canvas.elements, xml)
    return canvas
end

draw_ellipse!(canvas::IpeCanvas, e::Ellipse; kwargs...) =
    draw_ellipse!(canvas, e.center, e.r_major, e.r_minor; angle=e.angle, kwargs...)

"""
    draw_elliptic_arc!(canvas, center, r_major, r_minor, angle, a1, a2; ...)
    draw_arc!(canvas, arc::EllipticArc; ...)

Draw an elliptic arc from parametric angle `a1` to `a2` (in radians).
"""
function draw_elliptic_arc!(
    canvas::IpeCanvas, center::Point{2}, r_major::Real, r_minor::Real,
    angle::Real, a1::Real, a2::Real;
    stroke::Symbol = :black,
    pen::Symbol = :normal,
    arrow::Symbol = :none,
    rarrow::Symbol = :none
)
    c = _apply_tf(canvas, center)
    ca, sa = cos(Float64(angle)), sin(Float64(angle))
    a, b = Float64(r_major), Float64(r_minor)
    m11, m21 = a * ca, a * sa
    m12, m22 = -b * sa, b * ca
    p1 = Point(c[1] + m11 * cos(a1) + m12 * sin(a1), c[2] + m21 * cos(a1) + m22 * sin(a1))
    p2 = Point(c[1] + m11 * cos(a2) + m12 * sin(a2), c[2] + m21 * cos(a2) + m22 * sin(a2))
    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("pen", pen))$(_fmt_attr("arrow", arrow))$(_fmt_attr("rarrow", rarrow))"
    xml = "<path layer=\"$(canvas.active_layer)\"$attrs>\n" *
          @sprintf("%.4f %.4f %.4f %.4f %.3f %.3f %.3f %.3f %.3f %.3f arc\n</path>",
                  m11, m21, m12, m22, c[1], c[2], p1[1], p1[2], p2[1], p2[2])
    push!(canvas.elements, xml)
    return canvas
end

draw_arc!(canvas::IpeCanvas, arc::EllipticArc; kwargs...) =
    draw_elliptic_arc!(canvas, arc.ellipse.center, arc.ellipse.r_major, arc.ellipse.r_minor,
                       arc.ellipse.angle, arc.alpha1, arc.alpha2; kwargs...)

"""
    draw_bezier!(canvas, b::CubicBezier{2}; stroke=:black, pen=:normal, fill=:none, ...)

Draw a cubic Bézier curve.
"""
function draw_bezier!(
    canvas::IpeCanvas, b::CubicBezier{2};
    stroke::Symbol = :black,
    fill::Symbol = :none,
    pen::Symbol = :normal,
    dash::Symbol = :none,
    arrow::Symbol = :none,
    rarrow::Symbol = :none
)
    p0 = _apply_tf(canvas, b.p0)
    p1 = _apply_tf(canvas, b.p1)
    p2 = _apply_tf(canvas, b.p2)
    p3 = _apply_tf(canvas, b.p3)
    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("fill", fill))$(_fmt_attr("pen", pen))$(_fmt_attr("dash", dash))$(_fmt_attr("arrow", arrow))$(_fmt_attr("rarrow", rarrow))"
    xml = "<path layer=\"$(canvas.active_layer)\"$attrs>\n" *
          "$(_pt_str(p0)) m\n" *
          @sprintf("%.3f %.3f %.3f %.3f %.3f %.3f c\n</path>",
                  p1[1], p1[2], p2[1], p2[2], p3[1], p3[2])
    push!(canvas.elements, xml)
    return canvas
end

"""
    draw_spline!(canvas, spline::CubicSpline{2}; stroke=:black, pen=:normal, fill=:none, ...)
    draw_spline!(canvas, pts::AbsPntSeq{2}; method=:catmull_rom, closed=false, ...)

Draw a composite cubic spline curve.
"""
function draw_spline!(
    canvas::IpeCanvas, s::CubicSpline{2};
    stroke::Symbol = :black,
    fill::Symbol = :none,
    pen::Symbol = :normal,
    dash::Symbol = :none,
    arrow::Symbol = :none,
    rarrow::Symbol = :none
)
    isempty(s.segments) && return canvas
    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("fill", fill))$(_fmt_attr("pen", pen))$(_fmt_attr("dash", dash))$(_fmt_attr("arrow", arrow))$(_fmt_attr("rarrow", rarrow))"
    lines = String["<path layer=\"$(canvas.active_layer)\"$attrs>"]
    p0 = _apply_tf(canvas, s.segments[1].p0)
    push!(lines, "$(_pt_str(p0)) m")
    for seg in s.segments
        p1 = _apply_tf(canvas, seg.p1)
        p2 = _apply_tf(canvas, seg.p2)
        p3 = _apply_tf(canvas, seg.p3)
        push!(lines, @sprintf("%.3f %.3f %.3f %.3f %.3f %.3f c",
                              p1[1], p1[2], p2[1], p2[2], p3[1], p3[2]))
    end
    s.is_closed && push!(lines, "h")
    push!(lines, "</path>")
    push!(canvas.elements, join(lines, "\n"))
    return canvas
end

draw_spline!(canvas::IpeCanvas, pts::AbsPntSeq{2}; method::Symbol = :catmull_rom, closed::Bool = false, kwargs...) =
    draw_spline!(canvas, method == :natural ? interpolate_natural_spline(pts; closed=closed) :
                                             interpolate_catmull_rom(pts; closed=closed); kwargs...)

"""
    draw_bspline!(canvas, pts; closed=false, stroke=:black, pen=:normal, fill=:none, ...)

Draw an approximating cubic B-spline using native Ipe spline operators (`s` or `u`).
"""
function draw_bspline!(
    canvas::IpeCanvas, pts;
    closed::Bool = false,
    stroke::Symbol = :black,
    fill::Symbol = :none,
    pen::Symbol = :normal,
    dash::Symbol = :none
)
    n = length(pts)
    n >= 3 || error("B-spline requires at least 3 control points")
    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("fill", fill))$(_fmt_attr("pen", pen))$(_fmt_attr("dash", dash))"
    lines = String["<path layer=\"$(canvas.active_layer)\"$attrs>"]
    tf_pts = [_apply_tf(canvas, p) for p in pts]
    if closed
        pt_strs = [_pt_str(p) for p in tf_pts]
        push!(lines, "$(join(pt_strs, "\n")) u")
    else
        push!(lines, "$(_pt_str(tf_pts[1])) m")
        pt_strs = [_pt_str(p) for p in tf_pts[2:end]]
        push!(lines, "$(join(pt_strs, "\n")) s")
    end
    push!(lines, "</path>")
    push!(canvas.elements, join(lines, "\n"))
    return canvas
end

"""
    draw_polygon_with_holes!(canvas, outer, holes...; stroke=:black, fill=:gray7, ...)

Draw a filled planar region with interior holes using the even-odd winding rule (`fillrule="eofill"`).
"""
function draw_polygon_with_holes!(
    canvas::IpeCanvas, outer::AbsPntSeq{2}, holes::AbsPntSeq{2}...;
    stroke::Symbol = :black,
    fill::Symbol = :gray7,
    pen::Symbol = :normal,
    dash::Symbol = :none,
    opacity::Union{Symbol, Nothing} = nothing,
    tiling::Symbol = :none
)
    attrs = " fillrule=\"eofill\"$(_fmt_attr("stroke", stroke))$(_fmt_attr("fill", fill))$(_fmt_attr("pen", pen))$(_fmt_attr("dash", dash))$(_fmt_attr("opacity", opacity))$(_fmt_attr("tiling", tiling))"
    lines = String["<path layer=\"$(canvas.active_layer)\"$attrs>"]
    _append_subpath!(lines, canvas, outer)
    for hole in holes
        _append_subpath!(lines, canvas, hole)
    end
    push!(lines, "</path>")
    push!(canvas.elements, join(lines, "\n"))
    return canvas
end

function _append_subpath!(lines::Vector{String}, canvas::IpeCanvas, poly::AbsPntSeq{2})
    n = cardin(poly)
    n == 0 && return
    p0 = _apply_tf(canvas, poly[1])
    push!(lines, "$(_pt_str(p0)) m")
    for i in 2:n
        pi = _apply_tf(canvas, poly[i])
        push!(lines, "$(_pt_str(pi)) l")
    end
    push!(lines, "h")
end

"""
    ipe_group(f::Function, canvas::IpeCanvas; matrix=nothing, opacity=nothing)

Group elements emitted inside function `f(canvas)` under an Ipe `<group>` tag,
optionally applying an affine matrix transformation or group opacity.
"""
function ipe_group(
    f::Function, canvas::IpeCanvas;
    matrix::Union{Nothing, AbstractVector{<:Real}} = nothing,
    opacity::Union{Nothing, Symbol, String} = nothing
)
    attrs = ""
    if matrix !== nothing
        @assert length(matrix) == 6 "Matrix must have 6 values [m11, m21, m12, m22, dx, dy]"
        attrs *= @sprintf(" matrix=\"%.4f %.4f %.4f %.4f %.3f %.3f\"",
                          matrix[1], matrix[2], matrix[3], matrix[4], matrix[5], matrix[6])
    end
    if opacity !== nothing
        attrs *= " opacity=\"$(opacity)\""
    end
    push!(canvas.elements, "<group$attrs>")
    f(canvas)
    push!(canvas.elements, "</group>")
    return canvas
end


# -----------------------------------------------------------------------------
# High-Level Algorithmic & Conceptual Diagram Helpers
# -----------------------------------------------------------------------------

function _escape_ipe_xml(s::AbstractString)
    # Escape XML entities for Ipe text elements
    s = replace(s, "&" => "&amp;")
    s = replace(s, "<" => "&lt;")
    s = replace(s, ">" => "&gt;")
    return s
end

"""
    draw_label!(canvas, p, text; halign=:center, valign=:baseline, stroke=:black, size=:normal, style=:math)

Typeset a LaTeX text or math label with precise alignment. Accepts `LaTeXString` (`L"..."`), `String`, etc.
"""
function draw_label!(
    canvas::IpeCanvas, p::Point{2}, text;
    halign::Symbol = :center,
    valign::Symbol = :baseline,
    stroke::Symbol = :black,
    size::Symbol = :normal,
    style::Symbol = :math
)
    pt = _apply_tf(canvas, p)
    clean_text = string(text)
    # If math style and text not already wrapped in $, wrap it
    if style == :math && !startswith(strip(clean_text), "\$") && !startswith(strip(clean_text), "\\begin")
        clean_text = "\$" * clean_text * "\$"
    end
    xml_text = _escape_ipe_xml(clean_text)
    xml = "<text layer=\"$(canvas.active_layer)\" pos=\"$(_pt_str(pt))\" stroke=\"$stroke\" type=\"label\" size=\"$size\" halign=\"$halign\" valign=\"$valign\">$xml_text</text>"
    push!(canvas.elements, xml)
    return canvas
end

draw_label!(canvas::IpeCanvas, x::Real, y::Real, text; kwargs...) = draw_label!(canvas, Point(Float64(x), Float64(y)), text; kwargs...)

"""
    draw_bar!(canvas, x1, x2, y; height=18.0, stroke=:black, fill=:gray7, pen=:heavier, label_left=nothing, label_right=nothing)

Draw a horizontal array or domain bar representing a range [x1, x2], with optional endpoint ticks.
"""
function draw_bar!(
    canvas::IpeCanvas, x1::Real, x2::Real, y::Real;
    height::Real = 18.0,
    stroke::Symbol = :black,
    fill::Symbol = :gray7,
    pen::Symbol = :heavier,
    label_left = nothing,
    label_right = nothing
)
    y_bot = y - height / 2.0
    y_top = y + height / 2.0
    draw_box!(canvas, x1, y_bot, x2, y_top; stroke=stroke, fill=fill, pen=pen)

    if label_left !== nothing
        draw_label!(canvas, x1, y_bot - 12.0, string(label_left); size=:small, halign=:center, stroke=:gray2)
    end
    if label_right !== nothing
        draw_label!(canvas, x2, y_bot - 12.0, string(label_right); size=:small, halign=:center, stroke=:gray2)
    end
    return canvas
end

"""
    draw_span!(canvas, x1, x2, y; height=18.0, fill=:lightgreen, stroke=:darkgreen, pen=:heavier, dash=:solid)

Highlight an interval or active span within a bar.
"""
function draw_span!(
    canvas::IpeCanvas, x1::Real, x2::Real, y::Real;
    height::Real = 18.0,
    fill::Symbol = :lightgreen,
    stroke::Symbol = :darkgreen,
    pen::Symbol = :heavier,
    dash::Symbol = :solid
)
    y_bot = y - height / 2.0
    y_top = y + height / 2.0
    draw_box!(canvas, x1, y_bot, x2, y_top; fill=fill, stroke=stroke, pen=pen, dash=dash)
    return canvas
end

"""
    draw_dimension!(canvas, p1, p2; label=nothing, arrow=:both, stroke=:darkgreen, pen=:normal, label_offset=12.0)
    draw_dimension!(canvas, x1, x2, y; label=nothing, ...)

Draw a dimension line (double-arrow) with an optional text label.
"""
function draw_dimension!(
    canvas::IpeCanvas, p1::Point{2}, p2::Point{2};
    label = nothing,
    arrow::Symbol = :both,
    stroke::Symbol = :darkgreen,
    pen::Symbol = :normal,
    label_offset::Real = 12.0
)
    arr = arrow == :both ? :pointed : arrow
    rarr = arrow == :both ? :pointed : :none
    draw_segment!(canvas, p1, p2; stroke=stroke, pen=pen, arrow=arr, rarrow=rarr)

    if label !== nothing
        mid = Point((p1[1] + p2[1]) / 2.0, (p1[2] + p2[2]) / 2.0)
        dx = p2[1] - p1[1]
        dy = p2[2] - p1[2]
        len = hypot(dx, dy)
        normal = len > 0 ? Point(-dy / len, dx / len) : Point(0.0, 1.0)
        pos = Point(mid[1] + normal[1] * label_offset, mid[2] + normal[2] * label_offset)
        draw_label!(canvas, pos, label; stroke=stroke, size=:small, halign=:center)
    end
    return canvas
end

draw_dimension!(canvas::IpeCanvas, x1::Real, x2::Real, y::Real; label=nothing, label_offset::Real=-12.0, kwargs...) =
    draw_dimension!(canvas, Point(Float64(x1), Float64(y)), Point(Float64(x2), Float64(y));
                    label=label, label_offset=label_offset, kwargs...)

"""
    draw_arrow!(canvas, p1, p2; stroke=:darkgreen, pen=:fat, dash=:dashed, arrow=:pointed)

Draw a directed arrow (e.g. lifting arrow between domains).
"""
function draw_arrow!(
    canvas::IpeCanvas, p1::Point{2}, p2::Point{2};
    stroke::Symbol = :darkgreen,
    pen::Symbol = :fat,
    dash::Symbol = :dashed,
    arrow::Symbol = :pointed
)
    draw_segment!(canvas, p1, p2; stroke=stroke, pen=pen, dash=dash, arrow=arrow)
end

draw_arrow!(canvas::IpeCanvas, x1::Real, y1::Real, x2::Real, y2::Real; kwargs...) =
    draw_arrow!(canvas, Point(Float64(x1), Float64(y1)), Point(Float64(x2), Float64(y2)); kwargs...)

"""
    draw_curved_arrow!(canvas, p1, control_pt, p2; stroke=:blue, pen=:fat, arrow=:pointed)
    draw_curved_arrow!(canvas, p1, p2; bend=20.0, stroke=:blue, pen=:fat, arrow=:pointed)

Draw a curved Bézier arrow connecting two points.
"""
function draw_curved_arrow!(
    canvas::IpeCanvas, p1::Point{2}, control_pt::Point{2}, p2::Point{2};
    stroke::Symbol = :blue,
    pen::Symbol = :fat,
    arrow::Symbol = :pointed
)
    q1 = _apply_tf(canvas, p1)
    qc = _apply_tf(canvas, control_pt)
    q2 = _apply_tf(canvas, p2)
    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("pen", pen))$(_fmt_attr("arrow", arrow))"
    xml = "<path layer=\"$(canvas.active_layer)\"$attrs>\n" *
          "$(_pt_str(q1)) m\n$(_pt_str(qc))\n$(_pt_str(q2)) c\n</path>"
    push!(canvas.elements, xml)
    return canvas
end

function draw_curved_arrow!(
    canvas::IpeCanvas, p1::Point{2}, p2::Point{2};
    bend::Real = 20.0,
    stroke::Symbol = :blue,
    pen::Symbol = :fat,
    arrow::Symbol = :pointed
)
    mid = Point((p1[1] + p2[1]) / 2.0, (p1[2] + p2[2]) / 2.0)
    dx = p2[1] - p1[1]
    dy = p2[2] - p1[2]
    len = hypot(dx, dy)
    normal = len > 0 ? Point(-dy / len, dx / len) : Point(0.0, 1.0)
    ctrl = Point(mid[1] + normal[1] * bend, mid[2] + normal[2] * bend)
    return draw_curved_arrow!(canvas, p1, ctrl, p2; stroke=stroke, pen=pen, arrow=arrow)
end

# -----------------------------------------------------------------------------
# XML Generation & Compilation Pipelines
# -----------------------------------------------------------------------------

const DEFAULT_FALLBACK_IPESTYLE = raw"""<ipestyle name="basic">
<symbol name="mark/disk(sx)" transformations="translations">
<path fill="sym-stroke">
0.6 0 0 0.6 0 0 e
</path>
</symbol>
<symbol name="mark/circle(sx)" transformations="translations">
<path fill="sym-fill" stroke="sym-stroke">
0.6 0 0 0.6 0 0 e
</path>
</symbol>
<symbol name="mark/box(sx)" transformations="translations">
<path fill="sym-fill" stroke="sym-stroke">
-0.6 -0.6 m
0.6 -0.6 l
0.6 0.6 l
-0.6 0.6 l
h
</path>
</symbol>
<symbol name="arrow/pointed(spx)">
<path stroke="sym-stroke" fill="sym-stroke" pen="sym-pen">
0 0 m
-1 0.333 l
-0.8 0 l
-1 -0.333 l
h
</path>
</symbol>
<symbol name="arrow/normal(spx)">
<path stroke="sym-stroke" fill="sym-stroke" pen="sym-pen">
0 0 m
-1 0.333 l
-1 -0.333 l
h
</path>
</symbol>
<color name="red" value="1 0 0"/>
<color name="green" value="0 1 0"/>
<color name="blue" value="0 0 1"/>
<color name="yellow" value="1 1 0"/>
<color name="darkred" value="0.7 0 0"/>
<color name="darkgreen" value="0 0.5 0"/>
<color name="darkblue" value="0 0 0.6"/>
<color name="lightgreen" value="0.8 1 0.8"/>
<color name="lightblue" value="0.85 0.9 1"/>
<color name="lightred" value="1 0.85 0.85"/>
<color name="lightgray" value="0.9 0.9 0.9"/>
<color name="gray1" value="0.1 0.1 0.1"/>
<color name="gray2" value="0.2 0.2 0.2"/>
<color name="gray3" value="0.3 0.3 0.3"/>
<color name="gray4" value="0.4 0.4 0.4"/>
<color name="gray5" value="0.5 0.5 0.5"/>
<color name="gray6" value="0.6 0.6 0.6"/>
<color name="gray7" value="0.88 0.88 0.88"/>
<color name="gray8" value="0.8 0.8 0.8"/>
<color name="gray9" value="0.9 0.9 0.9"/>
<pen name="extra thin" value="0.05"/>
<pen name="thin" value="0.2"/>
<pen name="heavier" value="0.8"/>
<pen name="fat" value="1.2"/>
<pen name="ultrafat" value="2"/>
<pen name="ultrafat 4.0" value="4"/>
<pen name="ultrathin 0.5" value="0.5"/>
<pen name="ultrathin 0.25" value="0.25"/>
<opacity name="10%" value="0.1"/>
<opacity name="20%" value="0.2"/>
<opacity name="30%" value="0.3"/>
<opacity name="40%" value="0.4"/>
<opacity name="50%" value="0.5"/>
<opacity name="60%" value="0.6"/>
<opacity name="70%" value="0.7"/>
<opacity name="80%" value="0.8"/>
<opacity name="90%" value="0.9"/>
<textsize name="Huge" value="\Huge"/>
<textsize name="LARGE" value="\LARGE"/>
<textsize name="Large" value="\Large"/>
<textsize name="large" value="\large"/>
<textsize name="small" value="\small"/>
<textsize name="tiny" value="\tiny"/>
<textsize name="footnote" value="\footnotesize"/>
<tiling name="falling" angle="-60" step="4" width="1"/>
<tiling name="rising" angle="30" step="4" width="1"/>
<tiling name="hatch" angle="45" step="4" width="0.5"/>
<tiling name="crosshatch" angle="45" step="4" width="0.5"/>
<tiling name="vertical" angle="90" step="4" width="0.5"/>
<tiling name="horizontal" angle="0" step="4" width="0.5"/>
<layout paper="576 504" origin="0 0" frame="576 504" crop="yes"/>
</ipestyle>"""

"""
    to_xml(canvas)

Generate the complete, self-contained Ipe 7 XML document string.
"""
function to_xml(canvas::IpeCanvas)
    template = if canvas.stylesheet_template !== nothing && isfile(canvas.stylesheet_template)
        read(canvas.stylesheet_template, String)
    elseif isfile(DEFAULT_STYLE_FILE)
        read(DEFAULT_STYLE_FILE, String)
    else
        "<?xml version=\"1.0\"?>\n<!DOCTYPE ipe SYSTEM \"ipe.dtd\">\n<ipe version=\"70218\">\n<info tex=\"pdflatex\"/>\n" *
        DEFAULT_FALLBACK_IPESTYLE * "\n<page/>\n</ipe>"
    end

    # 1. Ensure bbox="cropbox" in info tag
    if occursin(r"<info\b", template)
        if !occursin("bbox=", template)
            template = replace(template, r"<info\b([^>]*)/>" => s"<info\1 bbox=\"cropbox\"/>")
        end
    end

    # 2. Ensure layout has crop="yes" and canvas paper sizing
    if occursin(r"<layout\b[^>]*/>", template)
        template = replace(template, r"<layout\b[^>]*/>" => "<layout paper=\"$(canvas.paper)\" origin=\"0 0\" frame=\"$(canvas.paper)\" crop=\"yes\"/>")
    elseif occursin(r"</ipestyle>\s*<page>", template)
        style_layout = "<ipestyle name=\"page_layout\">\n<layout paper=\"$(canvas.paper)\" origin=\"0 0\" frame=\"$(canvas.paper)\" crop=\"yes\"/>\n</ipestyle>"
        template = replace(template, r"</ipestyle>\s*<page>" => "</ipestyle>\n" * style_layout * "\n<page>")
    end

    # 3. Enhanced preamble
    preamble_xml = "<preamble>\n$(canvas.preamble)\n</preamble>"
    if occursin(r"<preamble>.*?</preamble>"s, template)
        template = replace(template, r"<preamble>.*?</preamble>"s => preamble_xml)
    end

    # 4. Build Page Content
    buf = IOBuffer()
    println(buf, "<page>")
    for layer in canvas.layers
        println(buf, "<layer name=\"$layer\"/>")
    end
    if isempty(canvas.views)
        println(buf, "<view layers=\"$(join(canvas.layers, " "))\" active=\"$(canvas.active_layer)\"/>")
    else
        for v in canvas.views
            println(buf, v)
        end
    end
    for elem in canvas.elements
        println(buf, elem)
    end
    println(buf, "</page>")
    page_xml = String(take!(buf))

    # Replace <page>...</page> block using a function to avoid backslash escape issues
    result = replace(template, r"<page>.*?</page>"s => (m -> page_xml))
    return result
end

"""
    save_ipe(canvas, filename)

Save the canvas content to an `.ipe` XML file.
"""
function save_ipe(canvas::IpeCanvas, filename::String)
    mkpath(dirname(abspath(filename)))
    xml_data = to_xml(canvas)
    write(filename, xml_data)
    return filename
end

"""
    compile_pdf(canvas_or_file, output_pdf=nothing)

Compile the canvas or an `.ipe` file to a cropped vector PDF using `ipetoipe -pdf`.
"""
function compile_pdf(canvas::IpeCanvas, output_pdf::String)
    if Sys.which("ipetoipe") === nothing
        @warn "ipetoipe command not found in PATH. Please install Ipe 7 to compile PDF figures."
        return false
    end
    mkpath(dirname(abspath(output_pdf)))
    ipe_temp = tempname() * ".ipe"
    save_ipe(canvas, ipe_temp)
    success = run(`ipetoipe -pdf $ipe_temp $output_pdf`).exitcode == 0
    rm(ipe_temp, force=true)
    return success
end

function compile_pdf(ipe_file::String, output_pdf::String=replace(ipe_file, r"\.ipe$" => ".pdf"))
    if Sys.which("ipetoipe") === nothing
        @warn "ipetoipe command not found in PATH. Please install Ipe 7 to compile PDF figures."
        return false
    end
    run(`ipetoipe -pdf $ipe_file $output_pdf`).exitcode == 0
end

"""
    save_figure_tex(filename; figure_name, caption="", label="")

Emit a standalone companion LaTeX figure fragment (`_fig.tex`) ready for `\\input`.
"""
function save_figure_tex(
    filename::String;
    figure_name::String = replace(basename(filename), r"(_fig)?\.tex$" => ""),
    caption::String = "",
    label::String = figure_name
)
    mkpath(dirname(abspath(filename)))
    tex = """
\\begin{figure}[t]
    \\centering
    \\IncludeGraphics{\\File{figs/$(figure_name)}}
    \\caption{$(caption)}
    \\figlab{$(label)}
\\end{figure}
"""
    write(filename, tex)
    return filename
end

"""
    export_figure(canvas, base_path; caption="", label="")

Export all three companion artifacts in one call:
`<base_path>.ipe`, `<base_path>.pdf`, and `<base_path>_fig.tex`.
"""
function export_figure(
    canvas::IpeCanvas, base_path::String;
    caption::String = "",
    label::String = replace(basename(base_path), "_" => ":")
)
    ipe_path = base_path * ".ipe"
    pdf_path = base_path * ".pdf"
    tex_path = base_path * "_fig.tex"

    save_ipe(canvas, ipe_path)
    compile_pdf(canvas, pdf_path)
    save_figure_tex(tex_path; figure_name=basename(base_path), caption=caption, label=label)

    return (ipe=ipe_path, pdf=pdf_path, tex=tex_path)
end

"""
    open_ipe(base_path; caption="", label="", kwargs...) do canvas
        ...
    end

Convenient block syntax to construct, save, compile, and generate LaTeX fragments in one go.
"""
function open_ipe(f::Function, base_path::String; caption::String="", label::String="", kwargs...)
    clean_base = replace(base_path, r"(\.ipe|\.pdf|_fig\.tex)$" => "")
    canvas = IpeCanvas(; kwargs...)
    f(canvas)
    lbl = isempty(label) ? replace(basename(clean_base), "_" => ":") : label
    return export_figure(canvas, clean_base; caption=caption, label=lbl)
end

"""
    edit_ipe(filename::String)

Launch the Ipe GUI editor on `filename` in the background.
"""
function edit_ipe(filename::String)
    if Sys.which("ipe") === nothing
        @warn "ipe executable not found in PATH. Please install Ipe 7."
        return nothing
    end
    run(`ipe $filename`, wait=false)
end

BasicCompGeometry.width(c::IpeCanvas) = c.width
BasicCompGeometry.height(c::IpeCanvas) = c.height

end # module IpeDraw
