"""
    IpeDraw

A clean, robust Julia interface for programmatically generating Ipe 7 vector figures.
Integrates natively with `BasicCompGeometry` geometric primitives (Point, Segment, BBox, Polygon)
and provides high-level helpers for algorithmic and conceptual diagrams.
"""
module IpeDraw

using ..BasicCompGeometry
using Printf

export IpeCanvas, Viewport, PageBox, Style, Theme, publication_theme
export open_ipe, figure, edit_ipe, add_preamble!
export draw_point!, draw_points!, draw_segment!, draw_box!, draw_polygon!
export draw_circle!, draw_arc!, draw_ellipse!, draw_elliptic_arc!
export draw_bezier!, draw_spline!, draw_bspline!, draw_polygon_with_holes!, ipe_group
export draw_bar!, draw_span!, draw_dimension!, draw_arrow!, draw_curved_arrow!
export draw_label!, set_layer!, add_layer!, add_view!, setup_transform!
export save_ipe, compile_pdf, save_figure_tex, export_figure
export draw!, mark!, label!, layer, with_style, fit!, page_space, inset, clip_to
export legend!, scale_bar!

const DEFAULT_STYLE_FILE = normpath(joinpath(@__DIR__, "..", "assets", "default.ipe"))

"""A uniform world-to-page transformation used by an `IpeCanvas`."""
struct Viewport
    scale::Float64
    tx::Float64
    ty::Float64
    flip_y::Bool
end

"""A rectangle in fixed Ipe page coordinates, used to place insets and overlays."""
struct PageBox
    x::Float64
    y::Float64
    width::Float64
    height::Float64

    function PageBox(x::Real, y::Real, width::Real, height::Real)
        width > 0 || throw(ArgumentError("page-box width must be positive"))
        height > 0 || throw(ArgumentError("page-box height must be positive"))
        new(Float64(x), Float64(y), Float64(width), Float64(height))
    end
end

"""Reusable keyword options for `draw!`, `mark!`, and `label!`."""
struct Style
    values::NamedTuple
end

Style(; kwargs...) = Style((; kwargs...))

"""A reusable collection of named `Style` objects."""
struct Theme
    styles::Dict{Symbol,Style}
end

function Theme(; kwargs...)
    styles = Dict{Symbol,Style}()
    for (name, value) in kwargs
        styles[name] = value isa Style ? value : Style(value)
    end
    return Theme(styles)
end

Base.getindex(theme::Theme, name::Symbol) = theme.styles[name]

function Base.getproperty(theme::Theme, name::Symbol)
    name === :styles && return getfield(theme, :styles)
    haskey(getfield(theme, :styles), name) && return getfield(theme, :styles)[name]
    return getfield(theme, name)
end

Base.propertynames(theme::Theme, private::Bool=false) =
    private ? (:styles, keys(theme.styles)...) : Tuple(keys(theme.styles))

"""Return restrained defaults suitable for papers and lecture notes."""
function publication_theme()
    return Theme(
        region=Style(fill_opacity=0.2, pen=:heavier),
        boundary=Style(stroke=:black, pen=:heavier),
        point=Style(stroke=:black, fill=:black, size=:normal),
        annotation=Style(stroke=:black, size=:small),
    )
end

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
    viewport::Union{Nothing, Viewport}
    active_style::Style

    function IpeCanvas(;
        width::Real = 576.0,
        height::Real = 504.0,
        paper::Union{String, Nothing} = nothing,
        bbox::String = "cropbox",
        preamble::String = "\\usepackage{amsmath,amssymb}\\def\\ipeMode{TRUE}\\def\\Sample{\\mathsf{R}}",
        template::Union{String, Nothing} = nothing,
        layer::String = "alpha"
    )
        paper_value = isnothing(paper) ? "$(Float64(width)) $(Float64(height))" : paper
        new(
            Float64(width),
            Float64(height),
            paper_value,
            bbox,
            preamble,
            template,
            layer,
            [layer],
            String[],
            String[],
            nothing,
            Style()
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
    vp = canvas.viewport
    vp === nothing && return p
    y_scale = vp.flip_y ? -vp.scale : vp.scale
    return Point(vp.tx + vp.scale * p.x, vp.ty + y_scale * p.y)
end

_apply_tf(canvas::IpeCanvas, x::Real, y::Real) = _apply_tf(canvas, Point(Float64(x), Float64(y)))
_apply_length(canvas::IpeCanvas, value::Real) =
    canvas.viewport === nothing ? Float64(value) : canvas.viewport.scale * Float64(value)

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
    return _setup_transform!(canvas, world_bb, PageBox(0, 0, canvas.width, canvas.height);
                             margin=margin, flip_y=flip_y)
end

function _setup_transform!(
    canvas::IpeCanvas, world_bb::BBox{2, T}, page_box::PageBox;
    margin::Real=0.0,
    flip_y::Bool=false,
) where {T}
    bl = bottom_left(world_bb)
    tr = top_right(world_bb)
    w_w = tr[1] - bl[1]
    w_h = tr[2] - bl[2]
    w_w <= 0 && (w_w = 1.0)
    w_h <= 0 && (w_h = 1.0)

    avail_w = page_box.width - 2 * margin
    avail_h = page_box.height - 2 * margin
    avail_w > 0 && avail_h > 0 || throw(ArgumentError("margin leaves no room in page box"))
    scale = min(avail_w / w_w, avail_h / w_h)

    x0 = page_box.x + margin + (avail_w - scale * w_w) / 2
    y0 = page_box.y + margin + (avail_h - scale * w_h) / 2
    tx = x0 - scale * bl.x
    ty = flip_y ? 2 * page_box.y + page_box.height - y0 + scale * bl.y : y0 - scale * bl.y
    canvas.viewport = Viewport(scale, tx, ty, flip_y)
    return canvas
end

"""Fit bounded geometry into `canvas` with a margin measured in page units."""
function fit!(canvas::IpeCanvas, objects; margin::Real = 30.0, flip_y::Bool = false)
    items = objects isa Tuple || objects isa AbstractVector ? objects : (objects,)
    setup_transform!(canvas, union_bbox(items...); margin=margin, flip_y=flip_y)
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
    dash::Symbol = :solid,
    opacity::Union{Symbol, String, Nothing} = nothing
)
    q1 = _apply_tf(canvas, x1, y1)
    q2 = _apply_tf(canvas, x2, y2)
    min_x, max_x = min(q1[1], q2[1]), max(q1[1], q2[1])
    min_y, max_y = min(q1[2], q2[2]), max(q1[2], q2[2])
    
    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("fill", fill))$(_fmt_attr("pen", pen))$(_fmt_attr("dash", dash))$(_fmt_attr("opacity", opacity))"
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
    dash::Symbol = :solid,
    opacity::Union{Symbol, String, Nothing} = nothing
)
    pts = [_apply_tf(canvas, p) for p in poly]
    length(pts) < 2 && return canvas

    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("fill", fill))$(_fmt_attr("pen", pen))$(_fmt_attr("dash", dash))$(_fmt_attr("opacity", opacity))"
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
    pen::Symbol = :normal,
    dash::Symbol = :solid,
    opacity::Union{Symbol, String, Nothing} = nothing
)
    c = _apply_tf(canvas, center)
    r = _apply_length(canvas, radius)
    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("fill", fill))$(_fmt_attr("pen", pen))$(_fmt_attr("dash", dash))$(_fmt_attr("opacity", opacity))"
    xml = "<path layer=\"$(canvas.active_layer)\"$attrs>\n" *
          @sprintf("%.3f 0 0 %.3f %.3f %.3f e\n</path>", r, r, c.x, c.y)
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
    r = _apply_length(canvas, radius)
    direction = canvas.viewport !== nothing && canvas.viewport.flip_y ? -1.0 : 1.0
    p1 = Point(c.x + r * cos(a1), c.y + direction * r * sin(a1))
    p2 = Point(c.x + r * cos(a2), c.y + direction * r * sin(a2))
    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("pen", pen))$(_fmt_attr("arrow", arrow))$(_fmt_attr("rarrow", rarrow))"
    xml = "<path layer=\"$(canvas.active_layer)\"$attrs>\n" *
          @sprintf("%.3f 0 0 %.3f %.3f %.3f %.3f %.3f %.3f %.3f arc\n</path>",
                  r, direction * r, c.x, c.y, p1.x, p1.y, p2.x, p2.y)
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
    opacity::Union{Symbol, String, Nothing} = nothing,
    tiling::Symbol = :none
)
    c = _apply_tf(canvas, center)
    ca, sa = cos(Float64(angle)), sin(Float64(angle))
    a, b = _apply_length(canvas, r_major), _apply_length(canvas, r_minor)
    direction = canvas.viewport !== nothing && canvas.viewport.flip_y ? -1.0 : 1.0
    m11, m21 = a * ca, a * sa
    m12, m22 = -b * sa, b * ca
    m21 *= direction
    m22 *= direction
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
    a, b = _apply_length(canvas, r_major), _apply_length(canvas, r_minor)
    direction = canvas.viewport !== nothing && canvas.viewport.flip_y ? -1.0 : 1.0
    m11, m21 = a * ca, a * sa
    m12, m22 = -b * sa, b * ca
    m21 *= direction
    m22 *= direction
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
    ipe_group(f::Function, canvas::IpeCanvas; matrix=nothing, opacity=nothing, clip=nothing)

Group elements emitted inside function `f(canvas)` under an Ipe `<group>` tag,
optionally applying an affine matrix, opacity, or Ipe path clip.
"""
function ipe_group(
    f::Function, canvas::IpeCanvas;
    matrix::Union{Nothing, AbstractVector{<:Real}} = nothing,
    opacity::Union{Nothing, Symbol, String} = nothing,
    clip::Union{Nothing,AbstractString} = nothing,
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
    clip === nothing || (attrs *= " clip=\"$(_escape_ipe_xml(clip))\"")
    push!(canvas.elements, "<group$attrs>")
    try
        f(canvas)
    finally
        push!(canvas.elements, "</group>")
    end
    return canvas
end

# -----------------------------------------------------------------------------
# Compact Geometry-Aware Drawing Interface
# -----------------------------------------------------------------------------

function _style_options(canvas::IpeCanvas, style::Union{Style,Nothing}, overrides::NamedTuple, allowed)
    unknown = filter(key -> !(key in allowed), keys(overrides))
    isempty(unknown) || throw(ArgumentError("unsupported drawing options: $(join(unknown, ", "))"))
    local_style = style === nothing ? NamedTuple() : style.values
    merged = merge(canvas.active_style.values, local_style, overrides)
    return (; (key => value for (key, value) in pairs(merged) if key in allowed)...)
end

function _opacity_name(value::Real)
    0 < value < 1 || throw(ArgumentError("named Ipe opacity must be strictly between zero and one"))
    percent = round(Int, 100 * value)
    percent % 10 == 0 ||
        throw(ArgumentError("Ipe opacity must be a multiple of 0.1"))
    return Symbol("$(percent)%")
end

_opacity_name(value::Union{Symbol,String}) = value

function _draw_fill_opacity!(renderer, canvas, object, options::NamedTuple)
    fill_opacity = get(options, :fill_opacity, nothing)
    base = (; (key => value for (key, value) in pairs(options) if key != :fill_opacity)...)
    fill_opacity === nothing && return renderer(canvas, object; base...)

    fill = get(base, :fill, :none)
    stroke = get(base, :stroke, :black)
    visible_fill = fill != :none && !(fill_opacity isa Real && iszero(fill_opacity))
    if visible_fill
        opacity = fill_opacity isa Real && isone(fill_opacity) ? nothing : _opacity_name(fill_opacity)
        fill_options = merge(base, (stroke=:none, opacity=opacity))
        renderer(canvas, object; fill_options...)
    end
    if stroke != :none
        stroke_options = merge(base, (fill=:none,))
        renderer(canvas, object; stroke_options...)
    end
    return canvas
end

const PATH_STYLE = (:stroke, :fill, :pen, :dash, :opacity, :fill_opacity)
const CURVE_STYLE = (:stroke, :pen, :dash, :arrow, :rarrow)

"""Draw a geometric object using multiple dispatch and optional reusable `Style`."""
function draw!(canvas::IpeCanvas, p::Point{2}; style::Union{Style,Nothing}=nothing, kwargs...)
    options = _style_options(canvas, style, (; kwargs...), (:stroke, :fill, :size, :shape))
    return draw_point!(canvas, p; options...)
end

function draw!(canvas::IpeCanvas, s::Segment{2}; style::Union{Style,Nothing}=nothing, kwargs...)
    options = _style_options(canvas, style, (; kwargs...), CURVE_STYLE)
    return draw_segment!(canvas, s; options...)
end

function draw!(canvas::IpeCanvas, bb::BBox{2}; style::Union{Style,Nothing}=nothing, kwargs...)
    options = _style_options(canvas, style, (; kwargs...), PATH_STYLE)
    return _draw_fill_opacity!(draw_box!, canvas, bb, options)
end

function draw!(canvas::IpeCanvas, circle::Sphere{2}; style::Union{Style,Nothing}=nothing, kwargs...)
    options = _style_options(canvas, style, (; kwargs...), PATH_STYLE)
    return _draw_fill_opacity!(draw_circle!, canvas, circle, options)
end

function draw!(canvas::IpeCanvas, arc::CircleArc; style::Union{Style,Nothing}=nothing, kwargs...)
    options = _style_options(canvas, style, (; kwargs...), CURVE_STYLE)
    return draw_arc!(canvas, arc; options...)
end

function draw!(canvas::IpeCanvas, ellipse::Ellipse; style::Union{Style,Nothing}=nothing, kwargs...)
    allowed = (PATH_STYLE..., :tiling)
    options = _style_options(canvas, style, (; kwargs...), allowed)
    return _draw_fill_opacity!(draw_ellipse!, canvas, ellipse, options)
end

function draw!(canvas::IpeCanvas, arc::EllipticArc; style::Union{Style,Nothing}=nothing, kwargs...)
    options = _style_options(canvas, style, (; kwargs...), CURVE_STYLE)
    return draw_arc!(canvas, arc; options...)
end

function draw!(canvas::IpeCanvas, curve::CubicBezier{2}; style::Union{Style,Nothing}=nothing, kwargs...)
    options = _style_options(canvas, style, (; kwargs...), (:stroke, :fill, :pen, :dash, :arrow, :rarrow))
    return draw_bezier!(canvas, curve; options...)
end

function draw!(canvas::IpeCanvas, curve::CubicSpline{2}; style::Union{Style,Nothing}=nothing, kwargs...)
    options = _style_options(canvas, style, (; kwargs...), (:stroke, :fill, :pen, :dash, :arrow, :rarrow))
    return draw_spline!(canvas, curve; options...)
end

function draw!(canvas::IpeCanvas, polygon::AbsPntSeq{2}; style::Union{Style,Nothing}=nothing, kwargs...)
    options = _style_options(canvas, style, (; kwargs...), (PATH_STYLE..., :close))
    return _draw_fill_opacity!(draw_polygon!, canvas, polygon, options)
end

function draw!(canvas::IpeCanvas, polygon::AbstractVector{<:Point{2}}; style::Union{Style,Nothing}=nothing, kwargs...)
    options = _style_options(canvas, style, (; kwargs...), (PATH_STYLE..., :close))
    return _draw_fill_opacity!(draw_polygon!, canvas, polygon, options)
end

function draw!(canvas::IpeCanvas, objects::AbstractVector; kwargs...)
    for object in objects
        draw!(canvas, object; kwargs...)
    end
    return canvas
end

"""Draw one point or a collection of points as fixed-size page marks."""
mark!(canvas::IpeCanvas, p::Point{2}; kwargs...) = draw!(canvas, p; kwargs...)

function mark!(canvas::IpeCanvas, points; style::Union{Style,Nothing}=nothing, kwargs...)
    options = _style_options(canvas, style, (; kwargs...), (:stroke, :fill, :size, :shape))
    return draw_points!(canvas, points; options...)
end

"""Run `f` on a layer and restore the previously active layer afterward."""
function layer(f::Function, canvas::IpeCanvas, name::Union{Symbol,AbstractString})
    previous = canvas.active_layer
    set_layer!(canvas, string(name))
    try
        f(canvas)
    finally
        set_layer!(canvas, previous)
    end
    return canvas
end

"""Apply style defaults within `f` and restore the previous defaults afterward."""
function with_style(f::Function, canvas::IpeCanvas, style::Style=Style(); kwargs...)
    previous = canvas.active_style
    canvas.active_style = Style(merge(previous.values, style.values, (; kwargs...)))
    try
        f(canvas)
    finally
        canvas.active_style = previous
    end
    return canvas
end

"""Draw in fixed page coordinates within `f`, restoring the world viewport afterward."""
function page_space(f::Function, canvas::IpeCanvas)
    previous = canvas.viewport
    canvas.viewport = nothing
    try
        f(canvas)
    finally
        canvas.viewport = previous
    end
    return canvas
end

function _page_clip_path(box::PageBox)
    x1, y1 = box.x, box.y
    x2, y2 = x1 + box.width, y1 + box.height
    return @sprintf("%.3f %.3f m %.3f %.3f l %.3f %.3f l %.3f %.3f l h",
                    x1, y1, x2, y1, x2, y2, x1, y2)
end

function _world_page_box(canvas::IpeCanvas, bb::BBox{2})
    p = _apply_tf(canvas, bottom_left(bb))
    q = _apply_tf(canvas, top_right(bb))
    return PageBox(min(p.x, q.x), min(p.y, q.y), abs(q.x - p.x), abs(q.y - p.y))
end

"""
    clip_to(f, canvas, box)

Clip everything emitted by `f` to a world-coordinate `BBox` or page-coordinate
`PageBox`. Ipe retains the clipped objects for later editing.
"""
function clip_to(f::Function, canvas::IpeCanvas, bb::BBox{2})
    return ipe_group(f, canvas; clip=_page_clip_path(_world_page_box(canvas, bb)))
end

function clip_to(f::Function, canvas::IpeCanvas, box::PageBox)
    return ipe_group(f, canvas; clip=_page_clip_path(box))
end

"""
    inset(f, canvas, page_box; fit, margin=8, flip_y=false, clip=true)

Draw a fitted world-coordinate scene inside a fixed page rectangle. The prior
viewport is restored after `f`, so insets compose with a fitted main figure.
"""
function inset(
    f::Function, canvas::IpeCanvas, page_box::PageBox;
    fit,
    margin::Real=8.0,
    flip_y::Bool=false,
    clip::Bool=true,
)
    items = fit isa Tuple || fit isa AbstractVector ? fit : (fit,)
    previous = canvas.viewport
    _setup_transform!(canvas, union_bbox(items...), page_box; margin=margin, flip_y=flip_y)
    try
        clip ? clip_to(f, canvas, page_box) : f(canvas)
    finally
        canvas.viewport = previous
    end
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
    style::Symbol = :math,
    offset::Tuple{<:Real,<:Real} = (0.0, 0.0),
    anchor::Union{Symbol,Nothing} = nothing
)
    pt = _apply_tf(canvas, p)
    pt += Point(Float64(offset[1]), Float64(offset[2]))
    if anchor !== nothing
        alignments = Dict(
            :center => (:center, :center),
            :north => (:center, :top), :south => (:center, :bottom),
            :east => (:right, :center), :west => (:left, :center),
            :northeast => (:right, :top), :northwest => (:left, :top),
            :southeast => (:right, :bottom), :southwest => (:left, :bottom),
        )
        haskey(alignments, anchor) || throw(ArgumentError("unknown label anchor: $anchor"))
        halign, valign = alignments[anchor]
    end
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

"""Place a LaTeX-aware label at a world point with an optional page-space offset."""
function label!(canvas::IpeCanvas, p::Point{2}, text; appearance::Union{Style,Nothing}=nothing, kwargs...)
    allowed = (:halign, :valign, :stroke, :size, :style, :offset, :anchor)
    options = _style_options(canvas, appearance, (; kwargs...), allowed)
    return draw_label!(canvas, p, text; options...)
end

label!(canvas::IpeCanvas, x::Real, y::Real, text; kwargs...) =
    label!(canvas, Point(Float64(x), Float64(y)), text; kwargs...)

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

function _anchored_box(canvas::IpeCanvas, width::Real, height::Real, position::Symbol, margin::Real)
    position in (:northwest, :north, :northeast, :west, :center, :east,
                 :southwest, :south, :southeast) ||
        throw(ArgumentError("unknown page position: $position"))
    horizontal = position in (:northwest, :west, :southwest) ? :west :
                 position in (:northeast, :east, :southeast) ? :east : :center
    vertical = position in (:northwest, :north, :northeast) ? :north :
               position in (:southwest, :south, :southeast) ? :south : :center
    x = horizontal == :west ? margin :
        horizontal == :east ? canvas.width - margin - width : (canvas.width - width) / 2
    y = vertical == :south ? margin :
        vertical == :north ? canvas.height - margin - height : (canvas.height - height) / 2
    return PageBox(x, y, width, height)
end

"""
    legend!(canvas, entries; position=:northeast, at=nothing, width=130, ...)

Draw a page-space legend. Each entry is a `label => Style` pair whose style is
shown as a rectangular swatch. Use `at=(x, y)` for an exact lower-left corner.
"""
function legend!(
    canvas::IpeCanvas, entries;
    position::Symbol=:northeast,
    at::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
    width::Real=130.0,
    row_height::Real=18.0,
    padding::Real=7.0,
    margin::Real=12.0,
    background::Symbol=:white,
    border::Symbol=:black,
    text_style::Style=Style(stroke=:black, size=:small, style=:normal),
)
    rows = collect(entries)
    all(row -> row isa Pair && last(row) isa Style, rows) ||
        throw(ArgumentError("legend entries must be label => Style pairs"))
    height = 2 * padding + row_height * length(rows)
    box = isnothing(at) ? _anchored_box(canvas, width, height, position, margin) :
                         PageBox(at[1], at[2], width, height)
    page_space(canvas) do cv
        draw_box!(cv, box.x, box.y, box.x + box.width, box.y + box.height;
                  fill=background, stroke=border)
        swatch_width = min(22.0, 0.22 * width)
        for (index, row) in enumerate(rows)
            cy = box.y + box.height - padding - (index - 0.5) * row_height
            swatch = BBox(point(box.x + padding, cy - 4),
                          point(box.x + padding + swatch_width, cy + 4))
            draw!(cv, swatch; style=last(row))
            label!(cv, point(box.x + padding + swatch_width + 7, cy), string(first(row));
                   appearance=text_style, halign=:left, valign=:center)
        end
    end
    return canvas
end

"""
    scale_bar!(canvas, length; position=:southwest, at=nothing, label=string(length), ...)

Draw a scale bar whose `length` is measured in current world units and whose
placement and tick size are fixed page units.
"""
function scale_bar!(
    canvas::IpeCanvas, length::Real;
    position::Symbol=:southwest,
    at::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
    label=string(length),
    margin::Real=12.0,
    tick::Real=4.0,
    style::Style=Style(stroke=:black, pen=:heavier),
    text_style::Style=Style(stroke=:black, size=:small, style=:normal),
)
    length > 0 || throw(ArgumentError("scale-bar length must be positive"))
    page_length = _apply_length(canvas, length)
    box = isnothing(at) ? _anchored_box(canvas, page_length, 2 * tick, position, margin) :
                         PageBox(at[1], at[2] - tick, page_length, 2 * tick)
    y = box.y + tick
    page_space(canvas) do cv
        draw!(cv, Segment(point(box.x, y), point(box.x + page_length, y)); style=style)
        draw!(cv, Segment(point(box.x, y - tick), point(box.x, y + tick)); style=style)
        draw!(cv, Segment(point(box.x + page_length, y - tick),
                          point(box.x + page_length, y + tick)); style=style)
        label === nothing || label!(cv, point(box.x + page_length / 2, y + tick + 3), label;
                                    appearance=text_style, halign=:center, valign=:bottom)
    end
    return canvas
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
    mkpath(dirname(abspath(output_pdf)))
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
    label::String = replace(basename(base_path), "_" => ":"),
    outputs = (:ipe, :pdf, :tex)
)
    requested = Set(Symbol.(outputs isa Symbol ? (outputs,) : outputs))
    all(output in (:ipe, :pdf, :tex) for output in requested) ||
        throw(ArgumentError("outputs may contain only :ipe, :pdf, and :tex"))
    ipe_path = base_path * ".ipe"
    pdf_path = base_path * ".pdf"
    tex_path = base_path * "_fig.tex"

    if :ipe in requested || :pdf in requested
        save_ipe(canvas, ipe_path)
    end
    compiled = :pdf in requested ? compile_pdf(ipe_path, pdf_path) : true
    :tex in requested &&
        save_figure_tex(tex_path; figure_name=basename(base_path), caption=caption, label=label)
    if !compiled
        @warn "PDF compilation failed; preserving Ipe source" path=ipe_path
        throw(ErrorException("PDF output was requested, but Ipe compilation failed"))
    elseif !(:ipe in requested)
        rm(ipe_path; force=true)
    end

    return (
        ipe=:ipe in requested ? ipe_path : nothing,
        pdf=:pdf in requested ? pdf_path : nothing,
        tex=:tex in requested ? tex_path : nothing,
    )
end

"""
    open_ipe(base_path; caption="", label="", preview=false, kwargs...) do canvas
        ...
    end

Convenient block syntax to construct, save, compile, and generate LaTeX fragments
in one go. With `preview=true`, open the retained `.ipe` source in Ipe.
"""
function open_ipe(
    f::Function, base_path::String;
    caption::String="",
    label::String="",
    outputs=(:ipe, :pdf, :tex),
    fit=nothing,
    margin::Real=30.0,
    flip_y::Bool=false,
    preview::Bool=false,
    kwargs...
)
    requested = Set(Symbol.(outputs isa Symbol ? (outputs,) : outputs))
    preview && !(:ipe in requested) &&
        throw(ArgumentError("preview=true requires :ipe in outputs"))
    clean_base = replace(base_path, r"(\.ipe|\.pdf|_fig\.tex)$" => "")
    canvas = IpeCanvas(; kwargs...)
    fit === nothing || fit!(canvas, fit; margin=margin, flip_y=flip_y)
    f(canvas)
    lbl = isempty(label) ? replace(basename(clean_base), "_" => ":") : label
    artifacts = export_figure(canvas, clean_base; caption=caption, label=lbl, outputs=outputs)
    if preview
        edit_ipe(artifacts.ipe)
    end
    return artifacts
end

"""
    figure(f, path; keep_source=true, source=:ipe, preview=false, kwargs...)

Generate a figure from geometry. A `.pdf` target is rendered through Ipe so that
geometry and LaTeX labels remain vector content. The editable `.ipe` source is
retained by default.
Set `preview=true` to retain and open the editable source after generation.
"""
function figure(
    f::Function, path::String;
    keep_source::Bool=true,
    source::Symbol=:ipe,
    preview::Bool=false,
    kwargs...
)
    source == :ipe || throw(ArgumentError("the only supported figure source is :ipe"))
    base, extension = splitext(path)
    if isempty(extension)
        base = path
        outputs = keep_source || preview ? (:ipe, :pdf) : (:pdf,)
    elseif extension == ".pdf"
        outputs = keep_source || preview ? (:ipe, :pdf) : (:pdf,)
    elseif extension == ".ipe"
        outputs = (:ipe,)
    else
        throw(ArgumentError("figure target must have extension .pdf or .ipe"))
    end
    return open_ipe(f, base; outputs=outputs, preview=preview, kwargs...)
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
