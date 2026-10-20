"""
    IpeDraw

A clean, robust Julia interface for programmatically generating Ipe 7 vector figures.
Integrates natively with `BasicCompGeometry` geometric primitives (Point, Segment, BBox, Polygon)
and provides high-level helpers for algorithmic and conceptual diagrams.
"""
module IpeDraw

using ..BasicCompGeometry
import ..BasicCompGeometry: cardin
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

Base.propertynames(theme::Theme, private::Bool = false) =
    private ? (:styles, keys(theme.styles)...) : Tuple(keys(theme.styles))

"""Return restrained defaults suitable for papers and lecture notes."""
function publication_theme()
    return Theme(
        region = Style(fill_opacity = 0.2, pen = :heavier),
        boundary = Style(stroke = :black, pen = :heavier),
        point = Style(stroke = :black, fill = :black, size = :normal),
        annotation = Style(stroke = :black, size = :small),
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
    stylesheet_template::Union{String,Nothing}
    active_layer::String
    layers::Vector{String}
    views::Vector{String}
    elements::Vector{String}
    viewport::Union{Nothing,Viewport}
    active_style::Style

    function IpeCanvas(;
        width::Real = 576.0,
        height::Real = 504.0,
        paper::Union{String,Nothing} = nothing,
        bbox::String = "cropbox",
        preamble::String = "\\usepackage{amsmath,amssymb}\\def\\ipeMode{TRUE}\\def\\Sample{\\mathsf{R}}",
        template::Union{String,Nothing} = nothing,
        layer::String = "alpha",
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
            Style(),
        )
    end
end

# -----------------------------------------------------------------------------
# Attribute Formatting Helpers
# -----------------------------------------------------------------------------

function _fmt_attr(name::String, val::Union{Symbol,String,Nothing})
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

_apply_tf(canvas::IpeCanvas, x::Real, y::Real) =
    _apply_tf(canvas, Point(Float64(x), Float64(y)))
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
function add_view!(canvas::IpeCanvas, layers::Vector{String}; active = nothing)
    act = active === nothing ? layers[1] : string(active)
    push!(canvas.views, "<view layers=\"$(join(layers, " "))\" active=\"$act\"/>")
    return canvas
end

"""
    add_preamble!(canvas, snippet)

Append LaTeX packages or macros to the canvas LaTeX preamble.
"""
function add_preamble!(canvas::IpeCanvas, snippet::AbstractString)
    canvas.preamble =
        isempty(canvas.preamble) ? String(snippet) :
        canvas.preamble * "\n" * String(snippet)
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
    canvas::IpeCanvas,
    world_bb::BBox{2,T};
    margin::Real = 30.0,
    flip_y::Bool = false,
) where {T}
    return _setup_transform!(
        canvas,
        world_bb,
        PageBox(0, 0, canvas.width, canvas.height);
        margin = margin,
        flip_y = flip_y,
    )
end

function _setup_transform!(
    canvas::IpeCanvas,
    world_bb::BBox{2,T},
    page_box::PageBox;
    margin::Real = 0.0,
    flip_y::Bool = false,
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
    setup_transform!(canvas, union_bbox(items...); margin = margin, flip_y = flip_y)
end

include("IpeDraw/primitives.jl")
include("IpeDraw/interface.jl")
include("IpeDraw/diagrams.jl")
include("IpeDraw/output.jl")

BasicCompGeometry.width(c::IpeCanvas) = c.width
BasicCompGeometry.height(c::IpeCanvas) = c.height

end # module IpeDraw
