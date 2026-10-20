# -----------------------------------------------------------------------------
# Compact Geometry-Aware Drawing Interface
# -----------------------------------------------------------------------------

function _style_options(
    canvas::IpeCanvas,
    style::Union{Style,Nothing},
    overrides::NamedTuple,
    allowed,
)
    unknown = filter(key -> !(key in allowed), keys(overrides))
    isempty(unknown) ||
        throw(ArgumentError("unsupported drawing options: $(join(unknown, ", "))"))
    local_style = style === nothing ? NamedTuple() : style.values
    merged = merge(canvas.active_style.values, local_style, overrides)
    return (; (key => value for (key, value) in pairs(merged) if key in allowed)...)
end

function _opacity_name(value::Real)
    0 < value < 1 ||
        throw(ArgumentError("named Ipe opacity must be strictly between zero and one"))
    percent = round(Int, 100 * value)
    percent % 10 == 0 || throw(ArgumentError("Ipe opacity must be a multiple of 0.1"))
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
        opacity =
            fill_opacity isa Real && isone(fill_opacity) ? nothing :
            _opacity_name(fill_opacity)
        fill_options = merge(base, (stroke = :none, opacity = opacity))
        renderer(canvas, object; fill_options...)
    end
    if stroke != :none
        stroke_options = merge(base, (fill = :none,))
        renderer(canvas, object; stroke_options...)
    end
    return canvas
end

const PATH_STYLE = (:stroke, :fill, :pen, :dash, :opacity, :fill_opacity)
const CURVE_STYLE = (:stroke, :pen, :dash, :arrow, :rarrow)

"""Draw a geometric object using multiple dispatch and optional reusable `Style`."""
function draw!(
    canvas::IpeCanvas,
    p::Point{2};
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    options = _style_options(canvas, style, (; kwargs...), (:stroke, :fill, :size, :shape))
    return draw_point!(canvas, p; options...)
end

function draw!(
    canvas::IpeCanvas,
    s::Segment{2};
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    options = _style_options(canvas, style, (; kwargs...), CURVE_STYLE)
    return draw_segment!(canvas, s; options...)
end

function draw!(
    canvas::IpeCanvas,
    bb::BBox{2};
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    options = _style_options(canvas, style, (; kwargs...), PATH_STYLE)
    return _draw_fill_opacity!(draw_box!, canvas, bb, options)
end

function draw!(
    canvas::IpeCanvas,
    circle::Sphere{2};
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    options = _style_options(canvas, style, (; kwargs...), PATH_STYLE)
    return _draw_fill_opacity!(draw_circle!, canvas, circle, options)
end

function draw!(
    canvas::IpeCanvas,
    arc::CircleArc;
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    options = _style_options(canvas, style, (; kwargs...), CURVE_STYLE)
    return draw_arc!(canvas, arc; options...)
end

function draw!(
    canvas::IpeCanvas,
    ellipse::Ellipse;
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    allowed = (PATH_STYLE..., :tiling)
    options = _style_options(canvas, style, (; kwargs...), allowed)
    return _draw_fill_opacity!(draw_ellipse!, canvas, ellipse, options)
end

function draw!(
    canvas::IpeCanvas,
    arc::EllipticArc;
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    options = _style_options(canvas, style, (; kwargs...), CURVE_STYLE)
    return draw_arc!(canvas, arc; options...)
end

function draw!(
    canvas::IpeCanvas,
    curve::CubicBezier{2};
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    options = _style_options(
        canvas,
        style,
        (; kwargs...),
        (:stroke, :fill, :pen, :dash, :arrow, :rarrow),
    )
    return draw_bezier!(canvas, curve; options...)
end

function draw!(
    canvas::IpeCanvas,
    curve::CubicSpline{2};
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    options = _style_options(
        canvas,
        style,
        (; kwargs...),
        (:stroke, :fill, :pen, :dash, :arrow, :rarrow),
    )
    return draw_spline!(canvas, curve; options...)
end

function draw!(
    canvas::IpeCanvas,
    polygon::AbsPntSeq{2};
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    options = _style_options(canvas, style, (; kwargs...), (PATH_STYLE..., :close))
    return _draw_fill_opacity!(draw_polygon!, canvas, polygon, options)
end

function draw!(
    canvas::IpeCanvas,
    polygon::AbstractVector{<:Point{2}};
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
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

function mark!(canvas::IpeCanvas, points; style::Union{Style,Nothing} = nothing, kwargs...)
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
function with_style(f::Function, canvas::IpeCanvas, style::Style = Style(); kwargs...)
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
    return @sprintf(
        "%.3f %.3f m %.3f %.3f l %.3f %.3f l %.3f %.3f l h",
        x1,
        y1,
        x2,
        y1,
        x2,
        y2,
        x1,
        y2
    )
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
    return ipe_group(f, canvas; clip = _page_clip_path(_world_page_box(canvas, bb)))
end

function clip_to(f::Function, canvas::IpeCanvas, box::PageBox)
    return ipe_group(f, canvas; clip = _page_clip_path(box))
end

"""
    inset(f, canvas, page_box; fit, margin=8, flip_y=false, clip=true)

Draw a fitted world-coordinate scene inside a fixed page rectangle. The prior
viewport is restored after `f`, so insets compose with a fitted main figure.
"""
function inset(
    f::Function,
    canvas::IpeCanvas,
    page_box::PageBox;
    fit,
    margin::Real = 8.0,
    flip_y::Bool = false,
    clip::Bool = true,
)
    items = fit isa Tuple || fit isa AbstractVector ? fit : (fit,)
    previous = canvas.viewport
    _setup_transform!(
        canvas,
        union_bbox(items...),
        page_box;
        margin = margin,
        flip_y = flip_y,
    )
    try
        clip ? clip_to(f, canvas, page_box) : f(canvas)
    finally
        canvas.viewport = previous
    end
    return canvas
end
