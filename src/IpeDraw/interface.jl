# -----------------------------------------------------------------------------
# Compact Geometry-Aware Drawing Interface
# -----------------------------------------------------------------------------

Base.@constprop :aggressive @inline _filter_keys(::Tuple{}, allowed, ::Val) = ()
Base.@constprop :aggressive @inline function _filter_keys(
    ks::Tuple{Symbol,Vararg{Symbol}},
    allowed,
    ::Val{keep},
) where {keep}
    k = ks[1]
    rest = _filter_keys(Base.tail(ks), allowed, Val(keep))
    return (k in allowed) == keep ? (k, rest...) : rest
end

function _style_options(
    active::Style{A},
    style::S,
    overrides::O,
    ::Val{allowed},
) where {A,S<:Union{Nothing,Style},O<:NamedTuple,allowed}
    unknown = _filter_keys(keys(overrides), allowed, Val(false))
    isempty(unknown) ||
        throw(ArgumentError("unsupported drawing options: $(join(unknown, ", "))"))
    local_style = style === nothing ? NamedTuple() : style.values
    merged = merge(active.values, local_style, overrides)
    ks = _filter_keys(keys(merged), allowed, Val(true))
    return NamedTuple{ks}(merged)
end

function _opacity_name(value::Real)
    0 < value < 1 ||
        throw(ArgumentError("named Ipe opacity must be strictly between zero and one"))
    percent = round(Int, 100 * value)
    percent % 10 == 0 || throw(ArgumentError("Ipe opacity must be a multiple of 0.1"))
    return Symbol("$(percent)%")
end

_opacity_name(value::Union{Symbol,String}) = value

function _draw_fill_opacity!(
    renderer::F,
    canvas::IpeCanvas,
    object,
    options::O,
) where {F,O<:NamedTuple}
    fill_opacity = get(options, :fill_opacity, nothing)
    base_keys = _filter_keys(keys(options), (:fill_opacity,), Val(false))
    base = NamedTuple{base_keys}(options)
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

const POINT_STYLE = (:stroke, :fill, :size, :shape)
const PATH_STYLE = (:stroke, :fill, :pen, :dash, :opacity, :fill_opacity)
const ELLIPSE_STYLE = (PATH_STYLE..., :tiling)
const POLYGON_STYLE = (PATH_STYLE..., :close)
const BEZIER_STYLE = (:stroke, :fill, :pen, :dash, :arrow, :rarrow)
const CURVE_STYLE = (:stroke, :pen, :dash, :arrow, :rarrow)

function _draw_styled!(
    renderer::F,
    canvas::IpeCanvas,
    object,
    active::Style{A},
    style::S,
    overrides::O,
    ::Val{allowed},
    ::Val{split_opacity},
) where {F,A,S<:Union{Nothing,Style},O<:NamedTuple,allowed,split_opacity}
    options = _style_options(active, style, overrides, Val(allowed))
    if split_opacity
        return _draw_fill_opacity!(renderer, canvas, object, options)
    end
    return renderer(canvas, object; options...)
end

@inline function _dispatch_active(
    fn::F,
    canvas::IpeCanvas,
    object,
    style::S,
    overrides::O,
) where {F,S<:Union{Nothing,Style},O<:NamedTuple}
    active = canvas.active_style
    if active isa Style{@NamedTuple{}}
        return fn(canvas, object, active, style, overrides)
    end
    return fn(canvas, object, active, style, overrides)
end

@inline _draw_spec(::Point{2}) = (draw_point!, Val(POINT_STYLE), Val(false))
@inline _draw_spec(::Segment{2}) = (draw_segment!, Val(CURVE_STYLE), Val(false))
@inline _draw_spec(::BBox{2}) = (draw_box!, Val(PATH_STYLE), Val(true))
@inline _draw_spec(::Sphere{2}) = (draw_circle!, Val(PATH_STYLE), Val(true))
@inline _draw_spec(::CircleArc) = (draw_arc!, Val(CURVE_STYLE), Val(false))
@inline _draw_spec(::Ellipse) = (draw_ellipse!, Val(ELLIPSE_STYLE), Val(true))
@inline _draw_spec(::EllipticArc) = (draw_arc!, Val(CURVE_STYLE), Val(false))
@inline _draw_spec(::CubicBezier{2}) = (draw_bezier!, Val(BEZIER_STYLE), Val(false))
@inline _draw_spec(::CubicSpline{2}) = (draw_spline!, Val(BEZIER_STYLE), Val(false))
@inline _draw_spec(::AbsPntSeq{2}) = (draw_polygon!, Val(POLYGON_STYLE), Val(true))
@inline _draw_spec(::AbstractVector{<:Point{2}}) =
    (draw_polygon!, Val(POLYGON_STYLE), Val(true))
@inline _draw_spec(::Any) = nothing

function _draw_each!(
    canvas::IpeCanvas,
    objects::AbstractVector,
    active::Style{A},
    style::S,
    overrides::O,
) where {A,S<:Union{Nothing,Style},O<:NamedTuple}
    for object in objects
        _draw_one!(canvas, object, active, style, overrides)
    end
    return canvas
end

function _draw_one!(
    canvas::IpeCanvas,
    object,
    active::Style{A},
    style::S,
    overrides::O,
) where {A,S<:Union{Nothing,Style},O<:NamedTuple}
    spec = _draw_spec(object)
    if spec !== nothing
        renderer, allowed, split_opacity = spec
        return _draw_styled!(
            renderer,
            canvas,
            object,
            active,
            style,
            overrides,
            allowed,
            split_opacity,
        )
    end
    object isa AbstractVector &&
        return _draw_each!(canvas, object, active, style, overrides)
    return style === nothing ? draw!(canvas, object; overrides...) :
           draw!(canvas, object; style = style, overrides...)
end

_mark_many!(
    canvas::IpeCanvas,
    points,
    active::Style{A},
    style::S,
    overrides::O,
) where {A,S<:Union{Nothing,Style},O<:NamedTuple} = _draw_styled!(
    draw_points!,
    canvas,
    points,
    active,
    style,
    overrides,
    Val(POINT_STYLE),
    Val(false),
)

"""Draw a geometric object using multiple dispatch and optional reusable `Style`."""
function draw!(
    canvas::IpeCanvas,
    p::Point{2};
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    return _dispatch_active(_draw_one!, canvas, p, style, (; kwargs...))
end

function draw!(
    canvas::IpeCanvas,
    s::Segment{2};
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    return _dispatch_active(_draw_one!, canvas, s, style, (; kwargs...))
end

function draw!(
    canvas::IpeCanvas,
    bb::BBox{2};
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    return _dispatch_active(_draw_one!, canvas, bb, style, (; kwargs...))
end

function draw!(
    canvas::IpeCanvas,
    circle::Sphere{2};
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    return _dispatch_active(_draw_one!, canvas, circle, style, (; kwargs...))
end

function draw!(
    canvas::IpeCanvas,
    arc::CircleArc;
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    return _dispatch_active(_draw_one!, canvas, arc, style, (; kwargs...))
end

function draw!(
    canvas::IpeCanvas,
    ellipse::Ellipse;
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    return _dispatch_active(_draw_one!, canvas, ellipse, style, (; kwargs...))
end

function draw!(
    canvas::IpeCanvas,
    arc::EllipticArc;
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    return _dispatch_active(_draw_one!, canvas, arc, style, (; kwargs...))
end

function draw!(
    canvas::IpeCanvas,
    curve::CubicBezier{2};
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    return _dispatch_active(_draw_one!, canvas, curve, style, (; kwargs...))
end

function draw!(
    canvas::IpeCanvas,
    curve::CubicSpline{2};
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    return _dispatch_active(_draw_one!, canvas, curve, style, (; kwargs...))
end

function draw!(
    canvas::IpeCanvas,
    polygon::AbsPntSeq{2};
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    return _dispatch_active(_draw_one!, canvas, polygon, style, (; kwargs...))
end

function draw!(
    canvas::IpeCanvas,
    polygon::AbstractVector{<:Point{2}};
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    return _dispatch_active(_draw_one!, canvas, polygon, style, (; kwargs...))
end

function draw!(
    canvas::IpeCanvas,
    objects::AbstractVector;
    style::Union{Style,Nothing} = nothing,
    kwargs...,
)
    return _dispatch_active(_draw_each!, canvas, objects, style, (; kwargs...))
end

"""Draw one point or a collection of points as fixed-size page marks."""
mark!(canvas::IpeCanvas, p::Point{2}; kwargs...) = draw!(canvas, p; kwargs...)

function mark!(canvas::IpeCanvas, points; style::Union{Style,Nothing} = nothing, kwargs...)
    return _dispatch_active(_mark_many!, canvas, points, style, (; kwargs...))
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

function _merge_style(
    previous::Style{A},
    style::Style{S},
    overrides::O,
) where {A,S,O<:NamedTuple}
    return Style(merge(previous.values, style.values, overrides))
end

"""Apply style defaults within `f` and restore the previous defaults afterward."""
function with_style(f::Function, canvas::IpeCanvas, style::Style = Style(); kwargs...)
    previous = canvas.active_style
    canvas.active_style = _merge_style(previous, style, (; kwargs...))
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
