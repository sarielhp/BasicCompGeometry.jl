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
    canvas::IpeCanvas,
    p::Point{2},
    text;
    halign::Symbol = :center,
    valign::Symbol = :baseline,
    stroke::Symbol = :black,
    size::Symbol = :normal,
    style::Symbol = :math,
    offset::Tuple{<:Real,<:Real} = (0.0, 0.0),
    anchor::Union{Symbol,Nothing} = nothing,
)
    pt = _apply_tf(canvas, p)
    pt += Point(Float64(offset[1]), Float64(offset[2]))
    if anchor !== nothing
        alignments = Dict(
            :center => (:center, :center),
            :north => (:center, :top),
            :south => (:center, :bottom),
            :east => (:right, :center),
            :west => (:left, :center),
            :northeast => (:right, :top),
            :northwest => (:left, :top),
            :southeast => (:right, :bottom),
            :southwest => (:left, :bottom),
        )
        haskey(alignments, anchor) || throw(ArgumentError("unknown label anchor: $anchor"))
        halign, valign = alignments[anchor]
    end
    clean_text = string(text)
    # If math style and text not already wrapped in $, wrap it
    if style == :math &&
       !startswith(strip(clean_text), "\$") &&
       !startswith(strip(clean_text), "\\begin")
        clean_text = "\$" * clean_text * "\$"
    end
    xml_text = _escape_ipe_xml(clean_text)
    xml = "<text layer=\"$(canvas.active_layer)\" pos=\"$(_pt_str(pt))\" stroke=\"$stroke\" type=\"label\" size=\"$size\" halign=\"$halign\" valign=\"$valign\">$xml_text</text>"
    push!(canvas.elements, xml)
    return canvas
end

draw_label!(canvas::IpeCanvas, x::Real, y::Real, text; kwargs...) =
    draw_label!(canvas, Point(Float64(x), Float64(y)), text; kwargs...)

"""Place a LaTeX-aware label at a world point with an optional page-space offset."""
function label!(
    canvas::IpeCanvas,
    p::Point{2},
    text;
    appearance::Union{Style,Nothing} = nothing,
    kwargs...,
)
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
    canvas::IpeCanvas,
    x1::Real,
    x2::Real,
    y::Real;
    height::Real = 18.0,
    stroke::Symbol = :black,
    fill::Symbol = :gray7,
    pen::Symbol = :heavier,
    label_left = nothing,
    label_right = nothing,
)
    y_bot = y - height / 2.0
    y_top = y + height / 2.0
    draw_box!(canvas, x1, y_bot, x2, y_top; stroke = stroke, fill = fill, pen = pen)

    if label_left !== nothing
        draw_label!(
            canvas,
            x1,
            y_bot - 12.0,
            string(label_left);
            size = :small,
            halign = :center,
            stroke = :gray2,
        )
    end
    if label_right !== nothing
        draw_label!(
            canvas,
            x2,
            y_bot - 12.0,
            string(label_right);
            size = :small,
            halign = :center,
            stroke = :gray2,
        )
    end
    return canvas
end

"""
    draw_span!(canvas, x1, x2, y; height=18.0, fill=:lightgreen, stroke=:darkgreen, pen=:heavier, dash=:solid)

Highlight an interval or active span within a bar.
"""
function draw_span!(
    canvas::IpeCanvas,
    x1::Real,
    x2::Real,
    y::Real;
    height::Real = 18.0,
    fill::Symbol = :lightgreen,
    stroke::Symbol = :darkgreen,
    pen::Symbol = :heavier,
    dash::Symbol = :solid,
)
    y_bot = y - height / 2.0
    y_top = y + height / 2.0
    draw_box!(
        canvas,
        x1,
        y_bot,
        x2,
        y_top;
        fill = fill,
        stroke = stroke,
        pen = pen,
        dash = dash,
    )
    return canvas
end

"""
    draw_dimension!(canvas, p1, p2; label=nothing, arrow=:both, stroke=:darkgreen, pen=:normal, label_offset=12.0)
    draw_dimension!(canvas, x1, x2, y; label=nothing, ...)

Draw a dimension line (double-arrow) with an optional text label.
"""
function draw_dimension!(
    canvas::IpeCanvas,
    p1::Point{2},
    p2::Point{2};
    label = nothing,
    arrow::Symbol = :both,
    stroke::Symbol = :darkgreen,
    pen::Symbol = :normal,
    label_offset::Real = 12.0,
)
    arr = arrow == :both ? :pointed : arrow
    rarr = arrow == :both ? :pointed : :none
    draw_segment!(canvas, p1, p2; stroke = stroke, pen = pen, arrow = arr, rarrow = rarr)

    if label !== nothing
        mid = Point((p1[1] + p2[1]) / 2.0, (p1[2] + p2[2]) / 2.0)
        dx = p2[1] - p1[1]
        dy = p2[2] - p1[2]
        len = hypot(dx, dy)
        normal = len > 0 ? Point(-dy / len, dx / len) : Point(0.0, 1.0)
        pos = Point(mid[1] + normal[1] * label_offset, mid[2] + normal[2] * label_offset)
        draw_label!(canvas, pos, label; stroke = stroke, size = :small, halign = :center)
    end
    return canvas
end

draw_dimension!(
    canvas::IpeCanvas,
    x1::Real,
    x2::Real,
    y::Real;
    label = nothing,
    label_offset::Real = -12.0,
    kwargs...,
) = draw_dimension!(
    canvas,
    Point(Float64(x1), Float64(y)),
    Point(Float64(x2), Float64(y));
    label = label,
    label_offset = label_offset,
    kwargs...,
)

"""
    draw_arrow!(canvas, p1, p2; stroke=:darkgreen, pen=:fat, dash=:dashed, arrow=:pointed)

Draw a directed arrow (e.g. lifting arrow between domains).
"""
function draw_arrow!(
    canvas::IpeCanvas,
    p1::Point{2},
    p2::Point{2};
    stroke::Symbol = :darkgreen,
    pen::Symbol = :fat,
    dash::Symbol = :dashed,
    arrow::Symbol = :pointed,
)
    draw_segment!(canvas, p1, p2; stroke = stroke, pen = pen, dash = dash, arrow = arrow)
end

draw_arrow!(canvas::IpeCanvas, x1::Real, y1::Real, x2::Real, y2::Real; kwargs...) =
    draw_arrow!(
        canvas,
        Point(Float64(x1), Float64(y1)),
        Point(Float64(x2), Float64(y2));
        kwargs...,
    )

"""
    draw_curved_arrow!(canvas, p1, control_pt, p2; stroke=:blue, pen=:fat, arrow=:pointed)
    draw_curved_arrow!(canvas, p1, p2; bend=20.0, stroke=:blue, pen=:fat, arrow=:pointed)

Draw a curved Bézier arrow connecting two points.
"""
function draw_curved_arrow!(
    canvas::IpeCanvas,
    p1::Point{2},
    control_pt::Point{2},
    p2::Point{2};
    stroke::Symbol = :blue,
    pen::Symbol = :fat,
    arrow::Symbol = :pointed,
)
    q1 = _apply_tf(canvas, p1)
    qc = _apply_tf(canvas, control_pt)
    q2 = _apply_tf(canvas, p2)
    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("pen", pen))$(_fmt_attr("arrow", arrow))"
    xml =
        "<path layer=\"$(canvas.active_layer)\"$attrs>\n" *
        "$(_pt_str(q1)) m\n$(_pt_str(qc))\n$(_pt_str(q2)) c\n</path>"
    push!(canvas.elements, xml)
    return canvas
end

function draw_curved_arrow!(
    canvas::IpeCanvas,
    p1::Point{2},
    p2::Point{2};
    bend::Real = 20.0,
    stroke::Symbol = :blue,
    pen::Symbol = :fat,
    arrow::Symbol = :pointed,
)
    mid = Point((p1[1] + p2[1]) / 2.0, (p1[2] + p2[2]) / 2.0)
    dx = p2[1] - p1[1]
    dy = p2[2] - p1[2]
    len = hypot(dx, dy)
    normal = len > 0 ? Point(-dy / len, dx / len) : Point(0.0, 1.0)
    ctrl = Point(mid[1] + normal[1] * bend, mid[2] + normal[2] * bend)
    return draw_curved_arrow!(
        canvas,
        p1,
        ctrl,
        p2;
        stroke = stroke,
        pen = pen,
        arrow = arrow,
    )
end

function _anchored_box(
    canvas::IpeCanvas,
    width::Real,
    height::Real,
    position::Symbol,
    margin::Real,
)
    position in (
        :northwest,
        :north,
        :northeast,
        :west,
        :center,
        :east,
        :southwest,
        :south,
        :southeast,
    ) || throw(ArgumentError("unknown page position: $position"))
    horizontal =
        position in (:northwest, :west, :southwest) ? :west :
        position in (:northeast, :east, :southeast) ? :east : :center
    vertical =
        position in (:northwest, :north, :northeast) ? :north :
        position in (:southwest, :south, :southeast) ? :south : :center
    x =
        horizontal == :west ? margin :
        horizontal == :east ? canvas.width - margin - width : (canvas.width - width) / 2
    y =
        vertical == :south ? margin :
        vertical == :north ? canvas.height - margin - height : (canvas.height - height) / 2
    return PageBox(x, y, width, height)
end

"""
    legend!(canvas, entries; position=:northeast, at=nothing, width=130, ...)

Draw a page-space legend. Each entry is a `label => Style` pair whose style is
shown as a rectangular swatch. Use `at=(x, y)` for an exact lower-left corner.
"""
function legend!(
    canvas::IpeCanvas,
    entries;
    position::Symbol = :northeast,
    at::Union{Nothing,Tuple{<:Real,<:Real}} = nothing,
    width::Real = 130.0,
    row_height::Real = 18.0,
    padding::Real = 7.0,
    margin::Real = 12.0,
    background::Symbol = :white,
    border::Symbol = :black,
    text_style::Style = Style(stroke = :black, size = :small, style = :normal),
)
    rows = collect(entries)
    all(row -> row isa Pair && last(row) isa Style, rows) ||
        throw(ArgumentError("legend entries must be label => Style pairs"))
    height = 2 * padding + row_height * length(rows)
    box =
        isnothing(at) ? _anchored_box(canvas, width, height, position, margin) :
        PageBox(at[1], at[2], width, height)
    page_space(canvas) do cv
        draw_box!(
            cv,
            box.x,
            box.y,
            box.x + box.width,
            box.y + box.height;
            fill = background,
            stroke = border,
        )
        swatch_width = min(22.0, 0.22 * width)
        for (index, row) in enumerate(rows)
            cy = box.y + box.height - padding - (index - 0.5) * row_height
            swatch = BBox(
                point(box.x + padding, cy - 4),
                point(box.x + padding + swatch_width, cy + 4),
            )
            draw!(cv, swatch; style = last(row))
            label!(
                cv,
                point(box.x + padding + swatch_width + 7, cy),
                string(first(row));
                appearance = text_style,
                halign = :left,
                valign = :center,
            )
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
    canvas::IpeCanvas,
    length::Real;
    position::Symbol = :southwest,
    at::Union{Nothing,Tuple{<:Real,<:Real}} = nothing,
    label = string(length),
    margin::Real = 12.0,
    tick::Real = 4.0,
    style::Style = Style(stroke = :black, pen = :heavier),
    text_style::Style = Style(stroke = :black, size = :small, style = :normal),
)
    length > 0 || throw(ArgumentError("scale-bar length must be positive"))
    page_length = _apply_length(canvas, length)
    box =
        isnothing(at) ? _anchored_box(canvas, page_length, 2 * tick, position, margin) :
        PageBox(at[1], at[2] - tick, page_length, 2 * tick)
    y = box.y + tick
    page_space(canvas) do cv
        draw!(cv, Segment(point(box.x, y), point(box.x + page_length, y)); style = style)
        draw!(cv, Segment(point(box.x, y - tick), point(box.x, y + tick)); style = style)
        draw!(
            cv,
            Segment(
                point(box.x + page_length, y - tick),
                point(box.x + page_length, y + tick),
            );
            style = style,
        )
        label === nothing || label!(
            cv,
            point(box.x + page_length / 2, y + tick + 3),
            label;
            appearance = text_style,
            halign = :center,
            valign = :bottom,
        )
    end
    return canvas
end
