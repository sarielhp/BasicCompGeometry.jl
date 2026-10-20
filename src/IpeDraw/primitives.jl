# -----------------------------------------------------------------------------
# Basic Geometry Drawing Primitives
# -----------------------------------------------------------------------------

"""
    draw_point!(canvas, p; stroke=:black, fill=:black, size=:normal, shape=:disk)

Draw a point as an Ipe symbol mark (e.g. `mark/disk(sx)`).
"""
function draw_point!(
    canvas::IpeCanvas,
    p::Point{2};
    stroke::Symbol = :black,
    fill::Symbol = :black,
    size::Symbol = :normal,
    shape::Symbol = :disk,
)
    pt = _apply_tf(canvas, p)
    mark_name =
        shape == :disk ? "mark/disk(sx)" :
        (shape == :circle ? "mark/circle(sx)" : "mark/box(sx)")
    xml = "<use layer=\"$(canvas.active_layer)\" name=\"$mark_name\" pos=\"$(_pt_str(pt))\" size=\"$size\"$(_fmt_attr("stroke", stroke))$(_fmt_attr("fill", fill))/>"
    push!(canvas.elements, xml)
    return canvas
end

draw_point!(canvas::IpeCanvas, x::Real, y::Real; kwargs...) =
    draw_point!(canvas, Point(Float64(x), Float64(y)); kwargs...)

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
    canvas::IpeCanvas,
    p1::Point{2},
    p2::Point{2};
    stroke::Symbol = :black,
    pen::Symbol = :normal,
    dash::Symbol = :solid,
    arrow::Symbol = :none,
    rarrow::Symbol = :none,
)
    q1 = _apply_tf(canvas, p1)
    q2 = _apply_tf(canvas, p2)
    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("pen", pen))$(_fmt_attr("dash", dash))$(_fmt_attr("arrow", arrow))$(_fmt_attr("rarrow", rarrow))"
    xml = "<path layer=\"$(canvas.active_layer)\"$attrs>\n$(_pt_str(q1)) m\n$(_pt_str(q2)) l\n</path>"
    push!(canvas.elements, xml)
    return canvas
end

draw_segment!(canvas::IpeCanvas, s::Segment{2}; kwargs...) =
    draw_segment!(canvas, s.p, s.q; kwargs...)
draw_segment!(canvas::IpeCanvas, x1::Real, y1::Real, x2::Real, y2::Real; kwargs...) =
    draw_segment!(
        canvas,
        Point(Float64(x1), Float64(y1)),
        Point(Float64(x2), Float64(y2));
        kwargs...,
    )

"""
    draw_box!(canvas, bb; stroke=:black, fill=:none, pen=:normal, dash=:solid)

Draw a rectangular box (axis-aligned bounding box).
"""
function draw_box!(
    canvas::IpeCanvas,
    x1::Real,
    y1::Real,
    x2::Real,
    y2::Real;
    stroke::Symbol = :black,
    fill::Symbol = :none,
    pen::Symbol = :normal,
    dash::Symbol = :solid,
    opacity::Union{Symbol,String,Nothing} = nothing,
)
    q1 = _apply_tf(canvas, x1, y1)
    q2 = _apply_tf(canvas, x2, y2)
    min_x, max_x = min(q1[1], q2[1]), max(q1[1], q2[1])
    min_y, max_y = min(q1[2], q2[2]), max(q1[2], q2[2])

    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("fill", fill))$(_fmt_attr("pen", pen))$(_fmt_attr("dash", dash))$(_fmt_attr("opacity", opacity))"
    xml =
        "<path layer=\"$(canvas.active_layer)\"$attrs>\n" * @sprintf(
            "%.3f %.3f m\n%.3f %.3f l\n%.3f %.3f l\n%.3f %.3f l\nh\n</path>",
            min_x,
            min_y,
            max_x,
            min_y,
            max_x,
            max_y,
            min_x,
            max_y
        )
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
    canvas::IpeCanvas,
    poly;
    close::Bool = true,
    stroke::Symbol = :black,
    fill::Symbol = :none,
    pen::Symbol = :normal,
    dash::Symbol = :solid,
    opacity::Union{Symbol,String,Nothing} = nothing,
)
    pts = [_apply_tf(canvas, p) for p in poly]
    length(pts) < 2 && return canvas

    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("fill", fill))$(_fmt_attr("pen", pen))$(_fmt_attr("dash", dash))$(_fmt_attr("opacity", opacity))"
    buf = IOBuffer()
    println(buf, "<path layer=\"$(canvas.active_layer)\"$attrs>")
    println(buf, "$(_pt_str(pts[1])) m")
    for i = 2:length(pts)
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
    canvas::IpeCanvas,
    center::Point{2},
    radius::Real;
    stroke::Symbol = :black,
    fill::Symbol = :none,
    pen::Symbol = :normal,
    dash::Symbol = :solid,
    opacity::Union{Symbol,String,Nothing} = nothing,
)
    c = _apply_tf(canvas, center)
    r = _apply_length(canvas, radius)
    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("fill", fill))$(_fmt_attr("pen", pen))$(_fmt_attr("dash", dash))$(_fmt_attr("opacity", opacity))"
    xml =
        "<path layer=\"$(canvas.active_layer)\"$attrs>\n" *
        @sprintf("%.3f 0 0 %.3f %.3f %.3f e\n</path>", r, r, c.x, c.y)
    push!(canvas.elements, xml)
    return canvas
end

draw_circle!(canvas::IpeCanvas, c::Sphere{2}; kwargs...) =
    draw_circle!(canvas, c.center, c.radius; kwargs...)
draw_circle!(canvas::IpeCanvas, x::Real, y::Real, radius::Real; kwargs...) =
    draw_circle!(canvas, Point(Float64(x), Float64(y)), radius; kwargs...)

"""
    draw_arc!(canvas, center, radius, a1, a2; stroke=:black, pen=:normal, arrow=:none, rarrow=:none)
    draw_arc!(canvas, arc::CircleArc; ...)

Draw a circular arc from angle `a1` to `a2` (in radians).
"""
function draw_arc!(
    canvas::IpeCanvas,
    center::Point{2},
    radius::Real,
    a1::Real,
    a2::Real;
    stroke::Symbol = :black,
    pen::Symbol = :normal,
    arrow::Symbol = :none,
    rarrow::Symbol = :none,
)
    c = _apply_tf(canvas, center)
    r = _apply_length(canvas, radius)
    direction = canvas.viewport !== nothing && canvas.viewport.flip_y ? -1.0 : 1.0
    p1 = Point(c.x + r * cos(a1), c.y + direction * r * sin(a1))
    p2 = Point(c.x + r * cos(a2), c.y + direction * r * sin(a2))
    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("pen", pen))$(_fmt_attr("arrow", arrow))$(_fmt_attr("rarrow", rarrow))"
    xml =
        "<path layer=\"$(canvas.active_layer)\"$attrs>\n" * @sprintf(
            "%.3f 0 0 %.3f %.3f %.3f %.3f %.3f %.3f %.3f arc\n</path>",
            r,
            direction * r,
            c.x,
            c.y,
            p1.x,
            p1.y,
            p2.x,
            p2.y
        )
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
    canvas::IpeCanvas,
    center::Point{2},
    r_major::Real,
    r_minor::Real;
    angle::Real = 0.0,
    stroke::Symbol = :black,
    fill::Symbol = :none,
    pen::Symbol = :normal,
    dash::Symbol = :none,
    opacity::Union{Symbol,String,Nothing} = nothing,
    tiling::Symbol = :none,
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
    xml =
        "<path layer=\"$(canvas.active_layer)\"$attrs>\n" *
        @sprintf("%.4f %.4f %.4f %.4f %.3f %.3f e\n</path>", m11, m21, m12, m22, c[1], c[2])
    push!(canvas.elements, xml)
    return canvas
end

draw_ellipse!(canvas::IpeCanvas, e::Ellipse; kwargs...) =
    draw_ellipse!(canvas, e.center, e.r_major, e.r_minor; angle = e.angle, kwargs...)

"""
    draw_elliptic_arc!(canvas, center, r_major, r_minor, angle, a1, a2; ...)
    draw_arc!(canvas, arc::EllipticArc; ...)

Draw an elliptic arc from parametric angle `a1` to `a2` (in radians).
"""
function draw_elliptic_arc!(
    canvas::IpeCanvas,
    center::Point{2},
    r_major::Real,
    r_minor::Real,
    angle::Real,
    a1::Real,
    a2::Real;
    stroke::Symbol = :black,
    pen::Symbol = :normal,
    arrow::Symbol = :none,
    rarrow::Symbol = :none,
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
    xml =
        "<path layer=\"$(canvas.active_layer)\"$attrs>\n" * @sprintf(
            "%.4f %.4f %.4f %.4f %.3f %.3f %.3f %.3f %.3f %.3f arc\n</path>",
            m11,
            m21,
            m12,
            m22,
            c[1],
            c[2],
            p1[1],
            p1[2],
            p2[1],
            p2[2]
        )
    push!(canvas.elements, xml)
    return canvas
end

draw_arc!(canvas::IpeCanvas, arc::EllipticArc; kwargs...) = draw_elliptic_arc!(
    canvas,
    arc.ellipse.center,
    arc.ellipse.r_major,
    arc.ellipse.r_minor,
    arc.ellipse.angle,
    arc.alpha1,
    arc.alpha2;
    kwargs...,
)

"""
    draw_bezier!(canvas, b::CubicBezier{2}; stroke=:black, pen=:normal, fill=:none, ...)

Draw a cubic Bézier curve.
"""
function draw_bezier!(
    canvas::IpeCanvas,
    b::CubicBezier{2};
    stroke::Symbol = :black,
    fill::Symbol = :none,
    pen::Symbol = :normal,
    dash::Symbol = :none,
    arrow::Symbol = :none,
    rarrow::Symbol = :none,
)
    p0 = _apply_tf(canvas, b.p0)
    p1 = _apply_tf(canvas, b.p1)
    p2 = _apply_tf(canvas, b.p2)
    p3 = _apply_tf(canvas, b.p3)
    attrs = "$(_fmt_attr("stroke", stroke))$(_fmt_attr("fill", fill))$(_fmt_attr("pen", pen))$(_fmt_attr("dash", dash))$(_fmt_attr("arrow", arrow))$(_fmt_attr("rarrow", rarrow))"
    xml =
        "<path layer=\"$(canvas.active_layer)\"$attrs>\n" *
        "$(_pt_str(p0)) m\n" *
        @sprintf(
            "%.3f %.3f %.3f %.3f %.3f %.3f c\n</path>",
            p1[1],
            p1[2],
            p2[1],
            p2[2],
            p3[1],
            p3[2]
        )
    push!(canvas.elements, xml)
    return canvas
end

"""
    draw_spline!(canvas, spline::CubicSpline{2}; stroke=:black, pen=:normal, fill=:none, ...)
    draw_spline!(canvas, pts::AbsPntSeq{2}; method=:catmull_rom, closed=false, ...)

Draw a composite cubic spline curve.
"""
function draw_spline!(
    canvas::IpeCanvas,
    s::CubicSpline{2};
    stroke::Symbol = :black,
    fill::Symbol = :none,
    pen::Symbol = :normal,
    dash::Symbol = :none,
    arrow::Symbol = :none,
    rarrow::Symbol = :none,
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
        push!(
            lines,
            @sprintf(
                "%.3f %.3f %.3f %.3f %.3f %.3f c",
                p1[1],
                p1[2],
                p2[1],
                p2[2],
                p3[1],
                p3[2]
            )
        )
    end
    s.is_closed && push!(lines, "h")
    push!(lines, "</path>")
    push!(canvas.elements, join(lines, "\n"))
    return canvas
end

draw_spline!(
    canvas::IpeCanvas,
    pts::AbsPntSeq{2};
    method::Symbol = :catmull_rom,
    closed::Bool = false,
    kwargs...,
) = draw_spline!(
    canvas,
    method == :natural ? interpolate_natural_spline(pts; closed = closed) :
    interpolate_catmull_rom(pts; closed = closed);
    kwargs...,
)

"""
    draw_bspline!(canvas, pts; closed=false, stroke=:black, pen=:normal, fill=:none, ...)

Draw an approximating cubic B-spline using native Ipe spline operators (`s` or `u`).
"""
function draw_bspline!(
    canvas::IpeCanvas,
    pts;
    closed::Bool = false,
    stroke::Symbol = :black,
    fill::Symbol = :none,
    pen::Symbol = :normal,
    dash::Symbol = :none,
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
    canvas::IpeCanvas,
    outer::AbsPntSeq{2},
    holes::AbsPntSeq{2}...;
    stroke::Symbol = :black,
    fill::Symbol = :gray7,
    pen::Symbol = :normal,
    dash::Symbol = :none,
    opacity::Union{Symbol,Nothing} = nothing,
    tiling::Symbol = :none,
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
    for i = 2:n
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
    f::Function,
    canvas::IpeCanvas;
    matrix::Union{Nothing,AbstractVector{<:Real}} = nothing,
    opacity::Union{Nothing,Symbol,String} = nothing,
    clip::Union{Nothing,AbstractString} = nothing,
)
    attrs = ""
    if matrix !== nothing
        @assert length(matrix) == 6 "Matrix must have 6 values [m11, m21, m12, m22, dx, dy]"
        attrs *= @sprintf(
            " matrix=\"%.4f %.4f %.4f %.4f %.3f %.3f\"",
            matrix[1],
            matrix[2],
            matrix[3],
            matrix[4],
            matrix[5],
            matrix[6]
        )
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
