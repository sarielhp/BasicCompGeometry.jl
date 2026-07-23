module CairoExt

using BasicCompGeometry
import Cairo

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

end # module