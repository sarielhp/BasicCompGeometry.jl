"""
    cairo_draw_setup(cr, bb, canvas_width, canvas_height, margin)

Set up Cairo context `cr` with a scale and translate transform so that
the bounding box `bb` is mapped to the canvas of size `canvas_width` x `canvas_height`,
with `margin` pixels of padding on each side.
The y-axis is flipped so that the coordinate system matches the standard
mathematical orientation (y increases upward).
"""
function cairo_draw_setup end

"""
    cairo_draw_points(cr, points, radius)

Draw points as filled circles of the given `radius` using the Cairo context `cr`.
`points` is any iterable of `Point{2, Float64}`.
"""
function cairo_draw_points end

"""
    cairo_draw_polygon(cr, poly, close)

Draw the edges of a polygon (or point sequence) `poly` using the Cairo context `cr`.
If `close` is true, the polygon is closed by connecting the last vertex back to the first.
"""
function cairo_draw_polygon end

export cairo_draw_setup, cairo_draw_points, cairo_draw_polygon