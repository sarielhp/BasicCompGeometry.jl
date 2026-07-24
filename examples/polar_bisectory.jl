#! /bin/env julia

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


#
# examples/2d_random_ch.jl
#
# Generate 20 random points in the disk of radius 0.5 centered at the origin,
# compute their convex hull, and output a PDF with four pages:
#   Page 1: the random point set and origin
#   Page 2: the convex hull polygon, points, and origin
#   Page 3: the original convex hull polygon, origin, polar points & dual polygon
#           of the supporting lines, and the polar lines corresponding to the convex polygon vertices.
#   Page 4: all elements from Page 3 plus 100 polar points of the bisector lines B(t)
#           between the origin and 100 uniformly sampled points P(t) on the polygon.
#
# The output is written to output/2d_random_ch.pdf (1000x1000 pixels).

using BasicCompGeometry
using LinearAlgebra
using Cairo

"""
    supporting_lines(poly)

Compute the 2D supporting lines (hyperplanes) of each edge of a 2D convex polygon `poly`.
Returns a vector of `Plane` objects.
"""
function supporting_lines(poly)
    n = length(poly)
    lines = Plane{2,Float64}[]
    for i = 1:n
        p1 = poly[i]
        p2 = poly[mod1(i + 1, n)]
        # Normal vector orthogonal to segment p1 -> p2
        normal = point(p2[2] - p1[2], -(p2[1] - p1[1]))
        push!(lines, Plane(p1, normal))
    end
    return lines
end

"""
    closed_polygon_at(poly, t::Real)

Evaluate point on the closed perimeter of convex polygon `poly` for normalized parameter `t` in `[0, 1]`.
`t=0` and `t=1` map to the starting vertex `poly[1]`.
"""
function closed_polygon_at(poly, t::Real)
    pnts = Points(poly)
    closed_pnts = pnts[1] == pnts[end] ? pnts : vcat(pnts, [pnts[1]])

    n = length(closed_pnts)
    lens = zeros(n)
    for i = 2:n
        lens[i] = lens[i-1] + norm(closed_pnts[i] - closed_pnts[i-1])
    end
    total_len = lens[end]
    target = clamp(t, 0.0, 1.0) * total_len

    if target <= 0.0
        return closed_pnts[1]
    elseif target >= total_len
        return closed_pnts[end-1]
    end

    i = searchsortedfirst(lens, target)
    i <= 1 && return closed_pnts[1]

    seg_len = lens[i] - lens[i-1]
    local_t = (target - lens[i-1]) / seg_len
    return closed_pnts[i-1] * (1.0 - local_t) + closed_pnts[i] * local_t
end

"""
    bisector_polar_point(p::Point{2,Float64})

Compute the polar point of the bisector line B(t) between origin (0,0) and point p.
The bisector plane is dot(p, x) = dot(p, p)/2, so its polar point is 2*p / dot(p, p).
"""
function bisector_polar_point(p::Point{2,Float64})
    d = dot(p, p)
    return (2.0 / d) * p
end

"""
    draw_origin(cr, s, sz=10)

Draw a small filled square centered at the origin (0, 0).
`s` is the scale factor returned by `cairo_draw_setup`. `sz` is the side length in pixels.
"""
function draw_origin(cr, s::Real, sz::Real=10)
    half_w = (sz / 2) / s
    Cairo.new_path(cr)
    Cairo.rectangle(cr, -half_w, -half_w, 2 * half_w, 2 * half_w)
    Cairo.fill(cr)
end

"""
    draw_plane_line(cr, pl::Union{Plane{2,Float64}, Line{2,Float64}}, bb::BBox{2,Float64})

Clip a 2D plane or line to bounding box `bb` and draw the line segment across `bb`.
"""
function draw_plane_line(cr, shape::Union{Plane{2,Float64}, Line{2,Float64}}, bb::BBox{2,Float64})
    pl = shape isa Line ? Plane(shape) : shape
    nx, ny = pl.n[1], pl.n[2]
    d = dot(pl.n, pl.p)
    xmin, xmax = bb.mini[1], bb.maxi[1]
    ymin, ymax = bb.mini[2], bb.maxi[2]

    pts = Point{2,Float64}[]
    if abs(ny) > 1e-12
        y1 = (d - nx * xmin) / ny
        ymin - 1e-9 <= y1 <= ymax + 1e-9 && push!(pts, point(xmin, clamp(y1, ymin, ymax)))
        y2 = (d - nx * xmax) / ny
        ymin - 1e-9 <= y2 <= ymax + 1e-9 && push!(pts, point(xmax, clamp(y2, ymin, ymax)))
    end
    if abs(nx) > 1e-12
        x1 = (d - ny * ymin) / nx
        xmin - 1e-9 <= x1 <= xmax + 1e-9 && push!(pts, point(clamp(x1, xmin, xmax), ymin))
        x2 = (d - ny * ymax) / nx
        xmin - 1e-9 <= x2 <= xmax + 1e-9 && push!(pts, point(clamp(x2, xmin, xmax), ymax))
    end

    if length(pts) >= 2
        Cairo.new_path(cr)
        Cairo.move_to(cr, pts[1][1], pts[1][2])
        Cairo.line_to(cr, pts[2][1], pts[2][2])
        Cairo.stroke(cr)
    end
end

function generate_bisectory_pdf(filename::String, ps::PntSeq{2,Float64})
    # Compute convex hull
    hull = convex_hull(ps)

    # Compute supporting lines of the convex hull and their polar points
    slines = supporting_lines(hull)
    polar_pts = PntSeq([polar(l) for l in slines])
    dual_hull = convex_hull(polar_pts)

    # Compute polar lines corresponding to vertices of the convex polygon
    polar_lines = [polar(p) for p in hull]

    # Compute inverted circular arcs for bisector segments (halved hull edge segments)
    n_hull = length(hull)
    inv_bisector_arcs = Union{CircleArc2F, Segment2F}[]
    for i = 1:n_hull
        p1 = hull[i] / 2.0
        p2 = hull[mod1(i + 1, n_hull)] / 2.0
        seg_half = Segment(p1, p2)
        push!(inv_bisector_arcs, invert(seg_half))
    end

    # Compute supporting line inversions (full circles) for Page 5 background
    inv_circles = Union{Circle2F, Line{2,Float64}}[]
    for i = 1:n_hull
        p1 = hull[i] / 2.0
        p2 = hull[mod1(i + 1, n_hull)] / 2.0
        l_half = Line(p1, p2 - p1)
        push!(inv_circles, invert(l_half))
    end

    # Compute bounding boxes for page setup
    bb12 = BBox(ps) + 0.2
    bound!(bb12, point(0.0, 0.0))

    # Page 3 bounding box covers original points, hull, polar points, dual hull, unit circle, and origin
    bb3 = BBox(ps) + 0.2
    bound!(bb3, polar_pts)
    bound!(bb3, point(0.0, 0.0))
    bound!(bb3, BBox(point(-1.0, -1.0), point(1.0, 1.0)))
    bb3 = bb3 + 0.2

    # Page 4 bounding box includes everything from Page 3 plus bisector arcs
    bb4 = BBox(ps) + 0.2
    bound!(bb4, polar_pts)
    bound!(bb4, point(0.0, 0.0))
    bound!(bb4, BBox(point(-1.0, -1.0), point(1.0, 1.0)))
    for arc_obj in inv_bisector_arcs
        if arc_obj isa CircleArc2F
            p1 = arc_obj.center + arc_obj.radius * point(cos(arc_obj.theta1), sin(arc_obj.theta1))
            p2 = arc_obj.center + arc_obj.radius * point(cos(arc_obj.theta2), sin(arc_obj.theta2))
            bound!(bb4, p1)
            bound!(bb4, p2)
            bound!(bb4, arc_obj.center)
        elseif arc_obj isa Segment2F
            bound!(bb4, arc_obj.p)
            bound!(bb4, arc_obj.q)
        end
    end
    bb4 = bb4 + 0.2

    # Page 5 bounding box includes everything from Page 4 plus the full inverted supporting circles
    bb5 = BBox(ps) + 0.2
    bound!(bb5, polar_pts)
    bound!(bb5, point(0.0, 0.0))
    bound!(bb5, BBox(point(-1.0, -1.0), point(1.0, 1.0)))

    for inv_obj in inv_circles
        if inv_obj isa Circle2F
            c_min = inv_obj.center - point(inv_obj.radius, inv_obj.radius)
            c_max = inv_obj.center + point(inv_obj.radius, inv_obj.radius)
            bound!(bb5, BBox(c_min, c_max))
        end
    end
    bb5 = bb5 + 0.2

    # Create PDF surface
    cw, ch = 1000, 1000
    outdir = joinpath(@__DIR__, "..", "output")
    mkpath(outdir)
    pdf_path = joinpath(outdir, filename)
    surface = CairoPDFSurface(pdf_path, cw, ch)
    cr = CairoContext(surface)

    # Page 1: draw the point set and origin
    Cairo.save(cr)
    s = cairo_draw_setup(cr, bb12, cw, ch, 20)
    set_source_rgb(cr, 0.0, 0.0, 1.0)
    cairo_draw_points(cr, Points(ps), 5 / s)
    set_source_rgb(cr, 0.0, 0.6, 0.0) # green square for origin
    draw_origin(cr, s, 10)
    Cairo.restore(cr)
    Cairo.show_page(cr)

    # Page 2: draw the convex hull polygon, points, and origin
    Cairo.save(cr)
    s = cairo_draw_setup(cr, bb12, cw, ch, 20)
    set_source_rgb(cr, 0.0, 0.0, 1.0)
    cairo_draw_points(cr, Points(ps), 5 / s)
    set_source_rgb(cr, 1.0, 0.0, 0.0)
    set_line_width(cr, 2.0)
    cairo_draw_polygon(cr, hull, true)
    set_source_rgb(cr, 0.0, 0.6, 0.0)
    draw_origin(cr, s, 10)
    Cairo.restore(cr)
    Cairo.show_page(cr)

    # Page 3: draw original convex hull, origin, polar points, dual polygon, polar lines, and unit inversion circle
    Cairo.save(cr)
    s = cairo_draw_setup(cr, bb3, cw, ch, 20)

    # Draw polar lines of the convex hull vertices
    set_source_rgba(cr, 0.4, 0.4, 0.8, 0.6)
    set_line_width(cr, 1.0)
    for pl in polar_lines
        draw_plane_line(cr, pl, bb3)
    end

    # Draw unit inversion circle (dashed charcoal)
    set_source_rgba(cr, 0.2, 0.2, 0.2, 0.7)
    set_line_width(cr, 1.5)
    Cairo.set_dash(cr, [6.0, 4.0], 0.0)
    Cairo.new_path(cr)
    Cairo.arc(cr, 0.0, 0.0, 1.0, 0, 2pi)
    Cairo.stroke(cr)
    Cairo.set_dash(cr, Float64[], 0.0)

    # Draw original convex hull in black
    set_source_rgb(cr, 0.0, 0.0, 0.0)
    set_line_width(cr, 2.0)
    cairo_draw_polygon(cr, hull, true)

    # Draw origin in green square
    set_source_rgb(cr, 0.0, 0.6, 0.0)
    draw_origin(cr, s, 10)

    # Draw polar points of supporting lines in dark purple
    set_source_rgb(cr, 0.5, 0.0, 0.5)
    cairo_draw_points(cr, Points(polar_pts), 6 / s)

    # Draw dual convex hull polygon of polar points in green
    set_source_rgb(cr, 0.1, 0.75, 0.2)
    set_line_width(cr, 2.0)
    cairo_draw_polygon(cr, dual_hull, true)

    Cairo.restore(cr)
    Cairo.show_page(cr)

    # Color selection function: 2 alternating colors, plus a 3rd color if n_hull is odd for the last element
    function get_edge_color(i::Int, n::Int)
        colors = [
            (0.0, 0.75, 0.85),   # Color 1: Cyan
            (0.95, 0.50, 0.0),   # Color 2: Amber / Orange
            (0.85, 0.15, 0.70),  # Color 3: Magenta / Purple (for last element if n is odd)
        ]
        if isodd(n) && i == n
            return colors[3]
        else
            return colors[mod1(i, 2)]
        end
    end

    # Page 4: draw Page 3 elements + inverted bisector circular arcs with alternating colors & 30% sector fills
    Cairo.save(cr)
    s = cairo_draw_setup(cr, bb4, cw, ch, 20)

    # 1. Draw 30% transparent circular sector fills FIRST before anything else
    for i = 1:n_hull
        col = get_edge_color(i, n_hull)
        set_source_rgba(cr, col[1], col[2], col[3], 0.30)
        arc_obj = inv_bisector_arcs[i]
        if arc_obj isa CircleArc2F
            cx, cy = arc_obj.center[1], arc_obj.center[2]
            r = arc_obj.radius
            Cairo.new_path(cr)
            Cairo.move_to(cr, cx, cy)
            Cairo.arc(cr, cx, cy, r, arc_obj.theta1, arc_obj.theta2)
            Cairo.close_path(cr)
            Cairo.fill(cr)
        end
    end

    # 2. Draw unit inversion circle (dashed charcoal)
    set_source_rgba(cr, 0.2, 0.2, 0.2, 0.7)
    set_line_width(cr, 1.5)
    Cairo.set_dash(cr, [6.0, 4.0], 0.0)
    Cairo.new_path(cr)
    Cairo.arc(cr, 0.0, 0.0, 1.0, 0, 2pi)
    Cairo.stroke(cr)
    Cairo.set_dash(cr, Float64[], 0.0)

    # 3. Draw polar lines of the convex hull vertices
    set_source_rgba(cr, 0.4, 0.4, 0.8, 0.6)
    set_line_width(cr, 1.0)
    for pl in polar_lines
        draw_plane_line(cr, pl, bb4)
    end

    # 4. Draw edges of the polar polygon (dual_hull) in Green
    set_source_rgb(cr, 0.1, 0.75, 0.2)
    set_line_width(cr, 2.0)
    cairo_draw_polygon(cr, dual_hull, true)

    # 5. Draw original convex hull polygon in Black
    set_source_rgb(cr, 0.0, 0.0, 0.0)
    set_line_width(cr, 2.0)
    cairo_draw_polygon(cr, hull, true)

    # 6. Draw inner halved polygon (hull / 2.0) with alternating colors
    set_line_width(cr, 2.5)
    for i = 1:n_hull
        col = get_edge_color(i, n_hull)
        set_source_rgb(cr, col...)
        p1 = hull[i] / 2.0
        p2 = hull[mod1(i + 1, n_hull)] / 2.0
        Cairo.new_path(cr)
        Cairo.move_to(cr, p1[1], p1[2])
        Cairo.line_to(cr, p2[1], p2[2])
        Cairo.stroke(cr)
    end

    # 7. Draw origin in green square
    set_source_rgb(cr, 0.0, 0.6, 0.0)
    draw_origin(cr, s, 10)

    # 8. Draw polar polygon vertices in matching alternating colors
    for i = 1:n_hull
        col = get_edge_color(i, n_hull)
        set_source_rgb(cr, col...)
        cairo_draw_points(cr, [polar_pts[i]], 6 / s)
    end

    # 9. Draw inverted bisector segment circular arcs with matching alternating colors
    set_line_width(cr, 3.0)
    for i = 1:n_hull
        col = get_edge_color(i, n_hull)
        set_source_rgb(cr, col...)
        arc_obj = inv_bisector_arcs[i]
        if arc_obj isa CircleArc2F
            Cairo.new_path(cr)
            Cairo.arc(cr, arc_obj.center[1], arc_obj.center[2], arc_obj.radius, arc_obj.theta1, arc_obj.theta2)
            Cairo.stroke(cr)
        elseif arc_obj isa Segment2F
            Cairo.new_path(cr)
            Cairo.move_to(cr, arc_obj.p[1], arc_obj.p[2])
            Cairo.line_to(cr, arc_obj.q[1], arc_obj.q[2])
            Cairo.stroke(cr)
        end
    end

    Cairo.restore(cr)
    Cairo.show_page(cr)

    # Page 5: inverted circles of halved polygon edges + all Page 4 elements
    Cairo.save(cr)
    s = cairo_draw_setup(cr, bb5, cw, ch, 20)

    # 1. Draw 30% transparent circular sector fills FIRST before anything else
    for i = 1:n_hull
        col = get_edge_color(i, n_hull)
        set_source_rgba(cr, col[1], col[2], col[3], 0.30)
        arc_obj = inv_bisector_arcs[i]
        if arc_obj isa CircleArc2F
            cx, cy = arc_obj.center[1], arc_obj.center[2]
            r = arc_obj.radius
            Cairo.new_path(cr)
            Cairo.move_to(cr, cx, cy)
            Cairo.arc(cr, cx, cy, r, arc_obj.theta1, arc_obj.theta2)
            Cairo.close_path(cr)
            Cairo.fill(cr)
        end
    end

    # 2. Draw full inverted circles of halved edges
    set_source_rgba(cr, 0.2, 0.6, 0.9, 0.4)
    set_line_width(cr, 1.5)
    for inv_obj in inv_circles
        if inv_obj isa Circle2F
            Cairo.new_path(cr)
            Cairo.arc(cr, inv_obj.center[1], inv_obj.center[2], inv_obj.radius, 0, 2pi)
            Cairo.stroke(cr)
        elseif inv_obj isa Line{2,Float64}
            draw_plane_line(cr, inv_obj, bb5)
        end
    end

    # 3. Draw unit inversion circle (dashed charcoal)
    set_source_rgba(cr, 0.2, 0.2, 0.2, 0.7)
    set_line_width(cr, 1.5)
    Cairo.set_dash(cr, [6.0, 4.0], 0.0)
    Cairo.new_path(cr)
    Cairo.arc(cr, 0.0, 0.0, 1.0, 0, 2pi)
    Cairo.stroke(cr)
    Cairo.set_dash(cr, Float64[], 0.0)

    # 4. Draw polar lines of the convex hull vertices
    set_source_rgba(cr, 0.4, 0.4, 0.8, 0.6)
    set_line_width(cr, 1.0)
    for pl in polar_lines
        draw_plane_line(cr, pl, bb5)
    end

    # 5. Draw edges of the polar polygon (dual_hull) in Green
    set_source_rgb(cr, 0.1, 0.75, 0.2)
    set_line_width(cr, 2.0)
    cairo_draw_polygon(cr, dual_hull, true)

    # 6. Draw original convex hull polygon in Black
    set_source_rgb(cr, 0.0, 0.0, 0.0)
    set_line_width(cr, 2.0)
    cairo_draw_polygon(cr, hull, true)

    # 7. Draw inner halved polygon (hull / 2.0) with alternating colors
    set_line_width(cr, 2.5)
    for i = 1:n_hull
        col = get_edge_color(i, n_hull)
        set_source_rgb(cr, col...)
        p1 = hull[i] / 2.0
        p2 = hull[mod1(i + 1, n_hull)] / 2.0
        Cairo.new_path(cr)
        Cairo.move_to(cr, p1[1], p1[2])
        Cairo.line_to(cr, p2[1], p2[2])
        Cairo.stroke(cr)
    end

    # 8. Draw origin in green square
    set_source_rgb(cr, 0.0, 0.6, 0.0)
    draw_origin(cr, s, 10)

    # 9. Draw polar polygon vertices in matching alternating colors
    for i = 1:n_hull
        col = get_edge_color(i, n_hull)
        set_source_rgb(cr, col...)
        cairo_draw_points(cr, [polar_pts[i]], 6 / s)
    end

    # 10. Draw inverted bisector segment circular arcs with matching alternating colors
    set_line_width(cr, 3.0)
    for i = 1:n_hull
        col = get_edge_color(i, n_hull)
        set_source_rgb(cr, col...)
        arc_obj = inv_bisector_arcs[i]
        if arc_obj isa CircleArc2F
            Cairo.new_path(cr)
            Cairo.arc(cr, arc_obj.center[1], arc_obj.center[2], arc_obj.radius, arc_obj.theta1, arc_obj.theta2)
            Cairo.stroke(cr)
        elseif arc_obj isa Segment2F
            Cairo.new_path(cr)
            Cairo.move_to(cr, arc_obj.p[1], arc_obj.p[2])
            Cairo.line_to(cr, arc_obj.q[1], arc_obj.q[2])
            Cairo.stroke(cr)
        end
    end

    Cairo.restore(cr)
    Cairo.show_page(cr)

    # Page 6: shape formed by circular arcs (filled light gray) with inner polar green polygon (filled white)
    Cairo.save(cr)
    s = cairo_draw_setup(cr, bb5, cw, ch, 20)

    # 1. Outer shape formed by circular arcs: filled in light gray, arcs stroked in same uniform color
    Cairo.new_path(cr)
    for i = 1:n_hull
        arc_obj = inv_bisector_arcs[i]
        if arc_obj isa CircleArc2F
            Cairo.arc(cr, arc_obj.center[1], arc_obj.center[2], arc_obj.radius, arc_obj.theta1, arc_obj.theta2)
        elseif arc_obj isa Segment2F
            Cairo.line_to(cr, arc_obj.q[1], arc_obj.q[2])
        end
    end
    Cairo.close_path(cr)
    set_source_rgb(cr, 0.9, 0.9, 0.9) # Light gray fill
    Cairo.fill_preserve(cr)
    set_source_rgb(cr, 0.15, 0.15, 0.5) # Uniform arc outline color
    set_line_width(cr, 3.0)
    Cairo.stroke(cr)

    # 2. Inner polar polygon (dual_hull): filled in white, edges in green
    Cairo.new_path(cr)
    Cairo.move_to(cr, dual_hull[1][1], dual_hull[1][2])
    for i = 2:length(dual_hull)
        Cairo.line_to(cr, dual_hull[i][1], dual_hull[i][2])
    end
    Cairo.close_path(cr)
    set_source_rgb(cr, 1.0, 1.0, 1.0) # White fill
    Cairo.fill_preserve(cr)
    set_source_rgb(cr, 0.1, 0.75, 0.2) # Green stroke
    set_line_width(cr, 2.0)
    Cairo.stroke(cr)

    Cairo.restore(cr)
    Cairo.show_page(cr)

    Cairo.finish(surface)
    println("Output written to ", pdf_path)
end

function main()
    # 1. Random points example (bisectory_random.pdf)
    n = 20
    pts = Point{2,Float64}[]
    c_disk = point(0.1, 0.0)
    for _ = 1:n
        angle = 2pi * rand()
        radius = 0.5 * sqrt(rand()) # uniform distribution in disk
        p = c_disk + point(radius * cos(angle), radius * sin(angle))
        push!(pts, p * 1.6)
    end
    ps_random = PntSeq(pts)
    generate_bisectory_pdf("bisectory_random.pdf", ps_random)

    # 2. Square of sidelength 0.8 centered at origin (bisectory_square.pdf)
    s_half = 0.8 / 2.0
    sq_pts = [
        point(-s_half, -s_half),
        point( s_half, -s_half),
        point( s_half,  s_half),
        point(-s_half,  s_half),
    ]
    ps_square = PntSeq(sq_pts)
    generate_bisectory_pdf("bisectory_square.pdf", ps_square)

    # 3. Regular 20-gon on circle of radius 0.8 centered at origin (bisectory_20.pdf)
    r20_pts = [point(0.8 * cos(2pi * i / 20), 0.8 * sin(2pi * i / 20)) for i = 0:19]
    ps_20 = PntSeq(r20_pts)
    generate_bisectory_pdf("bisectory_20.pdf", ps_20)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
