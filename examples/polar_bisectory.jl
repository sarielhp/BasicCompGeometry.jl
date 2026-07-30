#! /bin/env julia

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


#
# examples/polar_bisectory.jl
#
# Compute convex hull, polar duality, bisector inversion arcs, curve polygon representation,
# tangent line intersection functions, periodic angle plots, and staircase ray shooting.
# Outputs multi-page PDFs (bisectory_random.pdf, bisectory_square.pdf, bisectory_20.pdf) with 10 pages:
#   Page 1: random point set and origin
#   Page 2: convex hull polygon, points, and origin
#   Page 3: convex hull polygon, origin, polar points, dual polygon, vertex polar lines, unit circle
#   Page 4: bisector line inversion arcs (halved hull edges) with 30% sector fills
#   Page 5: inverted full supporting circles of halved polygon edges + Page 4 elements
#   Page 6: outer shape formed by circular arcs as CurvePolygon2D CP with inner dual hull polygon Q
#   Page 7: tangent line to Q at random angle alpha, tangency point q, and intersections u (right) & v (left) on CP
#   Page 8: plots of before_tangent_to_polygon (u) and after_tangent_to_polygon (v) over alpha in [0, 2pi]
#   Page 9: periodic curve plot with vertical gray interval segments between u and v curves
#   Page 10: staircase algorithm S(x_1) with horizontal ray-shooting and self-intersection termination
#
# Generates 200 PNG frame images (output/frames/%06d.png) and an MP4 video (output/staircase.mp4).
#

using BasicCompGeometry
using LinearAlgebra
using Cairo
using Printf

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

    # Compute outer shape as CurvePolygon2D
    CP = CurvePolygon2D(inv_bisector_arcs)

    # 1. Outer shape formed by circular arcs: filled in light gray, arcs stroked in same uniform color
    Cairo.new_path(cr)
    for c in CP.curves
        if c isa CircleArc2F
            Cairo.arc(cr, c.center[1], c.center[2], c.radius, c.theta1, c.theta2)
        elseif c isa Segment2F
            Cairo.line_to(cr, c.q[1], c.q[2])
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

    # Page 7: Tangent line to Q at random alpha, tangency point q, and intersections u & v with CP
    Cairo.save(cr)
    s = cairo_draw_setup(cr, bb5, cw, ch, 20)

    # Alias inner polar polygon as Q
    Q = dual_hull
    alpha_rand = 2pi * rand()
    line_tan, q_tan = tangent_line(Q, alpha_rand)
    q_pt, u_pt, u_t, v_pt, v_t = BasicCompGeometry.tangent_intersections_cp(CP, Q, alpha_rand)

    # 1. Draw outer shape CP filled in light gray
    Cairo.new_path(cr)
    for c in CP.curves
        if c isa CircleArc2F
            Cairo.arc(cr, c.center[1], c.center[2], c.radius, c.theta1, c.theta2)
        elseif c isa Segment2F
            Cairo.line_to(cr, c.q[1], c.q[2])
        end
    end
    Cairo.close_path(cr)
    set_source_rgb(cr, 0.9, 0.9, 0.9)
    Cairo.fill_preserve(cr)
    set_source_rgb(cr, 0.15, 0.15, 0.5)
    set_line_width(cr, 3.0)
    Cairo.stroke(cr)

    # 2. Draw inner polar polygon Q (dual_hull) filled in white, green stroke
    Cairo.new_path(cr)
    Cairo.move_to(cr, Q[1][1], Q[1][2])
    for i = 2:length(Q)
        Cairo.line_to(cr, Q[i][1], Q[i][2])
    end
    Cairo.close_path(cr)
    set_source_rgb(cr, 1.0, 1.0, 1.0)
    Cairo.fill_preserve(cr)
    set_source_rgb(cr, 0.1, 0.75, 0.2)
    set_line_width(cr, 2.0)
    Cairo.stroke(cr)

    # 3. Draw origin in green square
    set_source_rgb(cr, 0.0, 0.6, 0.0)
    draw_origin(cr, s, 10)

    # 4. Draw tangent line (dashed magenta/purple)
    set_source_rgb(cr, 0.8, 0.1, 0.8)
    set_line_width(cr, 2.0)
    Cairo.set_dash(cr, [6.0, 4.0], 0.0)
    draw_plane_line(cr, line_tan, bb5)
    Cairo.set_dash(cr, Float64[], 0.0)

    # 5. Draw tangency point q (red circle)
    set_source_rgb(cr, 1.0, 0.0, 0.0)
    cairo_draw_points(cr, [q_pt], 7 / s)

    # 6. Draw u (right turn, blue circle) and v (left turn, orange circle)
    set_source_rgb(cr, 0.0, 0.3, 0.9) # u: blue
    cairo_draw_points(cr, [u_pt], 8 / s)

    set_source_rgb(cr, 0.95, 0.5, 0.0) # v: orange
    cairo_draw_points(cr, [v_pt], 8 / s)

    Cairo.restore(cr)
    Cairo.show_page(cr)

    # Page 8: Plot of before_tangent_to_polygon and after_tangent_to_polygon vs alpha ∈ [0, 2π]
    Cairo.save(cr)
    # Custom 2D plot coordinate mapping: [0, 2π] × [0, 1] mapped to page canvas with padding
    pad = 80.0
    plot_w = cw - 2 * pad
    plot_h = ch - 2 * pad

    function map_plot(a_val, y_val)
        x_px = pad + (a_val / 2pi) * plot_w
        y_px = ch - (pad + y_val * plot_h)
        return x_px, y_px
    end

    # Background fill (white)
    set_source_rgb(cr, 1.0, 1.0, 1.0)
    Cairo.rectangle(cr, 0, 0, cw, ch)
    Cairo.fill(cr)

    # Grid lines and axes
    set_source_rgb(cr, 0.9, 0.9, 0.9)
    set_line_width(cr, 1.0)
    for y_grid in 0.0:0.2:1.0
        x1, y1 = map_plot(0.0, y_grid)
        x2, y2 = map_plot(2pi, y_grid)
        Cairo.move_to(cr, x1, y1)
        Cairo.line_to(cr, x2, y2)
        Cairo.stroke(cr)
    end
    for a_grid in 0.0:(pi/2):2pi
        x1, y1 = map_plot(a_grid, 0.0)
        x2, y2 = map_plot(a_grid, 1.0)
        Cairo.move_to(cr, x1, y1)
        Cairo.line_to(cr, x2, y2)
        Cairo.stroke(cr)
    end

    # Axes border
    set_source_rgb(cr, 0.2, 0.2, 0.2)
    set_line_width(cr, 2.0)
    x_bl, y_bl = map_plot(0.0, 0.0)
    x_tr, y_tr = map_plot(2pi, 1.0)
    Cairo.rectangle(cr, pad, y_tr, plot_w, plot_h)
    Cairo.stroke(cr)

    # Text Labels on Axes
    Cairo.select_font_face(cr, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_BOLD)
    Cairo.set_font_size(cr, 18.0)

    # Y-axis labels
    for y_val in 0.0:0.25:1.0
        x_p, y_p = map_plot(0.0, y_val)
        Cairo.move_to(cr, x_p - 45, y_p + 6)
        Cairo.show_text(cr, string(y_val))
    end

    # X-axis labels
    x_labels = [("0", 0.0), ("π/2", pi/2), ("π", pi), ("3π/2", 3pi/2), ("2π", 2pi)]
    for (lbl, a_val) in x_labels
        x_p, y_p = map_plot(a_val, 0.0)
        Cairo.move_to(cr, x_p - 15, y_p + 30)
        Cairo.show_text(cr, lbl)
    end

    # Axis titles
    Cairo.move_to(cr, cw / 2 - 40, ch - pad + 60)
    Cairo.show_text(cr, "α (radians)")

    # Sample 800 values (4 times 200) on [0, 2π]
    n_samples = 800
    alphas = range(0.0, 2pi, length=n_samples)
    u_vals = [before_tangent_to_polygon(CP, Q, a) for a in alphas]
    v_vals = [after_tangent_to_polygon(CP, Q, a) for a in alphas]

    # Helper function to plot values while breaking path on wrap-around jump discontinuity (> 0.5)
    function draw_discontinuous_plot(cr, alphas, vals, color)
        set_source_rgb(cr, color...)
        set_line_width(cr, 2.5)
        Cairo.new_path(cr)
        x0, y0 = map_plot(alphas[1], vals[1])
        Cairo.move_to(cr, x0, y0)
        for i = 2:length(alphas)
            xi, yi = map_plot(alphas[i], vals[i])
            if abs(vals[i] - vals[i-1]) > 0.5
                Cairo.stroke(cr)
                Cairo.new_path(cr)
                Cairo.move_to(cr, xi, yi)
            else
                Cairo.line_to(cr, xi, yi)
            end
        end
        Cairo.stroke(cr)
    end

    # Plot u_vals (before_tangent_to_polygon - Blue)
    draw_discontinuous_plot(cr, alphas, u_vals, (0.0, 0.3, 0.9))

    # Plot v_vals (after_tangent_to_polygon - Orange)
    draw_discontinuous_plot(cr, alphas, v_vals, (0.95, 0.5, 0.0))

    Cairo.restore(cr)
    Cairo.show_page(cr)

    # Page 9: Same graph as Page 8 (without labels), plus vertical gray segments between before (u) and after (v) curves
    Cairo.save(cr)

    # Background fill (white)
    set_source_rgb(cr, 1.0, 1.0, 1.0)
    Cairo.rectangle(cr, 0, 0, cw, ch)
    Cairo.fill(cr)

    # Grid lines and axes
    set_source_rgb(cr, 0.9, 0.9, 0.9)
    set_line_width(cr, 1.0)
    for y_grid in 0.0:0.2:1.0
        x1, y1 = map_plot(0.0, y_grid)
        x2, y2 = map_plot(2pi, y_grid)
        Cairo.move_to(cr, x1, y1)
        Cairo.line_to(cr, x2, y2)
        Cairo.stroke(cr)
    end
    for a_grid in 0.0:(pi/2):2pi
        x1, y1 = map_plot(a_grid, 0.0)
        x2, y2 = map_plot(a_grid, 1.0)
        Cairo.move_to(cr, x1, y1)
        Cairo.line_to(cr, x2, y2)
        Cairo.stroke(cr)
    end

    # Axes border
    set_source_rgb(cr, 0.2, 0.2, 0.2)
    set_line_width(cr, 2.0)
    Cairo.rectangle(cr, pad, y_tr, plot_w, plot_h)
    Cairo.stroke(cr)

    # Text Labels on Axes
    Cairo.select_font_face(cr, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_BOLD)
    Cairo.set_font_size(cr, 18.0)

    for y_val in 0.0:0.25:1.0
        x_p, y_p = map_plot(0.0, y_val)
        Cairo.move_to(cr, x_p - 45, y_p + 6)
        Cairo.show_text(cr, string(y_val))
    end

    for (lbl, a_val) in x_labels
        x_p, y_p = map_plot(a_val, 0.0)
        Cairo.move_to(cr, x_p - 15, y_p + 30)
        Cairo.show_text(cr, lbl)
    end

    # Axis title
    Cairo.move_to(cr, cw / 2 - 40, ch - pad + 60)
    Cairo.show_text(cr, "α (radians)")

    # Draw vertical gray segments between u and v for each sampled alpha
    set_source_rgba(cr, 0.5, 0.5, 0.5, 0.6)
    set_line_width(cr, 1.2)
    for i = 1:n_samples
        a_i = alphas[i]
        u_i = u_vals[i]
        v_i = v_vals[i]

        if u_i <= v_i
            # Continuous interval from u_i up to v_i
            x_p, y_u = map_plot(a_i, u_i)
            _, y_v = map_plot(a_i, v_i)
            Cairo.new_path(cr)
            Cairo.move_to(cr, x_p, y_u)
            Cairo.line_to(cr, x_p, y_v)
            Cairo.stroke(cr)
        else
            # Wrap-around: segment from u_i to 1.0, and from 0.0 to v_i
            x_p, y_u = map_plot(a_i, u_i)
            _, y_top = map_plot(a_i, 1.0)
            Cairo.new_path(cr)
            Cairo.move_to(cr, x_p, y_u)
            Cairo.line_to(cr, x_p, y_top)
            Cairo.stroke(cr)

            _, y_bot = map_plot(a_i, 0.0)
            _, y_v = map_plot(a_i, v_i)
            Cairo.new_path(cr)
            Cairo.move_to(cr, x_p, y_bot)
            Cairo.line_to(cr, x_p, y_v)
            Cairo.stroke(cr)
        end
    end

    Cairo.restore(cr)
    Cairo.show_page(cr)

    # Draw Page 10 (PDF page) with random x1
    x1_pdf = 2pi * rand()
    draw_page_10!(cr, x1_pdf, cw, ch, pad, plot_w, plot_h, y_tr, alphas, u_vals, v_vals, n_samples, map_plot, draw_discontinuous_plot, x_labels)

    Cairo.restore(cr)
    Cairo.show_page(cr)

    Cairo.finish(surface)
    println("Output written to ", pdf_path)
end

function draw_page_10!(cr, x1::Float64, cw, ch, pad, plot_w, plot_h, y_tr, alphas, u_vals, v_vals, n_samples, map_plot, draw_discontinuous_plot, x_labels)
    # Background fill (white)
    set_source_rgb(cr, 1.0, 1.0, 1.0)
    Cairo.rectangle(cr, 0, 0, cw, ch)
    Cairo.fill(cr)

    # Grid lines and axes
    set_source_rgb(cr, 0.9, 0.9, 0.9)
    set_line_width(cr, 1.0)
    for y_grid in 0.0:0.2:1.0
        x1_p, y1_p = map_plot(0.0, y_grid)
        x2_p, y2_p = map_plot(2pi, y_grid)
        Cairo.move_to(cr, x1_p, y1_p)
        Cairo.line_to(cr, x2_p, y2_p)
        Cairo.stroke(cr)
    end
    for a_grid in 0.0:(pi/2):2pi
        x1_p, y1_p = map_plot(a_grid, 0.0)
        x2_p, y2_p = map_plot(a_grid, 1.0)
        Cairo.move_to(cr, x1_p, y1_p)
        Cairo.line_to(cr, x2_p, y2_p)
        Cairo.stroke(cr)
    end

    # Axes border
    set_source_rgb(cr, 0.2, 0.2, 0.2)
    set_line_width(cr, 2.0)
    Cairo.rectangle(cr, pad, y_tr, plot_w, plot_h)
    Cairo.stroke(cr)

    # Text Labels on Axes
    Cairo.select_font_face(cr, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_BOLD)
    Cairo.set_font_size(cr, 18.0)

    for y_val in 0.0:0.25:1.0
        x_p, y_p = map_plot(0.0, y_val)
        Cairo.move_to(cr, x_p - 45, y_p + 6)
        Cairo.show_text(cr, string(y_val))
    end

    for (lbl, a_val) in x_labels
        x_p, y_p = map_plot(a_val, 0.0)
        Cairo.move_to(cr, x_p - 15, y_p + 30)
        Cairo.show_text(cr, lbl)
    end

    # Axis title
    Cairo.move_to(cr, cw / 2 - 40, ch - pad + 60)
    Cairo.show_text(cr, "α (radians)")

    # Draw vertical gray segments between u and v for Page 10 background
    set_source_rgba(cr, 0.5, 0.5, 0.5, 0.6)
    set_line_width(cr, 1.2)
    for i = 1:n_samples
        a_i = alphas[i]
        u_i = u_vals[i]
        v_i = v_vals[i]

        if u_i <= v_i
            x_p, y_u = map_plot(a_i, u_i)
            _, y_v = map_plot(a_i, v_i)
            Cairo.new_path(cr)
            Cairo.move_to(cr, x_p, y_u)
            Cairo.line_to(cr, x_p, y_v)
            Cairo.stroke(cr)
        else
            x_p, y_u = map_plot(a_i, u_i)
            _, y_top = map_plot(a_i, 1.0)
            Cairo.new_path(cr)
            Cairo.move_to(cr, x_p, y_u)
            Cairo.line_to(cr, x_p, y_top)
            Cairo.stroke(cr)

            _, y_bot = map_plot(a_i, 0.0)
            _, y_v = map_plot(a_i, v_i)
            Cairo.new_path(cr)
            Cairo.move_to(cr, x_p, y_bot)
            Cairo.line_to(cr, x_p, y_v)
            Cairo.stroke(cr)
        end
    end

    # Plot u_vals (before - Blue)
    draw_discontinuous_plot(cr, alphas, u_vals, (0.0, 0.3, 0.9))

    # Plot v_vals (after - Orange)
    draw_discontinuous_plot(cr, alphas, v_vals, (0.95, 0.5, 0.0))

    # --- Staircase Algorithm ---
    function eval_fg(x_val)
        rel_x = mod(x_val, 2pi)
        idx = searchsortedfirst(alphas, rel_x)
        if idx <= 1
            return u_vals[1], v_vals[1]
        elseif idx > n_samples
            return u_vals[end], v_vals[end]
        end
        t_loc = (rel_x - alphas[idx-1]) / (alphas[idx] - alphas[idx-1])
        f_x = (1 - t_loc) * u_vals[idx-1] + t_loc * u_vals[idx]
        g_x = (1 - t_loc) * v_vals[idx-1] + t_loc * v_vals[idx]
        return f_x, g_x
    end

    function ray_shoot_f(x_curr, y_target)
        x_curr_mod = mod(x_curr, 2pi)
        idx_curr = searchsortedfirst(alphas, x_curr_mod)
        idx_curr = clamp(idx_curr, 1, n_samples)

        for k = 1:n_samples
            i1 = mod1(idx_curr + k - 1, n_samples)
            i2 = mod1(i1 + 1, n_samples)

            f1 = u_vals[i1]
            f2 = u_vals[i2]
            a1 = alphas[i1]
            a2 = i2 == 1 ? 2pi : alphas[i2]

            if (f1 <= y_target <= f2) || (f2 <= y_target <= f1)
                if f1 == f2
                    t_cross = 0.0
                else
                    t_cross = (y_target - f1) / (f2 - f1)
                end
                a_cross = a1 + t_cross * (a2 - a1)

                delta_a = mod(a_cross - x_curr_mod, 2pi)
                if k == 1 && delta_a < 1e-6
                    continue
                end
                return x_curr + delta_a
            end
        end
        return x_curr + 0.1
    end

    # Generate Staircase starting at x1
    stair_v_segments = Tuple{Float64,Float64,Float64}[]
    stair_h_segments = Tuple{Float64,Float64,Float64}[]

    f_x1, g_x1 = eval_fg(x1)
    push!(stair_v_segments, (x1, f_x1, g_x1))

    x_curr = x1
    y_curr_top = g_x1
    max_steps = 100

    for step = 1:max_steps
        x_next = ray_shoot_f(x_curr, y_curr_top)
        
        earliest_intersection_x = x_next
        found_intersection = false
        
        for (x_v, y_start_v, y_end_v) in stair_v_segments
            dx = mod(x_v - x_curr, 2pi)
            if dx < 1e-6
                continue
            end
            if dx <= x_next - x_curr + 1e-9
                in_interval = false
                if y_start_v <= y_end_v
                    in_interval = (y_start_v - 1e-9 <= y_curr_top <= y_end_v + 1e-9)
                else
                    in_interval = (y_curr_top >= y_start_v - 1e-9) || (y_curr_top <= y_end_v + 1e-9)
                end
                
                if in_interval
                    x_intersect = x_curr + dx
                    if x_intersect < earliest_intersection_x
                        earliest_intersection_x = x_intersect
                        found_intersection = true
                    end
                end
            end
        end

        push!(stair_h_segments, (y_curr_top, x_curr, earliest_intersection_x))

        if found_intersection
            break
        end

        f_next, g_next = eval_fg(earliest_intersection_x)
        push!(stair_v_segments, (earliest_intersection_x, f_next, g_next))

        x_curr = earliest_intersection_x
        y_curr_top = g_next
    end

    # Draw Staircase in Bold Crimson
    set_source_rgb(cr, 0.8, 0.0, 0.0)
    set_line_width(cr, 3.0)

    for (x_val, y_start, y_end) in stair_v_segments
        x_draw = mod(x_val, 2pi)
        if y_start <= y_end
            x_p, y1_p = map_plot(x_draw, y_start)
            _, y2_p = map_plot(x_draw, y_end)
            Cairo.new_path(cr)
            Cairo.move_to(cr, x_p, y1_p)
            Cairo.line_to(cr, x_p, y2_p)
            Cairo.stroke(cr)
        else
            x_p, y1_p = map_plot(x_draw, y_start)
            _, y_top = map_plot(x_draw, 1.0)
            Cairo.new_path(cr)
            Cairo.move_to(cr, x_p, y1_p)
            Cairo.line_to(cr, x_p, y_top)
            Cairo.stroke(cr)

            _, y_bot = map_plot(x_draw, 0.0)
            _, y2_p = map_plot(x_draw, y_end)
            Cairo.new_path(cr)
            Cairo.move_to(cr, x_p, y_bot)
            Cairo.line_to(cr, x_p, y2_p)
            Cairo.stroke(cr)
        end
    end

    for (y_val, x_start, x_end) in stair_h_segments
        k_start = floor(Int, x_start / (2pi))
        k_end = floor(Int, x_end / (2pi))

        if k_start == k_end
            a1 = mod(x_start, 2pi)
            a2 = mod(x_end, 2pi)
            if a2 == 0.0 && x_end > x_start
                a2 = 2pi
            end
            x1_p, y_p = map_plot(a1, y_val)
            x2_p, _ = map_plot(a2, y_val)
            Cairo.new_path(cr)
            Cairo.move_to(cr, x1_p, y_p)
            Cairo.line_to(cr, x2_p, y_p)
            Cairo.stroke(cr)
        else
            curr_x = x_start
            while curr_x < x_end
                next_boundary = (floor(curr_x / (2pi)) + 1) * 2pi
                seg_end = min(x_end, next_boundary)

                a1 = mod(curr_x, 2pi)
                a2 = seg_end == next_boundary ? 2pi : mod(seg_end, 2pi)

                x1_p, y_p = map_plot(a1, y_val)
                x2_p, _ = map_plot(a2, y_val)
                Cairo.new_path(cr)
                Cairo.move_to(cr, x1_p, y_p)
                Cairo.line_to(cr, x2_p, y_p)
                Cairo.stroke(cr)

                curr_x = seg_end
            end
        end
    end

    # Draw starting point (x1, f(x1)) as a red dot
    x1_draw = mod(x1, 2pi)
    x1_p, y1_p = map_plot(x1_draw, f_x1)
    set_source_rgb(cr, 1.0, 0.0, 0.0)
    Cairo.new_path(cr)
    Cairo.arc(cr, x1_p, y1_p, 6.0, 0, 2pi)
    Cairo.fill(cr)
end

"""
    generate_staircase_movie(ps::PntSeq{2,Float64}; num_frames=200, fps=10)

Generate 200 PNG frames of Page 10 for alpha1 in [0, 2π], saved in output/frames/%06d.png,
and render an MP4 movie output/staircase.mp4 using ImageMagick / convert at 10 fps.
"""
function generate_staircase_movie(ps::PntSeq{2,Float64}; num_frames=200, fps=10)
    # Pre-calculate convex hull, CP, Q, and functions before/after
    hull = convex_hull(ps)
    n_hull = length(hull)

    inv_bisector_arcs = Union{CircleArc2F, Segment2F}[]
    for i = 1:n_hull
        p1 = hull[i] / 2.0
        p2 = hull[mod1(i + 1, n_hull)] / 2.0
        seg_half = Segment(p1, p2)
        push!(inv_bisector_arcs, invert(seg_half))
    end
    CP = CurvePolygon2D(inv_bisector_arcs)

    slines = supporting_lines(hull)
    polar_pts = PntSeq([polar(l) for l in slines])
    dual_hull = convex_hull(polar_pts)
    Q = dual_hull

    # Pre-sample f and g with 800 values
    n_samples = 800
    alphas = range(0.0, 2pi, length=n_samples)
    u_vals = [before_tangent_to_polygon(CP, Q, a) for a in alphas]
    v_vals = [after_tangent_to_polygon(CP, Q, a) for a in alphas]

    cw, ch = 1000, 1000
    pad = 80.0
    plot_w = cw - 2 * pad
    plot_h = ch - 2 * pad

    function map_plot(a_val, y_val)
        x_px = pad + (a_val / 2pi) * plot_w
        y_px = ch - (pad + y_val * plot_h)
        return x_px, y_px
    end
    _, y_tr = map_plot(2pi, 1.0)

    function draw_discontinuous_plot(cr, alphas, vals, color)
        set_source_rgb(cr, color...)
        set_line_width(cr, 2.5)
        Cairo.new_path(cr)
        x0, y0 = map_plot(alphas[1], vals[1])
        Cairo.move_to(cr, x0, y0)
        for i = 2:length(alphas)
            xi, yi = map_plot(alphas[i], vals[i])
            if abs(vals[i] - vals[i-1]) > 0.5
                Cairo.stroke(cr)
                Cairo.new_path(cr)
                Cairo.move_to(cr, xi, yi)
            else
                Cairo.line_to(cr, xi, yi)
            end
        end
        Cairo.stroke(cr)
    end

    x_labels = [("0", 0.0), ("π/2", pi/2), ("π", pi), ("3π/2", 3pi/2), ("2π", 2pi)]

    outdir = joinpath(@__DIR__, "..", "output")
    frames_dir = joinpath(outdir, "frames")
    mkpath(frames_dir)

    alpha_list = range(0.0, 2pi, length=num_frames)

    println("Generating ", num_frames, " PNG frames in ", frames_dir, "...")
    for (idx, a1_val) in enumerate(alpha_list)
        img_surface = CairoARGBSurface(cw, ch)
        cr = CairoContext(img_surface)

        draw_page_10!(cr, a1_val, cw, ch, pad, plot_w, plot_h, y_tr, alphas, u_vals, v_vals, n_samples, map_plot, draw_discontinuous_plot, x_labels)

        frame_file = joinpath(frames_dir, @sprintf("%06d.png", idx))
        write_to_png(img_surface, frame_file)
        Cairo.finish(img_surface)
    end

    mp4_path = joinpath(outdir, "staircase.mp4")
    gif_path = joinpath(outdir, "staircase.gif")
    println("Assembling movie (fps=", fps, ")...")
    
    # Try ffmpeg for MP4 video assembly
    ffmpeg_cmd = `ffmpeg -y -framerate $fps -i $(joinpath(frames_dir, "%06d.png")) -c:v libx264 -pix_fmt yuv420p $mp4_path`
    try
        run(ffmpeg_cmd)
        println("MP4 movie generated successfully at: ", mp4_path)
    catch
        # Fallback to ImageMagick convert if ffmpeg is not installed
        try
            delay_val = round(Int, 100 / fps)
            frame_files = sort(readdir(frames_dir, join=true))
            run(`convert -delay $delay_val $frame_files $gif_path`)
            println("Animated GIF movie generated at: ", gif_path)
        catch e
            println("Could not assemble movie: ", e)
        end
    end
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

    # 4. Generate PNG frames and MP4 movie of Page 10 for random points
    generate_staircase_movie(ps_random; num_frames=200, fps=10)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
