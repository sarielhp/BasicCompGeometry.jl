#! /bin/env julia

using Pkg
Pkg.activate(@__DIR__)

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
    draw_plane_line(cr, pl::Plane{2,Float64}, bb::BBox{2,Float64})

Clip a 2D plane (line) to bounding box `bb` and draw the line segment across `bb`.
"""
function draw_plane_line(cr, pl::Plane{2,Float64}, bb::BBox{2,Float64})
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

function main(n::Int=20)
    # Generate n random points inside circle of radius 0.5 centered at the origin
    pts = Point{2,Float64}[]
    for _ = 1:n
        angle = 2pi * rand()
        radius = 0.5 * sqrt(rand()) # uniform distribution in disk
        push!(pts, point(radius * cos(angle), radius * sin(angle)))
    end
    ps = PntSeq(pts)

    # Compute convex hull
    hull = convex_hull(ps)

    # Compute supporting lines of the convex hull and their polar points
    slines = supporting_lines(hull)
    polar_pts = PntSeq([polar(l) for l in slines])
    dual_hull = convex_hull(polar_pts)

    # Compute polar lines corresponding to vertices of the convex polygon
    polar_lines = [polar(p) for p in hull]
    # Sample 500 uniform parameter values t in [0, 1] along hull boundary
    ts = range(0.0, 1.0, length=500)
    bisector_polars = PntSeq([bisector_polar_point(closed_polygon_at(hull, t)) for t in ts])

    # Compute bounding boxes for page setup
    bb12 = BBox(ps) + 0.2
    bound!(bb12, point(0.0, 0.0))

    # Page 3 bounding box covers original points, hull, polar points, dual hull, and origin
    bb3 = BBox(ps) + 0.2
    bound!(bb3, polar_pts)
    bound!(bb3, point(0.0, 0.0))
    bb3 = bb3 + 0.2

    # Page 4 bounding box includes everything from Page 3 plus all 500 bisector polar points
    bb4 = BBox(ps) + 0.2
    bound!(bb4, polar_pts)
    bound!(bb4, bisector_polars)
    bound!(bb4, point(0.0, 0.0))
    bb4 = bb4 + 0.2

    # Create PDF surface
    cw, ch = 1000, 1000
    outdir = joinpath(@__DIR__, "..", "output")
    mkpath(outdir)
    surface = CairoPDFSurface(joinpath(outdir, "2d_random_ch.pdf"), cw, ch)
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

    # Page 3: draw original convex hull, origin, polar points, dual polygon, and polar lines
    Cairo.save(cr)
    s = cairo_draw_setup(cr, bb3, cw, ch, 20)

    # Draw polar lines of the convex hull vertices
    set_source_rgba(cr, 0.4, 0.4, 0.8, 0.6)
    set_line_width(cr, 1.0)
    for pl in polar_lines
        draw_plane_line(cr, pl, bb3)
    end

    # Draw original convex hull in red
    set_source_rgb(cr, 1.0, 0.0, 0.0)
    set_line_width(cr, 2.0)
    cairo_draw_polygon(cr, hull, true)

    # Draw origin in green square
    set_source_rgb(cr, 0.0, 0.6, 0.0)
    draw_origin(cr, s, 10)

    # Draw polar points of supporting lines in dark purple
    set_source_rgb(cr, 0.5, 0.0, 0.5)
    cairo_draw_points(cr, Points(polar_pts), 6 / s)

    # Draw dual convex hull polygon of polar points in orange
    set_source_rgb(cr, 0.9, 0.4, 0.0)
    set_line_width(cr, 2.0)
    cairo_draw_polygon(cr, dual_hull, true)

    Cairo.restore(cr)
    Cairo.show_page(cr)

    # Page 4: draw Page 3 elements + 100 bisector polar points
    Cairo.save(cr)
    s = cairo_draw_setup(cr, bb4, cw, ch, 20)

    # Draw polar lines of the convex hull vertices
    set_source_rgba(cr, 0.4, 0.4, 0.8, 0.6)
    set_line_width(cr, 1.0)
    for pl in polar_lines
        draw_plane_line(cr, pl, bb4)
    end

    # Draw original convex hull in red
    set_source_rgb(cr, 1.0, 0.0, 0.0)
    set_line_width(cr, 2.0)
    cairo_draw_polygon(cr, hull, true)

    # Draw origin in green square
    set_source_rgb(cr, 0.0, 0.6, 0.0)
    draw_origin(cr, s, 10)

    # Draw polar points of supporting lines in dark purple
    set_source_rgb(cr, 0.5, 0.0, 0.5)
    cairo_draw_points(cr, Points(polar_pts), 6 / s)

    # Draw dual convex hull polygon of polar points in orange
    set_source_rgb(cr, 0.9, 0.4, 0.0)
    set_line_width(cr, 2.0)
    cairo_draw_polygon(cr, dual_hull, true)

    # Draw 100 bisector polar points in teal/cyan
    set_source_rgb(cr, 0.0, 0.7, 0.8)
    cairo_draw_points(cr, Points(bisector_polars), 4 / s)

    Cairo.restore(cr)
    Cairo.show_page(cr)

    Cairo.finish(surface)
    println("Output written to ", joinpath(outdir, "2d_random_ch.pdf"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(20)
end
