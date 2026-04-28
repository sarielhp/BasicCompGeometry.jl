#! /bin/env julia

"""
    Example: Longest Convex Subset

This example demonstrates how to find the largest subset of points that form 
a convex polygon from a set of random points in the plane.
"""

using BasicCompGeometry
using BasicCompGeometry.LongestConvexSubset
using Printf
using Cairo

"""
    draw_polygon(cr, P, f_close::Bool = false)

Draws the edges of a polygon `P` using the Cairo context `cr`.
"""
function draw_polygon(cr, P, f_close::Bool = false)
    nv::Int = length(P)
    (nv == 0) && return

    for i = 2:nv
        p = P[i-1]
        q = P[i]
        move_to(cr, p[1], p[2])
        line_to(cr, q[1], q[2])
    end

    if f_close
        q = P[1]
        line_to(cr, q[1], q[2])
    end

    Cairo.stroke(cr)
end

"""
    draw_points(cr, P)

Draws points as small circles using the Cairo context `cr`.
"""
function draw_points(cr, P)
    nv::Int = length(P)
    rad::Float64 = 0.005

    for i = 1:nv
        p = P[i]
        move_to(cr, p[1] - rad, p[2])
        arc(cr, p[1], p[2], rad, 0.0, 2pi)
    end

    Cairo.stroke(cr)
end

"""
    draw_polygon_vertices(cr, P, r::Float64)

Draws the vertices of polygon `P` as filled circles with radius `r`.
"""
function draw_polygon_vertices(cr, P, r::Float64)
    nv::Int = length(P)

    for i = 1:nv
        p = P[i]
        Cairo.arc(cr, p[1], p[2], r, 0.0, 2pi)
        Cairo.fill(cr)
    end

    Cairo.stroke(cr)
end

function square_delta(delta::Float64)
    q = Polygon2F()
    push!(
        q,
        point(-delta, -delta),
        point(-delta, 1.0 + delta),
        point(1.0 + delta, 1.0 + delta),
        point(1.0 + delta, -delta),
        point(-delta, -delta),
    )
    return q
end

"""
    output_to_pdf(ps, sol, filename)

Generates a PDF visualization of the points and the longest convex subset.
"""
function output_to_pdf(ps, sol, filename::String)
    # Bounding box setup
    bb = BBox2F()
    bound!(bb, ps)

    # Expand for margin
    # Using the project's BBox operations
    center = point((bb.min[1] + bb.max[1])/2, (bb.min[2] + bb.max[2])/2)
    half_w = (bb.max[1] - bb.min[1]) * 1.1 / 2
    half_h = (bb.max[2] - bb.min[2]) * 1.1 / 2
    bbo = BBox(center - point(half_w, half_h), center + point(half_w, half_h))

    iwidth = 800
    iheight = Int(round(iwidth * (bbo.max[2] - bbo.min[2]) / (bbo.max[1] - bbo.min[1])))

    c = Cairo.CairoPDFSurface(filename, iwidth, iheight)
    cr = CairoContext(c)

    # Transform setup
    xcal = iwidth / (bbo.max[1] - bbo.min[1])
    Cairo.scale(cr, xcal, xcal)
    Cairo.translate(cr, -bbo.min[1], -bbo.min[2])

    # Draw the solution polygon
    set_line_width(cr, 0.005)
    set_source_rgb(cr, 0.5, 0.0, 1.0)
    draw_polygon(cr, sol, true)

    # Draw the vertices of the solution
    set_source_rgb(cr, 0.5, 0.5, 1.0)
    draw_polygon_vertices(cr, sol, 0.01)

    # Draw all input points
    set_source_rgb(cr, 1.0, 0.0, 0.0)
    draw_points(cr, ps)

    Cairo.finish(c)
    println("Output written to: ", filename)
end

function main(n::Int)
    # Generate random points
    ps = rand_polygon(2, Float64, n)

    println("Computing largest convex subset for $n points...")
    t = @timed sol = compute_largest_convex_subset(ps)

    @printf("Time: %.3f seconds, Largest convex subset size: %d\n", t.time, length(sol))

    # Create output directory if it doesn't exist
    if !isdir("output")
        mkdir("output")
    end

    filename = joinpath("output", "longest_convex_subset_$n.pdf")
    output_to_pdf(ps, sol, filename)
end

# Run the example
if abspath(PROGRAM_FILE) == @__FILE__
    main(50)
    main(100)
end
