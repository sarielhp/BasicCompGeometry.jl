#!/usr/bin/env julia
# Bisector cell example:
# 10 random lines tangent to the unit square, their bisectors with the origin,
# and the cell in their arrangement containing the origin.

using BasicCompGeometry
using Cairo
using LinearAlgebra
using Printf
using Random

function tangent_lines_to_unit_square(n::Int)
    square = PntSeq([point(-0.5, -0.5), point(0.5, -0.5), point(0.5, 0.5), point(-0.5, 0.5)])
    lines = Line{2,Float64}[]
    for _ in 1:n
        alpha = rand() * 2 * pi
        line, _ = tangent_line(square, alpha)
        push!(lines, line)
    end
    return lines
end

function angle_bisectors(l1::Line{2,T}, l2::Line{2,T}) where {T}
    n1 = Point{2,T}(-l1.u[2], l1.u[1])
    n1 = n1 / norm(n1)
    n2 = Point{2,T}(-l2.u[2], l2.u[1])
    n2 = n2 / norm(n2)
    c1 = dot(l1.p, n1)
    c2 = dot(l2.p, n2)
    nd1 = n1 - n2
    nd2 = n1 + n2
    d1 = c1 - c2
    d2 = c1 + c2
    p1 = Point{2,T}(0.0, 0.0)
    p2 = Point{2,T}(0.0, 0.0)
    if abs(nd1[1]) > 1e-12
        p1 = Point{2,T}(d1 / nd1[1], 0.0)
    elseif abs(nd1[2]) > 1e-12
        p1 = Point{2,T}(0.0, d1 / nd1[2])
    end
    if abs(nd2[1]) > 1e-12
        p2 = Point{2,T}(d2 / nd2[1], 0.0)
    elseif abs(nd2[2]) > 1e-12
        p2 = Point{2,T}(0.0, d2 / nd2[2])
    end
    dir1 = Point{2,T}(-nd1[2], nd1[1])
    dir2 = Point{2,T}(-nd2[2], nd2[1])
    return Line(p1, dir1), Line(p2, dir2)
end

function cell_containing_origin(lines, origin=point(0.0, 0.0))
    n = length(lines)
    bisectors = [Hyperbola(line, origin) for line in lines]
    candidates = Tuple{Point2F,Int,Int}[]
    for i in 1:n
        for j in (i+1):n
            b1, b2 = angle_bisectors(lines[i], lines[j])
            for bis in (b1, b2)
                pts = intersect_line_curve(bis, bisectors[i])
                for (pt, _) in pts
                    d_origin = dist(pt, origin)
                    on_boundary = true
                    for k in 1:n
                        if distance(pt, lines[k]) < d_origin - 1e-6
                            on_boundary = false
                            break
                        end
                    end
                    if on_boundary && abs(distance(pt, lines[i]) - d_origin) < 1e-6
                        push!(candidates, (pt, i, j))
                    end
                end
            end
        end
    end
    isempty(candidates) && return CurvePolygon2D(Curve2D{Float64}[]), bisectors
    unique_pts = Point2F[]
    seen = Set{Int}()
    for (pt, _, _) in candidates
        h = hash(round.(pt, digits=6))
        if !(h in seen)
            push!(unique_pts, pt)
            push!(seen, h)
        end
    end
    sort!(unique_pts, by=p -> atan(p[2], p[1]))
    m = length(unique_pts)
    curves = Curve2D{Float64}[]
    for idx in 1:m
        p_cur = unique_pts[idx]
        p_next = unique_pts[mod1(idx + 1, m)]
        mid = (p_cur + p_next) / 2
        best_k = 0
        best_diff = Inf
        for k in 1:n
            d = abs(dist(mid, origin) - distance(mid, lines[k]))
            if d < best_diff
                best_diff = d
                best_k = k
            end
        end
        push!(curves, bisectors[best_k])
    end
    return CurvePolygon2D(curves), bisectors
end

function draw_bisector_cell(lines, origin=point(0.0, 0.0))
    cell, bisectors = cell_containing_origin(lines, origin)
    bb = BBox(point(-3.0, -3.0), point(3.0, 3.0))
    canvas_path = joinpath(@__DIR__, "bisector_cell.pdf")
    open_canvas(canvas_path, 800, 800; title="Bisector Cell") do canvas
        cairo_draw_setup(canvas, bb, 800, 800, 40)
        Cairo.set_source_rgb(canvas, 1, 1, 1)
        Cairo.paint(canvas)
        Cairo.set_source_rgb(canvas, 0.9, 0.9, 0.9)
        for x in -3:0.5:3
            Cairo.move_to(canvas, x, -3)
            Cairo.line_to(canvas, x, 3)
            Cairo.stroke(canvas)
            Cairo.move_to(canvas, -3, x)
            Cairo.line_to(canvas, 3, x)
            Cairo.stroke(canvas)
        end
        colors = [(0.2, 0.2, 0.8), (0.8, 0.2, 0.2), (0.2, 0.8, 0.2),
                  (0.8, 0.5, 0.1), (0.5, 0.2, 0.8), (0.1, 0.7, 0.7),
                  (0.7, 0.1, 0.5), (0.5, 0.5, 0.5), (0.9, 0.6, 0.2), (0.3, 0.6, 0.3)]
        for (idx, h) in enumerate(bisectors)
            c = colors[mod1(idx, length(colors))]
            Cairo.set_source_rgb(canvas, c[1], c[2], c[3])
            Cairo.set_line_width(canvas, 0.02)
            pts = Point2F[]
            for branch in (1, 2)
                for t in range(-3.0, 3.0, length=200)
                    p = at(h, Float64(t); branch=branch)
                    if all(-3 .<= p .<= 3)
                        push!(pts, p)
                    end
                end
            end
            if !isempty(pts)
                Cairo.move_to(canvas, pts[1][1], pts[1][2])
                for p in pts[2:end]
                    Cairo.line_to(canvas, p[1], p[2])
                end
                Cairo.stroke(canvas)
            end
        end
        Cairo.set_source_rgb(canvas, 0, 0, 0)
        Cairo.set_line_width(canvas, 0.03)
        n_curves = length(cell.curves)
        for i in 1:n_curves
            h = cell.curves[i]
            p_start = point_on(cell, Float64(i-1) / n_curves)
            p_end = point_on(cell, Float64(i) / n_curves)
            bv = max(b_coeff(h), 1e-12)
            t_start = asinh((p_start[1] * sin(rotation_angle(h)) - p_start[2] * cos(rotation_angle(h))) / bv)
            t_end = asinh((p_end[1] * sin(rotation_angle(h)) - p_end[2] * cos(rotation_angle(h))) / bv)
            branch = 1
            if dist(at(h, 0.0; branch=1), origin) > dist(at(h, 0.0; branch=2), origin)
                branch = 2
            end
            Cairo.move_to(canvas, p_start[1], p_start[2])
            for t in range(t_start, t_end, length=50)
                p = at(h, Float64(t); branch=branch)
                Cairo.line_to(canvas, p[1], p[2])
            end
            Cairo.line_to(canvas, p_end[1], p_end[2])
            Cairo.stroke(canvas)
        end
        Cairo.set_source_rgb(canvas, 0, 0.5, 0)
        Cairo.set_line_width(canvas, 0.01)
        square = [point(-0.5, -0.5), point(0.5, -0.5), point(0.5, 0.5), point(-0.5, 0.5), point(-0.5, -0.5)]
        Cairo.move_to(canvas, square[1][1], square[1][2])
        for p in square[2:end]
            Cairo.line_to(canvas, p[1], p[2])
        end
        Cairo.stroke(canvas)
        Cairo.set_source_rgb(canvas, 1, 0, 0)
        Cairo.arc(canvas, origin[1], origin[2], 0.05, 0, 2pi)
        Cairo.fill(canvas)
        description(canvas, "Bisector cell of 10 random lines tangent to the unit square with the origin")
    end
    return get_file_path(canvas_path)
end

function main()
    Random.seed!(42)
    lines = tangent_lines_to_unit_square(10)
    origin = point(0.0, 0.0)
    path = draw_bisector_cell(lines, origin)
    println("Output: $path")
end

main()