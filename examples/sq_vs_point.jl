#!/usr/bin/env julia
# Voronoi cell of a point w.r.t. the boundary lines of the unit square.

using QuickEnv

using BasicCompGeometry
using Cairo
using LinearAlgebra
using Printf
using Random

function square_boundary_lines()
    left   = Line(point(0.0, 0.0), point(0.0, 1.0))
    right  = Line(point(1.0, 0.0), point(0.0, 1.0))
    bottom = Line(point(0.0, 0.0), point(1.0, 0.0))
    top    = Line(point(0.0, 1.0), point(1.0, 0.0))
    return [left, right, bottom, top]
end

function distance_to_line(x::Point{2,T}, line::Line{2,T}) where {T}
    n = Point{2,T}(-line.u[2], line.u[1])
    n = n / norm(n)
    return abs(dot(x - line.p, n))
end

function voronoi_cell(p::Point{2,T}) where {T}
    lines = square_boundary_lines()
    n = 4
    parabolas = [Parabola(p, lines[i]) for i in 1:n]
    interior = point(T(0.5), T(0.5))
    normals = Vector{Point{2,T}}(undef, n)
    for i in 1:n
        ni = Point{2,T}(-lines[i].u[2], lines[i].u[1])
        ni = ni / norm(ni)
        if dot(interior - lines[i].p, ni) < 0
            ni = -ni
        end
        normals[i] = ni
    end
    candidates = Point{2,T}[]
    for i in 1:n
        for j in (i+1):n
            nd = normals[i] - normals[j]
            if norm(nd) < 1e-12; continue; end
            ci = dot(lines[i].p, normals[i])
            cj = dot(lines[j].p, normals[j])
            d = ci - cj
            dir = Point{2,T}(-nd[2], nd[1])
            pt0 = Point{2,T}(0.0, 0.0)
            if abs(nd[1]) > 1e-12
                pt0 = Point{2,T}(d / nd[1], 0.0)
            elseif abs(nd[2]) > 1e-12
                pt0 = Point{2,T}(0.0, d / nd[2])
            end
            bis = Line(pt0, dir)
            pts = intersect_line_curve(bis, parabolas[i])
            for (pt, _) in pts
                d_p = dist(pt, p)
                ok = true
                for kk in 1:n
                    sd = dot(pt - lines[kk].p, normals[kk])
                    if sd < d_p - 1e-8
                        ok = false; break
                    end
                end
                if ok && abs(dot(pt - lines[i].p, normals[i]) - d_p) < 1e-8
                    push!(candidates, pt)
                end
            end
        end
    end
    isempty(candidates) && return parabolas, Point{2,T}[]
    sort!(candidates, by=pt -> atan(pt[2] - p[2], pt[1] - p[1]))
    unique_pts = Point{2,T}[]
    for pt in candidates
        if isempty(unique_pts) || dist(pt, unique_pts[end]) > 1e-8
            push!(unique_pts, pt)
        end
    end
    return parabolas, unique_pts
end

function cell_area(p::Point{2,T}) where {T}
    _, verts = voronoi_cell(p)
    length(verts) < 3 && return 0.0
    m = length(verts)
    samples = Point{2,T}[]
    for i in 1:m
        push!(samples, verts[i])
        if i < m
            mid = (verts[i] + verts[i+1]) / 2
            push!(samples, mid)
        end
    end
    mid = (verts[end] + verts[1]) / 2
    push!(samples, mid)
    area = 0.0
    for i in 1:length(samples)
        j = mod1(i + 1, length(samples))
        area += samples[i][1] * samples[j][2]
        area -= samples[j][1] * samples[i][2]
    end
    return abs(area) / 2
end

function sample_cell_boundary(p::Point{2,T}, n_samples=200) where {T}
    _, verts = voronoi_cell(p)
    length(verts) < 3 && return Point{2,T}[]
    lines = square_boundary_lines()
    parabolas = [Parabola(p, lines[i]) for i in 1:4]
    m = length(verts)
    pts = Point{2,T}[]
    for i in 1:m
        p_cur = verts[i]
        p_next = verts[mod1(i + 1, m)]
        best_k = 0
        best_diff = Inf
        mid = (p_cur + p_next) / 2
        for k in 1:4
            d = abs(dist(mid, p) - distance_to_line(mid, lines[k]))
            if d < best_diff
                best_diff = d; best_k = k
            end
        end
        h = parabolas[best_k]
        t_cur = (p_cur[1] - vertex(h)[1]) * (-axis_direction(h)[2]) + (p_cur[2] - vertex(h)[2]) * axis_direction(h)[1]
        t_next = (p_next[1] - vertex(h)[1]) * (-axis_direction(h)[2]) + (p_next[2] - vertex(h)[2]) * axis_direction(h)[1]
        for s in range(0.0, 1.0, length=n_samples ÷ m)
            t = t_cur + s * (t_next - t_cur)
            push!(pts, at(h, t))
        end
    end
    return pts
end

function draw_page1(canvas, p)
    lines = square_boundary_lines()
    parabolas, verts = voronoi_cell(p)
    bb = BBox(point(-0.2, -0.2), point(1.2, 1.2))
    cairo_draw_setup(canvas, bb, 800, 800, 40)
    Cairo.set_source_rgb(canvas, 1, 1, 1)
    Cairo.paint(canvas)
    if length(verts) >= 3
        boundary = sample_cell_boundary(p, 300)
        Cairo.set_source_rgb(canvas, 0.8, 0.9, 1.0)
        Cairo.move_to(canvas, boundary[1][1], boundary[1][2])
        for pt in boundary[2:end]
            Cairo.line_to(canvas, pt[1], pt[2])
        end
        Cairo.close_path(canvas)
        Cairo.fill(canvas)
    end
    Cairo.set_source_rgb(canvas, 0.5, 0.5, 0.5)
    Cairo.set_line_width(canvas, 0.005)
    for h in parabolas
        pts = Point2F[]
        for t in range(-3.0, 3.0, length=200)
            pt = at(h, t)
            if -0.2 <= pt[1] <= 1.2 && -0.2 <= pt[2] <= 1.2
                push!(pts, pt)
            end
        end
        if !isempty(pts)
            Cairo.move_to(canvas, pts[1][1], pts[1][2])
            for pt in pts[2:end]
                Cairo.line_to(canvas, pt[1], pt[2])
            end
            Cairo.stroke(canvas)
        end
    end
    Cairo.set_source_rgb(canvas, 0, 0, 0)
    Cairo.set_line_width(canvas, 0.01)
    Cairo.move_to(canvas, 0, 0); Cairo.line_to(canvas, 1, 0)
    Cairo.line_to(canvas, 1, 1); Cairo.line_to(canvas, 0, 1)
    Cairo.close_path(canvas); Cairo.stroke(canvas)
    if length(verts) >= 3
        Cairo.set_source_rgb(canvas, 0, 0, 0.8)
        Cairo.set_line_width(canvas, 0.008)
        boundary = sample_cell_boundary(p, 300)
        Cairo.move_to(canvas, boundary[1][1], boundary[1][2])
        for pt in boundary[2:end]
            Cairo.line_to(canvas, pt[1], pt[2])
        end
        Cairo.close_path(canvas); Cairo.stroke(canvas)
        for v in verts
            Cairo.arc(canvas, v[1], v[2], 0.015, 0, 2pi)
            Cairo.fill(canvas)
        end
    end
    Cairo.set_source_rgb(canvas, 1, 0, 0)
    Cairo.arc(canvas, p[1], p[2], 0.02, 0, 2pi)
    Cairo.fill(canvas)
    area_val = cell_area(p)
    description(canvas, @sprintf("Voronoi cell of point (%.2f, %.2f) — area = %.4f", p[1], p[2], area_val))
end

function draw_page2(canvas)
    n = 40
    xs = range(0.0, 1.0, length=n)
    ys = range(0.0, 1.0, length=n)
    vals = zeros(n, n)
    max_val = 0.0
    for (i, x) in enumerate(xs)
        for (j, y) in enumerate(ys)
            a = cell_area(point(x, y))
            vals[i, j] = a
            if a > max_val; max_val = a; end
        end
    end
    margin = 60
    cw, ch = 800, 800
    plot_l = margin
    plot_r = cw - margin
    plot_b = ch - margin
    plot_t = margin
    pw = plot_r - plot_l
    ph = plot_b - plot_t
    Cairo.set_source_rgb(canvas, 1, 1, 1)
    Cairo.paint(canvas)
    azim = pi / 4
    elev = pi / 6
    scale_3d = 0.6
    function project(x, y, z)
        cx = (x - 0.5) * scale_3d
        cy = (y - 0.5) * scale_3d
        cz = z / max_val * 0.4
        rx = cx * cos(azim) - cy * sin(azim)
        ry = cx * sin(azim) * sin(elev) + cy * cos(azim) * sin(elev) + cz * cos(elev)
        px = plot_l + pw * (0.5 + rx)
        py = plot_b - ph * (0.5 + ry)
        return px, py
    end
    Cairo.set_source_rgb(canvas, 0.9, 0.9, 0.9)
    for i in 1:n
        for j in 1:n
            if i < n
                x1, y1 = project(xs[i], ys[j], vals[i, j])
                x2, y2 = project(xs[i+1], ys[j], vals[i+1, j])
                Cairo.move_to(canvas, x1, y1); Cairo.line_to(canvas, x2, y2); Cairo.stroke(canvas)
            end
            if j < n
                x1, y1 = project(xs[i], ys[j], vals[i, j])
                x2, y2 = project(xs[i], ys[j+1], vals[i, j+1])
                Cairo.move_to(canvas, x1, y1); Cairo.line_to(canvas, x2, y2); Cairo.stroke(canvas)
            end
        end
    end
    for i in 1:n
        for j in 1:n
            if i < n && j < n
                z = (vals[i, j] + vals[i+1, j] + vals[i, j+1] + vals[i+1, j+1]) / 4
                intensity = 0.3 + 0.7 * z / max_val
                x1, y1 = project(xs[i], ys[j], vals[i, j])
                x2, y2 = project(xs[i+1], ys[j], vals[i+1, j])
                x3, y3 = project(xs[i+1], ys[j+1], vals[i+1, j+1])
                x4, y4 = project(xs[i], ys[j+1], vals[i, j+1])
                Cairo.set_source_rgb(canvas, 0.2, 0.2 * intensity, 0.8 * intensity)
                Cairo.move_to(canvas, x1, y1)
                Cairo.line_to(canvas, x2, y2)
                Cairo.line_to(canvas, x3, y3)
                Cairo.line_to(canvas, x4, y4)
                Cairo.close_path(canvas); Cairo.fill(canvas)
            end
        end
    end
    Cairo.set_source_rgb(canvas, 0, 0, 0)
    Cairo.set_line_width(canvas, 1.5)
    px0, py0 = project(0.0, 0.0, 0.0)
    px1, py1 = project(1.0, 0.0, 0.0)
    px2, py2 = project(0.0, 1.0, 0.0)
    px3, py3 = project(0.0, 0.0, max_val)
    Cairo.move_to(canvas, px0, py0); Cairo.line_to(canvas, px1, py1); Cairo.stroke(canvas)
    Cairo.move_to(canvas, px0, py0); Cairo.line_to(canvas, px2, py2); Cairo.stroke(canvas)
    Cairo.move_to(canvas, px0, py0); Cairo.line_to(canvas, px3, py3); Cairo.stroke(canvas)
    Cairo.set_font_size(canvas, 12)
    Cairo.move_to(canvas, px1 + 5, py1 + 15); Cairo.show_text(canvas, "x")
    Cairo.move_to(canvas, px2 + 5, py2 + 15); Cairo.show_text(canvas, "y")
    Cairo.move_to(canvas, px3 + 5, py3 - 5); Cairo.show_text(canvas, "f(p)")
    description(canvas, @sprintf("Surface plot of f(p) = area of Voronoi cell. Max = %.4f", max_val))
end

function draw_page3(canvas)
    n = 100
    xs = range(0.0, 1.0, length=n)
    ys = range(0.0, 1.0, length=n)
    bb = BBox(point(-0.1, -0.1), point(1.1, 1.1))
    cairo_draw_setup(canvas, bb, 800, 800, 40)
    Cairo.set_source_rgb(canvas, 1, 1, 1)
    Cairo.paint(canvas)
    for (i, x) in enumerate(xs)
        for (j, y) in enumerate(ys)
            a = cell_area(point(x, y))
            if a > 1/6
                Cairo.set_source_rgb(canvas, 0.8, 0.9, 1.0)
                Cairo.rectangle(canvas, x - 0.005, y - 0.005, 0.01, 0.01)
                Cairo.fill(canvas)
            end
        end
    end
    Cairo.set_source_rgb(canvas, 0, 0, 0)
    Cairo.set_line_width(canvas, 0.01)
    Cairo.move_to(canvas, 0, 0); Cairo.line_to(canvas, 1, 0)
    Cairo.line_to(canvas, 1, 1); Cairo.line_to(canvas, 0, 1)
    Cairo.close_path(canvas); Cairo.stroke(canvas)
    description(canvas, "Region where f(p) > 1/6 (shaded)")
end

function main()
    Random.seed!(42)
    p = point(0.3, 0.4)
    output_dir = joinpath(@__DIR__, "..", "output")
    isdir(output_dir) || mkpath(output_dir)
    canvas_path = joinpath(output_dir, "sq_vs_point.pdf")
    open_canvas(canvas_path, 800, 800; title="Point vs Square Voronoi Cell") do canvas
        draw_page1(canvas, p)
        Cairo.show_page(canvas)
        draw_page2(canvas)
        Cairo.show_page(canvas)
        draw_page3(canvas)
    end
    println("Output: $(get_file_path(canvas_path))")
end

main()