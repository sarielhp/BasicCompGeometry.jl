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
    n = Point{2,T}(-line.u.y, line.u.x)
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
        ni = Point{2,T}(-lines[i].u.y, lines[i].u.x)
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
            dir = Point{2,T}(-nd.y, nd.x)
            pt0 = Point{2,T}(0.0, 0.0)
            if abs(nd.x) > 1e-12
                pt0 = Point{2,T}(d / nd.x, 0.0)
            elseif abs(nd.y) > 1e-12
                pt0 = Point{2,T}(0.0, d / nd.y)
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
    sort!(candidates, by=pt -> atan(pt.y - p.y, pt.x - p.x))
    unique_pts = Point{2,T}[]
    for pt in candidates
        if isempty(unique_pts) || dist(pt, unique_pts[end]) > 1e-8
            push!(unique_pts, pt)
        end
    end
    return parabolas, unique_pts
end

"""
    voronoi_curved_polygon(p::Point{2,T})

Compute the Voronoi cell of point `p` w.r.t. the unit square boundary as a `CurvePolygon2D`.
"""
function voronoi_curved_polygon(p::Point{2,T}) where {T}
    lines = square_boundary_lines()
    parabolas, verts = voronoi_cell(p)
    length(verts) < 2 && return CurvePolygon2D(ParabolaArc{T}[])
    if length(verts) == 2
        p1, p2 = verts[1], verts[2]
        active_k = Int[]
        for k in 1:4
            d1 = abs(dist(p1, p) - distance_to_line(p1, lines[k]))
            d2 = abs(dist(p2, p) - distance_to_line(p2, lines[k]))
            if d1 < 1e-6 && d2 < 1e-6
                push!(active_k, k)
            end
        end
        if length(active_k) == 2
            arc1 = ParabolaArc(parabolas[active_k[1]], p1, p2)
            arc2 = ParabolaArc(parabolas[active_k[2]], p2, p1)
            mid1 = at(arc1, 0.5)
            if turn_sign(p1, mid1, p) > 0
                arc1 = ParabolaArc(parabolas[active_k[2]], p1, p2)
                arc2 = ParabolaArc(parabolas[active_k[1]], p2, p1)
            end
            return CurvePolygon2D([arc1, arc2])
        end
    end
    m = length(verts)
    arcs = ParabolaArc{T}[]
    for idx in 1:m
        p_cur = verts[idx]
        p_next = verts[mod1(idx + 1, m)]
        mid = (p_cur + p_next) / 2
        best_k = 0
        best_diff = Inf
        for k in 1:4
            d = abs(dist(mid, p) - distance_to_line(mid, lines[k]))
            if d < best_diff
                best_diff = d; best_k = k
            end
        end
        push!(arcs, ParabolaArc(parabolas[best_k], p_cur, p_next))
    end
    return CurvePolygon2D(arcs)
end

"""
    point_in_curved_polygon(pt::Point{2,T}, p::Point{2,T})

Check if a point `pt` is inside the Voronoi cell of point `p` using `is_inside`.
"""
function point_in_curved_polygon(pt::Point{2,T}, p::Point{2,T}) where {T}
    cp = voronoi_curved_polygon(p)
    return is_inside(pt, cp)
end

function sample_cell_boundary(p::Point{2,T}, n_samples=300) where {T}
    cp = voronoi_curved_polygon(p)
    isempty(cp.curves) && return Point{2,T}[]
    m = length(cp.curves)
    pts = Point{2,T}[]
    samples_per_curve = max(2, n_samples ÷ m)
    for c in cp.curves
        for s in range(0.0, 1.0, length=samples_per_curve)
            push!(pts, at(c, s))
        end
    end
    return pts
end

function cell_area(p::Point{2,T}) where {T}
    boundary = sample_cell_boundary(p, 400)
    isempty(boundary) && return 0.0
    area = 0.0
    for i in 1:length(boundary)
        j = mod1(i + 1, length(boundary))
        area += boundary[i].x * boundary[j].y
        area -= boundary[j].x * boundary[i].y
    end
    return abs(area) / 2
end

function draw_cell_diagram(canvas, p::Point{2,Float64})
    lines = square_boundary_lines()
    parabolas, verts = voronoi_cell(p)
    bb = BBox(point(-0.2, -0.2), point(1.2, 1.2))
    cairo_draw_setup(canvas, bb, 800, 800, 40)
    Cairo.set_source_rgb(canvas, 1, 1, 1)
    Cairo.paint(canvas)
    
    cp = voronoi_curved_polygon(p)
    if !isempty(cp.curves)
        boundary = sample_cell_boundary(p, 300)
        if !isempty(boundary)
            Cairo.set_source_rgb(canvas, 0.8, 0.9, 1.0)
            Cairo.move_to(canvas, boundary[1].x, boundary[1].y)
            for pt in boundary[2:end]
                Cairo.line_to(canvas, pt.x, pt.y)
            end
            Cairo.close_path(canvas)
            Cairo.fill(canvas)
        end
    end
    
    Cairo.set_source_rgb(canvas, 0.65, 0.65, 0.65)
    Cairo.set_line_width(canvas, 0.003)
    for h in parabolas
        pts = Point2F[]
        for t in range(-3.0, 3.0, length=300)
            pt = at(h, t)
            if -0.2 <= pt.x <= 1.2 && -0.2 <= pt.y <= 1.2
                push!(pts, pt)
            end
        end
        if !isempty(pts)
            Cairo.move_to(canvas, pts[1].x, pts[1].y)
            for pt in pts[2:end]
                Cairo.line_to(canvas, pt.x, pt.y)
            end
            Cairo.stroke(canvas)
        end
    end
    
    Cairo.set_source_rgb(canvas, 0, 0, 0)
    Cairo.set_line_width(canvas, 0.008)
    Cairo.move_to(canvas, 0, 0); Cairo.line_to(canvas, 1, 0)
    Cairo.line_to(canvas, 1, 1); Cairo.line_to(canvas, 0, 1)
    Cairo.close_path(canvas); Cairo.stroke(canvas)
    
    if !isempty(cp.curves)
        boundary = sample_cell_boundary(p, 300)
        if !isempty(boundary)
            Cairo.set_source_rgb(canvas, 0.0, 0.0, 0.8)
            Cairo.set_line_width(canvas, 0.007)
            Cairo.move_to(canvas, boundary[1].x, boundary[1].y)
            for pt in boundary[2:end]
                Cairo.line_to(canvas, pt.x, pt.y)
            end
            Cairo.close_path(canvas); Cairo.stroke(canvas)
            
            for v in verts
                Cairo.arc(canvas, v.x, v.y, 0.0015, 0, 2pi)
                Cairo.fill(canvas)
            end
        end
    end
    
    Cairo.set_source_rgb(canvas, 1, 0, 0)
    Cairo.arc(canvas, p.x, p.y, 0.02, 0, 2pi)
    Cairo.fill(canvas)
    
    area_val = cell_area(p)
    description(canvas, @sprintf("Voronoi cell of point (%.2f, %.2f) — area = %.4f", p.x, p.y, area_val))
end

draw_page1(canvas, p=point(0.3, 0.4)) = draw_cell_diagram(canvas, p)

function draw_page2(canvas)
    Cairo.reset_transform(canvas)
    n = 40
    xs = range(0.01, 0.99, length=n)
    ys = range(0.01, 0.99, length=n)
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
                z = (vals[i, j] + vals[i+1, j] + vals[i+1, j+1] + vals[i, j+1]) / 4
                intensity = clamp(0.3 + 0.7 * (z / max(max_val, 1e-9)), 0.0, 1.0)
                isnan(intensity) && (intensity = 0.3)
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

function draw_page4(canvas, p=point(0.05, 0.05))
    lines = square_boundary_lines()
    parabolas, verts = voronoi_cell(p)
    
    # 1. Main View (Full Unit Square)
    bb = BBox(point(-0.15, -0.15), point(1.15, 1.15))
    cairo_draw_setup(canvas, bb, 800, 800, 40)
    Cairo.set_source_rgb(canvas, 1, 1, 1)
    Cairo.paint(canvas)
    
    cp = voronoi_curved_polygon(p)
    
    # Fill Voronoi cell
    if !isempty(cp.curves)
        boundary = sample_cell_boundary(p, 300)
        if !isempty(boundary)
            Cairo.set_source_rgb(canvas, 0.8, 0.9, 1.0)
            Cairo.move_to(canvas, boundary[1].x, boundary[1].y)
            for pt in boundary[2:end]
                Cairo.line_to(canvas, pt.x, pt.y)
            end
            Cairo.close_path(canvas)
            Cairo.fill(canvas)
        end
    end
    
    # Draw parabolas
    Cairo.set_source_rgb(canvas, 0.65, 0.65, 0.65)
    Cairo.set_line_width(canvas, 0.004)
    for h in parabolas[1:2:3] # Active parabolas 1 (x=0) and 3 (y=0)
        pts = Point2F[]
        for t in range(-3.0, 3.0, length=300)
            pt = at(h, t)
            if -0.15 <= pt.x <= 1.15 && -0.15 <= pt.y <= 1.15
                push!(pts, pt)
            end
        end
        if !isempty(pts)
            Cairo.move_to(canvas, pts[1].x, pts[1].y)
            for pt in pts[2:end]
                Cairo.line_to(canvas, pt.x, pt.y)
            end
            Cairo.stroke(canvas)
        end
    end
    
    # Draw unit square
    Cairo.set_source_rgb(canvas, 0, 0, 0)
    Cairo.set_line_width(canvas, 0.008)
    Cairo.move_to(canvas, 0, 0); Cairo.line_to(canvas, 1, 0)
    Cairo.line_to(canvas, 1, 1); Cairo.line_to(canvas, 0, 1)
    Cairo.close_path(canvas); Cairo.stroke(canvas)
    
    # Draw cell boundary
    if !isempty(cp.curves)
        boundary = sample_cell_boundary(p, 300)
        if !isempty(boundary)
            Cairo.set_source_rgb(canvas, 0.0, 0.0, 0.8)
            Cairo.set_line_width(canvas, 0.007)
            Cairo.move_to(canvas, boundary[1].x, boundary[1].y)
            for pt in boundary[2:end]
                Cairo.line_to(canvas, pt.x, pt.y)
            end
            Cairo.close_path(canvas); Cairo.stroke(canvas)
            
            for v in verts
                Cairo.set_source_rgb(canvas, 0.0, 0.0, 0.8)
                Cairo.arc(canvas, v.x, v.y, 0.006, 0, 2pi)
                Cairo.fill(canvas)
            end
        end
    end
    
    # Draw site point
    Cairo.set_source_rgb(canvas, 1, 0, 0)
    Cairo.arc(canvas, p.x, p.y, 0.012, 0, 2pi)
    Cairo.fill(canvas)
    
    # 2. Magnified Inset in Top-Right
    Cairo.reset_transform(canvas)
    ix, iy, iw, ih = 390.0, 70.0, 350.0, 350.0
    
    Cairo.set_source_rgb(canvas, 1.0, 1.0, 1.0)
    Cairo.rectangle(canvas, ix, iy, iw, ih)
    Cairo.fill(canvas)
    
    Cairo.set_source_rgb(canvas, 0.2, 0.3, 0.5)
    Cairo.set_line_width(canvas, 1.5)
    Cairo.rectangle(canvas, ix, iy, iw, ih)
    Cairo.stroke(canvas)
    
    w_min = -0.015
    w_max = 0.215
    w_span = w_max - w_min
    
    function map_inset(wx, wy)
        sx = ix + (wx - w_min) / w_span * iw
        sy = (iy + ih) - (wy - w_min) / w_span * ih
        return sx, sy
    end
    
    # Inset: Axes
    sx0, sy0 = map_inset(0.0, 0.0)
    sx1, sy1 = map_inset(0.20, 0.0)
    sx2, sy2 = map_inset(0.0, 0.20)
    Cairo.set_source_rgb(canvas, 0.0, 0.0, 0.0)
    Cairo.set_line_width(canvas, 2.0)
    Cairo.move_to(canvas, sx0, sy0); Cairo.line_to(canvas, sx1, sy0); Cairo.stroke(canvas)
    Cairo.move_to(canvas, sx0, sy0); Cairo.line_to(canvas, sx0, sy2); Cairo.stroke(canvas)
    
    # Inset: Symmetry line y = x
    sd1_x, sd1_y = map_inset(0.0, 0.0)
    sd2_x, sd2_y = map_inset(0.20, 0.20)
    Cairo.set_source_rgb(canvas, 0.7, 0.7, 0.7)
    Cairo.set_line_width(canvas, 1.0)
    Cairo.move_to(canvas, sd1_x, sd1_y); Cairo.line_to(canvas, sd2_x, sd2_y); Cairo.stroke(canvas)
    
    # Inset: Shaded cell
    boundary = sample_cell_boundary(p, 300)
    if !isempty(boundary)
        Cairo.set_source_rgb(canvas, 0.8, 0.9, 1.0)
        bx, by = map_inset(boundary[1].x, boundary[1].y)
        Cairo.move_to(canvas, bx, by)
        for pt in boundary[2:end]
            bx, by = map_inset(pt.x, pt.y)
            Cairo.line_to(canvas, bx, by)
        end
        Cairo.close_path(canvas); Cairo.fill(canvas)
        
        # Stroke boundary
        Cairo.set_source_rgb(canvas, 0.0, 0.0, 0.8)
        Cairo.set_line_width(canvas, 2.5)
        bx, by = map_inset(boundary[1].x, boundary[1].y)
        Cairo.move_to(canvas, bx, by)
        for pt in boundary[2:end]
            bx, by = map_inset(pt.x, pt.y)
            Cairo.line_to(canvas, bx, by)
        end
        Cairo.close_path(canvas); Cairo.stroke(canvas)
    end
    
    # Inset: Vertices
    for (i, v) in enumerate(verts)
        vx, vy = map_inset(v.x, v.y)
        Cairo.set_source_rgb(canvas, 0.0, 0.2, 0.8)
        Cairo.arc(canvas, vx, vy, 4.0, 0, 2pi)
        Cairo.fill(canvas)
        
        Cairo.select_font_face(canvas, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_BOLD)
        Cairo.set_font_size(canvas, 11)
        Cairo.move_to(canvas, vx + 7, vy - 4)
        Cairo.show_text(canvas, @sprintf("V%d (%.3f, %.3f)", i, v.x, v.y))
    end
    
    # Inset: Site p
    px, py = map_inset(p.x, p.y)
    Cairo.set_source_rgb(canvas, 0.9, 0.1, 0.1)
    Cairo.arc(canvas, px, py, 6.0, 0, 2pi)
    Cairo.fill(canvas)
    Cairo.select_font_face(canvas, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_BOLD)
    Cairo.set_source_rgb(canvas, 0.8, 0.0, 0.0)
    Cairo.set_font_size(canvas, 12)
    Cairo.move_to(canvas, px + 8, py + 12)
    Cairo.show_text(canvas, "p (0.05, 0.05)")
    
    # Inset: Title badge
    Cairo.select_font_face(canvas, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_BOLD)
    Cairo.set_source_rgb(canvas, 0.2, 0.3, 0.5)
    Cairo.set_font_size(canvas, 12)
    Cairo.move_to(canvas, ix + 12, iy + 22)
    Cairo.show_text(canvas, "Magnified Corner View [0, 0.20]²")
    
    area_val = cell_area(p)
    description(canvas, @sprintf("Voronoi cell of corner site (0.05, 0.05) — 2 parabolic arcs, area = %.5f", area_val))
end

function draw_page5(canvas)
    Cairo.reset_transform(canvas)
    cw, ch = 800.0, 800.0
    
    # Clean background
    Cairo.set_source_rgb(canvas, 0.98, 0.99, 1.0)
    Cairo.paint(canvas)
    
    # Top header bar
    Cairo.set_source_rgb(canvas, 0.12, 0.23, 0.45)
    Cairo.rectangle(canvas, 0, 0, cw, 70)
    Cairo.fill(canvas)
    
    Cairo.select_font_face(canvas, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_BOLD)
    Cairo.set_source_rgb(canvas, 1.0, 1.0, 1.0)
    Cairo.set_font_size(canvas, 20)
    Cairo.move_to(canvas, 35, 44)
    Cairo.show_text(canvas, "Page 4 Analysis: Corner Voronoi Cell at p = (0.05, 0.05)")
    
    # Body text setup
    left_margin = 40.0
    y = 110.0
    
    function draw_section_heading(title)
        Cairo.select_font_face(canvas, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_BOLD)
        Cairo.set_source_rgb(canvas, 0.12, 0.25, 0.55)
        Cairo.set_font_size(canvas, 14)
        Cairo.move_to(canvas, left_margin, y)
        Cairo.show_text(canvas, title)
        y += 24.0
    end
    
    function draw_body_text(lines; is_code=false, color=(0.15, 0.15, 0.18))
        Cairo.set_source_rgb(canvas, color...)
        if is_code
            Cairo.select_font_face(canvas, "Monospace", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
            Cairo.set_font_size(canvas, 11)
        else
            Cairo.select_font_face(canvas, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
            Cairo.set_font_size(canvas, 12)
        end
        for line in lines
            Cairo.move_to(canvas, left_margin + 12.0, y)
            Cairo.show_text(canvas, line)
            y += 18.0
        end
        y += 12.0
    end
    
    draw_section_heading("1. Active Boundaries in the Corner Regime")
    draw_body_text([
        "For interior sites like p = (0.3, 0.4), all four square edges participate in the cell boundary (4 arcs).",
        "When the site is placed near a corner at p = (0.05, 0.05), the distances to the left edge (x = 0) and",
        "bottom edge (y = 0) are only 0.05, whereas the distances to the right edge (x = 1) and top edge (y = 1)",
        "are 0.95. Every point in the Voronoi cell is much closer to p than to the distant edges.",
        "Consequently, the right and top boundary lines are completely inactive."
    ])
    
    draw_section_heading("2. Two-Parabola Boundary & Exact Equations")
    draw_body_text([
        "The cell boundary is formed exclusively by 2 parabolic bisectors:",
        "  • Left edge  (x = 0):  dist(q, p) = x  =>  x = 0.025 + 10*(y - 0.05)^2",
        "  • Bottom edge (y = 0): dist(q, p) = y  =>  y = 0.025 + 10*(x - 0.05)^2"
    ], is_code=true, color=(0.1, 0.2, 0.35))
    
    draw_section_heading("3. Two Intersection Vertices Along the Diagonal")
    draw_body_text([
        "Because px = py = 0.05, the two parabolas intersect symmetrically along the diagonal y = x.",
        "Equating x = y yields the quadratic equation 10*x^2 - 2*x + 0.05 = 0 (or 200*x^2 - 40*x + 1 = 0).",
        "This produces two distinct real intersection points:",
        "  • Lower vertex V1 = ( (2 - sqrt(2))/20, (2 - sqrt(2))/20 )  ≈  (0.02929, 0.02929)",
        "  • Upper vertex V2 = ( (2 + sqrt(2))/20, (2 + sqrt(2))/20 )  ≈  (0.17071, 0.17071)"
    ], is_code=true, color=(0.1, 0.2, 0.35))
    
    draw_section_heading("4. The Self-Enclosing Lens & Vanishing Area")
    draw_body_text([
        "The two parabolic arcs between V1 and V2 close directly upon each other, forming a 2-arc",
        "convex \"leaf\" / \"lens\" that completely encloses the site p = (0.05, 0.05).",
        "",
        "  • Perimeter (arc length): 0.4388 (compared to 1.6242 for p = (0.3, 0.4))",
        "  • Cell Area:              0.0094 (compared to 0.1984 for p = (0.3, 0.4) and max 0.2189 at center)",
        "",
        "As p -> (0, 0), the vertices collapse as O(||p||) and the cell area vanishes quadratically as O(||p||^2)."
    ])
    
    # Bottom info box
    Cairo.set_source_rgb(canvas, 0.90, 0.94, 0.98)
    Cairo.rectangle(canvas, 40, y + 5, 720, 50)
    Cairo.fill(canvas)
    Cairo.set_source_rgb(canvas, 0.2, 0.4, 0.7)
    Cairo.set_line_width(canvas, 1.0)
    Cairo.rectangle(canvas, 40, y + 5, 720, 50)
    Cairo.stroke(canvas)
    
    Cairo.select_font_face(canvas, "Sans", Cairo.FONT_SLANT_ITALIC, Cairo.FONT_WEIGHT_NORMAL)
    Cairo.set_source_rgb(canvas, 0.1, 0.2, 0.4)
    Cairo.set_font_size(canvas, 11)
    Cairo.move_to(canvas, 55, y + 34)
    Cairo.show_text(canvas, "Summary: Page 4 illustrates the transition from 4-arc interior cells to 2-arc corner lenses.")
    
    description(canvas, "Mathematical analysis and explanation of the corner Voronoi cell on Page 4")
end

draw_page6(canvas, p=point(0.5, 0.10)) = draw_cell_diagram(canvas, p)

function draw_edge_analysis(canvas, p::Point{2,Float64}, page_diag_num::Int)
    Cairo.reset_transform(canvas)
    cw, ch = 800.0, 800.0
    
    lines = square_boundary_lines()
    parabolas, verts = voronoi_cell(p)
    cp = voronoi_curved_polygon(p)
    area_val = cell_area(p)
    perim = geom_length(cp)
    
    # Clean background
    Cairo.set_source_rgb(canvas, 0.98, 0.99, 1.0)
    Cairo.paint(canvas)
    
    # Top header bar
    Cairo.set_source_rgb(canvas, 0.12, 0.23, 0.45)
    Cairo.rectangle(canvas, 0, 0, cw, 70)
    Cairo.fill(canvas)
    
    Cairo.select_font_face(canvas, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_BOLD)
    Cairo.set_source_rgb(canvas, 1.0, 1.0, 1.0)
    Cairo.set_font_size(canvas, 20)
    Cairo.move_to(canvas, 35, 44)
    Cairo.show_text(canvas, @sprintf("Page %d Analysis: Edge-Adjacent Voronoi Cell at p = (%.2f, %.2f)", page_diag_num, p.x, p.y))
    
    # Body text setup
    left_margin = 40.0
    y = 110.0
    
    function draw_section_heading(title)
        Cairo.select_font_face(canvas, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_BOLD)
        Cairo.set_source_rgb(canvas, 0.12, 0.25, 0.55)
        Cairo.set_font_size(canvas, 14)
        Cairo.move_to(canvas, left_margin, y)
        Cairo.show_text(canvas, title)
        y += 24.0
    end
    
    function draw_body_text(text_lines; is_code=false, color=(0.15, 0.15, 0.18))
        Cairo.set_source_rgb(canvas, color...)
        if is_code
            Cairo.select_font_face(canvas, "Monospace", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
            Cairo.set_font_size(canvas, 11)
        else
            Cairo.select_font_face(canvas, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
            Cairo.set_font_size(canvas, 12)
        end
        for line in text_lines
            Cairo.move_to(canvas, left_margin + 12.0, y)
            Cairo.show_text(canvas, line)
            y += 18.0
        end
        y += 12.0
    end
    
    d_bot = p.y
    d_top = 1.0 - p.y
    k_bot = 1.0 / (2 * p.y)
    k_top = 1.0 / (2 * (1.0 - p.y))
    v_bot = p.y / 2
    v_top = 1.0 - (1.0 - p.y) / 2
    
    draw_section_heading("1. Edge-Adjacent Regime & Boundary Participation")
    draw_body_text([
        @sprintf("When the site is positioned at p = (%.2f, %.2f), it lies near the bottom edge (y = 0, distance %.2f),", p.x, p.y, d_bot),
        @sprintf("equidistant to the left and right edges (distance %.2f), and at distance %.2f from the top edge.", p.x, d_top),
        "Because the side and top edges are sufficiently close relative to the site, ALL FOUR square",
        "boundaries contribute parabolic arcs to the Voronoi cell boundary."
    ])
    
    draw_section_heading("2. Four Parabolic Bisectors & Exact Equations")
    draw_body_text([
        "The cell is bounded by 4 distinct parabolic segments:",
        @sprintf("  • Bottom edge (y = 0): y = %.4f + %.4f*(x - 0.5)^2", v_bot, k_bot),
        @sprintf("  • Left edge   (x = 0): x = 0.2500 + (y - %.2f)^2", p.y),
        @sprintf("  • Right edge  (x = 1): x = 0.7500 - (y - %.2f)^2", p.y),
        @sprintf("  • Top edge    (y = 1): y = %.4f - %.4f*(x - 0.5)^2", v_top, k_top)
    ], is_code=true, color=(0.1, 0.2, 0.35))
    
    draw_section_heading("3. Four Vertices & Bilateral Symmetry")
    if length(verts) == 4
        draw_body_text([
            "Because p is centered horizontally at x = 0.5, the cell possesses exact mirror symmetry across x = 0.5.",
            "The 4 Voronoi vertices form two symmetric pairs:",
            @sprintf("  • Bottom-Right V1 = (%.4f, %.4f)  and  Bottom-Left V4 = (%.4f, %.4f)", verts[1].x, verts[1].y, verts[4].x, verts[4].y),
            @sprintf("  • Top-Right    V2 = (%.4f, %.4f)  and  Top-Left    V3 = (%.4f, %.4f)", verts[2].x, verts[2].y, verts[3].x, verts[3].y)
        ], is_code=true, color=(0.1, 0.2, 0.35))
    else
        draw_body_text([
            "The cell vertices are computed from pairwise parabolic bisector intersections."
        ])
    end
    
    draw_section_heading("4. The Asymmetric Dome / Teardrop Shape")
    draw_body_text([
        "The resulting Voronoi cell forms an asymmetric dome shape flattened along the bottom and",
        "extending upward towards the center of the square:",
        @sprintf("  • Bottom arc: sharply curved, dipping to apex y = %.4f (halfway between p and bottom edge)", v_bot),
        @sprintf("  • Top arc:    broad, shallow crest with apex at (0.50, %.4f)", v_top),
        "",
        @sprintf("  • Perimeter (arc length): %.4f", perim),
        @sprintf("  • Cell Area:              %.4f (approx %.1f%% of maximum center area 0.2189)", area_val, area_val / 0.2189 * 100)
    ])
    
    # Bottom info box
    Cairo.set_source_rgb(canvas, 0.90, 0.94, 0.98)
    Cairo.rectangle(canvas, 40, y + 5, 720, 50)
    Cairo.fill(canvas)
    Cairo.set_source_rgb(canvas, 0.2, 0.4, 0.7)
    Cairo.set_line_width(canvas, 1.0)
    Cairo.rectangle(canvas, 40, y + 5, 720, 50)
    Cairo.stroke(canvas)
    
    Cairo.select_font_face(canvas, "Sans", Cairo.FONT_SLANT_ITALIC, Cairo.FONT_WEIGHT_NORMAL)
    Cairo.set_source_rgb(canvas, 0.1, 0.2, 0.4)
    Cairo.set_font_size(canvas, 11)
    Cairo.move_to(canvas, 55, y + 34)
    Cairo.show_text(canvas, @sprintf("Summary: Page %d shows the 4-arc asymmetric dome geometry for site p = (%.2f, %.2f).", page_diag_num, p.x, p.y))
    
    description(canvas, @sprintf("Mathematical analysis and explanation of the edge-adjacent Voronoi cell on Page %d", page_diag_num))
end

draw_page7(canvas) = draw_edge_analysis(canvas, point(0.5, 0.10), 6)
draw_page8(canvas, p=point(0.5, 0.02)) = draw_cell_diagram(canvas, p)
draw_page9(canvas) = draw_edge_analysis(canvas, point(0.5, 0.02), 8)

function main()
    Random.seed!(42)
    p = point(0.3, 0.4)
    
    computed_area = cell_area(p)
    println("Computed area: ", @sprintf("%.6f", computed_area))
    
    output_dir = joinpath(@__DIR__, "..", "output")
    isdir(output_dir) || mkpath(output_dir)
    canvas_path = joinpath(output_dir, "sq_vs_point.pdf")
    open_canvas(canvas_path, 800, 800; title="Point vs Square Voronoi Cell") do canvas
        draw_page1(canvas, p)
        Cairo.show_page(canvas)
        draw_page2(canvas)
        Cairo.show_page(canvas)
        draw_page3(canvas)
        Cairo.show_page(canvas)
        draw_page4(canvas, point(0.05, 0.05))
        Cairo.show_page(canvas)
        draw_page5(canvas)
        Cairo.show_page(canvas)
        draw_page6(canvas, point(0.5, 0.10))
        Cairo.show_page(canvas)
        draw_page7(canvas)
        Cairo.show_page(canvas)
        draw_page8(canvas, point(0.5, 0.02))
        Cairo.show_page(canvas)
        draw_page9(canvas)
    end
    println("Output: ", relpath(get_file_path(canvas_path), normpath(joinpath(@__DIR__, ".."))))
end

main()