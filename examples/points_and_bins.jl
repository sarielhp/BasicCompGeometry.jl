#!/usr/bin/env julia

using QuickEnv

using BasicCompGeometry
using Cairo
using LinearAlgebra, Printf, Random

# ==============================================================================
# Halfplane clipping (from prophet.jl)
# ==============================================================================

function clip_polygon_halfplane(poly::Vector{Point2F}, M::Point2F, n::Point2F; eps::Float64=1e-11)
    isempty(poly) && return Point2F[]
    out = Point2F[]
    nv = length(poly)
    for i in 1:nv
        A = poly[i]
        B = poly[mod1(i + 1, nv)]
        sA = dot(A - M, n)
        sB = dot(B - M, n)
        inA = (sA <= eps)
        inB = (sB <= eps)
        if inB
            if !inA
                denom = dot(B - A, n)
                t = abs(denom) > 1e-14 ? clamp(dot(M - A, n) / denom, 0.0, 1.0) : 0.5
                push!(out, A + t * (B - A))
            end
            push!(out, B)
        else
            if inA
                denom = dot(B - A, n)
                t = abs(denom) > 1e-14 ? clamp(dot(M - A, n) / denom, 0.0, 1.0) : 0.5
                push!(out, A + t * (B - A))
            end
        end
    end
    return out
end

# ==============================================================================
# Voronoi cell computation
# ==============================================================================

function compute_voronoi_cell(pts::Vector{Point2F}, j::Int, outer_box::Vector{Point2F})
    p_j = pts[j]
    poly = copy(outer_box)
    order = sort(collect(1:length(pts)), by = k -> dist_sq(pts[k], p_j))
    for k in order
        k == j && continue
        p_k = pts[k]
        M = (p_j + p_k) / 2.0
        n = p_k - p_j
        poly = clip_polygon_halfplane(poly, M, n)
        isempty(poly) && break
    end
    return poly
end

# ==============================================================================
# Curved polygon (bisector cell) computation — numerical sampling
# ==============================================================================

function sample_curved_cell_boundary(p::Point2F, square::Vector{Point2F}, n_samples::Int=200)
    lines = Line{2,Float64}[]
    for i in 1:length(square)
        p1 = square[i]
        p2 = square[mod1(i + 1, length(square))]
        push!(lines, Line(p1, p2 - p1))
    end
    max_r = maximum(dist(p, v) for v in square) * 2.0
    boundary = Point2F[]
    for i in 0:(n_samples-1)
        theta = 2 * pi * i / n_samples
        d = point(cos(theta), sin(theta))
        lo, hi = 0.0, max_r
        for _ in 1:60
            mid = (lo + hi) / 2
            x = p + mid * d
            d_p = dist(x, p)
            min_d_line = minimum(distance(x, l) for l in lines)
            if d_p <= min_d_line + 1e-10
                lo = mid
            else
                hi = mid
            end
        end
        push!(boundary, p + lo * d)
    end
    return boundary
end

# ==============================================================================
# Drawing helpers
# ==============================================================================

function draw_curved_cell(cr, boundary::Vector{Point2F})
    n = length(boundary)
    n < 3 && return
    Cairo.new_path(cr)
    Cairo.move_to(cr, boundary[1][1], boundary[1][2])
    for i in 2:n
        Cairo.line_to(cr, boundary[i][1], boundary[i][2])
    end
    Cairo.close_path(cr)
end

function draw_page_number(cr, page_num, cw, ch, margin; show=true)
    show || return
    Cairo.save(cr)
    Cairo.set_source_rgb(cr, 0.4, 0.4, 0.4)
    Cairo.select_font_face(cr, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
    Cairo.set_font_size(cr, 10.0)
    Cairo.move_to(cr, margin, ch - margin * 0.3)
    Cairo.show_text(cr, string(page_num))
    Cairo.restore(cr)
end

function draw_grid(cr, grid_n, cell_size, lw=1.0)
    Cairo.set_source_rgb(cr, 0, 0, 0)
    Cairo.set_line_width(cr, lw)
    for i in 0:grid_n
        Cairo.move_to(cr, i * cell_size, 0)
        Cairo.line_to(cr, i * cell_size, 1)
        Cairo.stroke(cr)
        Cairo.move_to(cr, 0, i * cell_size)
        Cairo.line_to(cr, 1, i * cell_size)
        Cairo.stroke(cr)
    end
end

function fill_cell(cr, ix, iy, cell_size, r, g, b)
    cx0 = (ix - 1) * cell_size
    cy0 = (iy - 1) * cell_size
    Cairo.set_source_rgb(cr, r, g, b)
    Cairo.rectangle(cr, cx0, cy0, cell_size, cell_size)
    Cairo.fill(cr)
end

function draw_points(cr, pts, radius)
    for p in pts
        Cairo.arc(cr, p[1], p[2], radius, 0, 2pi)
        Cairo.fill(cr)
    end
end

# ==============================================================================
# Main
# ==============================================================================

function main(; show_page_numbers::Bool=true)
    N = 300
    grid_n = 10
    cell_size = 1.0 / grid_n

    Random.seed!(42)
    pts = [point(rand(), rand()) for _ in 1:N]

    cell_counts = zeros(Int, grid_n, grid_n)
    cell_points = [Point2F[] for _ in 1:grid_n, _ in 1:grid_n]
    for p in pts
        ix = min(floor(Int, p[1] * grid_n) + 1, grid_n)
        iy = min(floor(Int, p[2] * grid_n) + 1, grid_n)
        cell_counts[ix, iy] += 1
        push!(cell_points[ix, iy], p)
    end

    all_singletons = [(ix, iy) for ix in 1:grid_n, iy in 1:grid_n if cell_counts[ix, iy] == 1]
    k = length(all_singletons)
    @printf("Total singleton cells: %d\n", k)

    singleton_cell = argmin(all_singletons) do (ix, iy)
        (ix - 0.5 * (grid_n + 1))^2 + (iy - 0.5 * (grid_n + 1))^2
    end

    ix, iy = singleton_cell
    singleton_pt = cell_points[ix, iy][1]
    cx0, cy0 = (ix - 1) * cell_size, (iy - 1) * cell_size
    cx1, cy1 = ix * cell_size, iy * cell_size

    @printf("Center singleton cell: (%d, %d) = [%.4f, %.4f] x [%.4f, %.4f]\n", ix, iy, cx0, cx1, cy0, cy1)
    @printf("Singleton point: (%.6f, %.6f)\n", singleton_pt[1], singleton_pt[2])

    outer_box = [point(-0.1, -0.1), point(1.1, -0.1), point(1.1, 1.1), point(-0.1, 1.1)]
    pt_idx = findfirst(p -> p == singleton_pt, pts)
    voronoi_cell = compute_voronoi_cell(pts, pt_idx, outer_box)

    square = [point(cx0, cy0), point(cx1, cy0), point(cx1, cy1), point(cx0, cy1)]
    curved_boundary = sample_curved_cell_boundary(singleton_pt, square, 300)

    tiled_pts = copy(pts)
    for dy in -1:1, dx in -1:1
        (dx == 0 && dy == 0) && continue
        for p in pts
            push!(tiled_pts, p + point(dx, dy))
        end
    end
    singleton_indices = [findfirst(p -> p == cell_points[sx, sy][1], pts) for (sx, sy) in all_singletons]
    singleton_voronoi_cells = [compute_voronoi_cell(tiled_pts, idx, outer_box) for idx in singleton_indices]

    n_zoom = 20
    zoom_boxes = [BBox(
        point((1 - t) * 0.0 + t * cx0, (1 - t) * 0.0 + t * cy0),
        point((1 - t) * 1.0 + t * cx1, (1 - t) * 1.0 + t * cy1)
    ) for t in range(0.0, 1.0, length=n_zoom)]

    n_zoom_out = 20
    p15_x0 = max(0.0, (ix - 2) * cell_size)
    p15_x1 = min(1.0, (ix + 1) * cell_size)
    p15_y0 = max(0.0, (iy - 2) * cell_size)
    p15_y1 = min(1.0, (iy + 1) * cell_size)
    zoom_out_boxes = [BBox(
        point((1 - t) * p15_x0 + t * 0.0, (1 - t) * p15_y0 + t * 0.0),
        point((1 - t) * p15_x1 + t * 1.0, (1 - t) * p15_y1 + t * 1.0)
    ) for t in range(0.0, 1.0, length=n_zoom_out)]

    output_dir = joinpath(@__DIR__, "..", "output")
    mkpath(output_dir)
    pdf_path = joinpath(output_dir, "points_and_bins.pdf")
    cw, ch = 800, 800
    margin = 40.0

function draw_all_pages(canvas, pts, grid_n, cell_size, all_singletons, k, ix, iy,
            singleton_pt, cell_points, cx0, cx1, cy0, cy1, curved_boundary, voronoi_cell,
            singleton_voronoi_cells, n_zoom, zoom_boxes, n_zoom_out, zoom_out_boxes,
            cw, ch, margin; show_page_numbers=true)
        page_num = 1
        bg_r, bg_g, bg_b = 1.0, 0.995, 0.95

        # Page 1: Empty grid
        Cairo.save(canvas)
        Cairo.set_source_rgb(canvas, bg_r, bg_g, bg_b)
        Cairo.paint(canvas)
        cairo_draw_setup(canvas, BBox(point(0, 0), point(1, 1)), cw, ch, margin)
        draw_grid(canvas, grid_n, cell_size)
        Cairo.restore(canvas)
        draw_page_number(canvas, page_num, cw, ch, margin; show=show_page_numbers)
        description(canvas, "Page $page_num: Empty 10x10 grid")
        page_num += 1
        Cairo.show_page(canvas)

        # Page 2: Grid with points
        Cairo.save(canvas)
        Cairo.set_source_rgb(canvas, bg_r, bg_g, bg_b)
        Cairo.paint(canvas)
        cairo_draw_setup(canvas, BBox(point(0, 0), point(1, 1)), cw, ch, margin)
        fill_cell(canvas, ix, iy, cell_size, 0.85, 0.85, 0.85)
        draw_grid(canvas, grid_n, cell_size)
        Cairo.set_source_rgb(canvas, 0.85, 0.15, 0.15)
        draw_points(canvas, pts, 0.005)
        Cairo.set_source_rgb(canvas, 0.0, 0.3, 0.9)
        Cairo.arc(canvas, singleton_pt[1], singleton_pt[2], 0.008, 0, 2pi)
        Cairo.fill(canvas)
        Cairo.restore(canvas)
        draw_page_number(canvas, page_num, cw, ch, margin; show=show_page_numbers)
        description(canvas, "Page $page_num: $N random points in 10x10 grid (blue = singleton)")
        page_num += 1
        Cairo.show_page(canvas)

        # Pages 3-(3+k-1): Incrementally show singleton cells
        for (ci, (sx, sy)) in enumerate(all_singletons)
            Cairo.save(canvas)
            Cairo.set_source_rgb(canvas, bg_r, bg_g, bg_b)
            Cairo.paint(canvas)
            cairo_draw_setup(canvas, BBox(point(0, 0), point(1, 1)), cw, ch, margin)
            for (cj, (ssx, ssy)) in enumerate(all_singletons)
                if cj <= ci
                    fill_cell(canvas, ssx, ssy, cell_size, 0.85, 0.85, 0.85)
                end
            end
            draw_grid(canvas, grid_n, cell_size)
            Cairo.set_source_rgb(canvas, 0.85, 0.15, 0.15)
            draw_points(canvas, pts, 0.005)
            Cairo.set_source_rgb(canvas, 0.0, 0.3, 0.9)
            sp = cell_points[sx, sy][1]
            Cairo.arc(canvas, sp[1], sp[2], 0.008, 0, 2pi)
            Cairo.fill(canvas)
            Cairo.restore(canvas)
            draw_page_number(canvas, page_num, cw, ch, margin; show=show_page_numbers)
            description(canvas, "Page $page_num: Singleton cell $ci of $k")
            page_num += 1
            Cairo.show_page(canvas)
        end

        # Alternating green frames
        for fi in 1:6
            is_light = isodd(fi)
            r, g, b = is_light ? (0.8, 0.95, 0.8) : (0.3, 0.8, 0.3)
            Cairo.save(canvas)
            Cairo.set_source_rgb(canvas, bg_r, bg_g, bg_b)
            Cairo.paint(canvas)
            cairo_draw_setup(canvas, BBox(point(0, 0), point(1, 1)), cw, ch, margin)
            for (ssx, ssy) in all_singletons
                fill_cell(canvas, ssx, ssy, cell_size, r, g, b)
            end
            draw_grid(canvas, grid_n, cell_size)
            Cairo.set_source_rgb(canvas, 0.85, 0.15, 0.15)
            draw_points(canvas, pts, 0.005)
            Cairo.set_source_rgb(canvas, 0.0, 0.3, 0.9)
            Cairo.arc(canvas, singleton_pt[1], singleton_pt[2], 0.008, 0, 2pi)
            Cairo.fill(canvas)
            Cairo.restore(canvas)
            draw_page_number(canvas, page_num, cw, ch, margin; show=show_page_numbers)
            shade = is_light ? "light green" : "green"
            description(canvas, "Page $page_num: Singleton cells highlighted in $shade ($fi/6)")
            page_num += 1
            Cairo.show_page(canvas)
        end

        # Zoom into singleton cell
        for zi in 1:n_zoom
            bb = zoom_boxes[zi]
            bl = bottom_left(bb)
            tr = top_right(bb)
            bb_w = tr[1] - bl[1]
            max_r = bb_w / 200.0
            pt_r = min(0.005, max_r)
            sing_r = min(0.008, max_r * 1.6)
            Cairo.save(canvas)
            Cairo.set_source_rgb(canvas, bg_r, bg_g, bg_b)
            Cairo.paint(canvas)
            cairo_draw_setup(canvas, bb, cw, ch, margin)
            fill_cell(canvas, ix, iy, cell_size, 0.85, 0.85, 0.85)
            Cairo.rectangle(canvas, bl[1], bl[2], bb_w, tr[2] - bl[2])
            Cairo.clip(canvas.cr)
            draw_grid(canvas, grid_n, cell_size)
            Cairo.set_source_rgb(canvas, 0.85, 0.15, 0.15)
            draw_points(canvas, pts, pt_r)
            Cairo.set_source_rgb(canvas, 0.0, 0.3, 0.9)
            Cairo.arc(canvas, singleton_pt[1], singleton_pt[2], sing_r, 0, 2pi)
            Cairo.fill(canvas)
            Cairo.restore(canvas)
            draw_page_number(canvas, page_num, cw, ch, margin; show=show_page_numbers)
            description(canvas, "Page $page_num: Zoom into cell ($ix, $iy) — step $zi of $n_zoom")
            page_num += 1
            Cairo.show_page(canvas)
        end

        # Curved polygon
        p13_w = cx1 - cx0
        p13_max_r = p13_w / 200.0
        p13_sing_r = min(0.008, p13_max_r * 1.6)
        Cairo.save(canvas)
        Cairo.set_source_rgb(canvas, bg_r, bg_g, bg_b)
        Cairo.paint(canvas)
        cairo_draw_setup(canvas, BBox(point(cx0, cy0), point(cx1, cy1)), cw, ch, margin)
        fill_cell(canvas, ix, iy, cell_size, 0.85, 0.85, 0.85)
        draw_grid(canvas, grid_n, cell_size)
        Cairo.set_source_rgb(canvas, 0.7, 0.7, 0.7)
        Cairo.set_line_width(canvas, 0.002)
        Cairo.rectangle(canvas, cx0, cy0, cell_size, cell_size)
        Cairo.stroke(canvas)
        if !isempty(curved_boundary)
            Cairo.set_source_rgba(canvas, 0.2, 0.6, 0.2, 0.3)
            draw_curved_cell(canvas, curved_boundary)
            Cairo.fill(canvas)
            Cairo.set_source_rgb(canvas, 0.0, 0.5, 0.0)
            Cairo.set_line_width(canvas, 0.005)
            draw_curved_cell(canvas, curved_boundary)
            Cairo.stroke(canvas)
        end
        Cairo.set_source_rgb(canvas, 0.0, 0.3, 0.9)
        Cairo.arc(canvas, singleton_pt[1], singleton_pt[2], p13_sing_r, 0, 2pi)
        Cairo.fill(canvas)
        Cairo.restore(canvas)
        draw_page_number(canvas, page_num, cw, ch, margin; show=show_page_numbers)
        description(canvas, "Page $page_num: Curved polygon (bisector cell against grid cell boundary)")
        page_num += 1
        Cairo.show_page(canvas)

        # Zoom out to show neighbors
        nx0 = max(0.0, cx0 - cell_size)
        nx1 = min(1.0, cx1 + cell_size)
        ny0 = max(0.0, cy0 - cell_size)
        ny1 = min(1.0, cy1 + cell_size)
        p14_w = nx1 - nx0
        p14_max_r = p14_w / 200.0
        p14_pt_r = min(0.005, p14_max_r)
        p14_sing_r = min(0.008, p14_max_r * 1.6)
        Cairo.save(canvas)
        Cairo.set_source_rgb(canvas, bg_r, bg_g, bg_b)
        Cairo.paint(canvas)
        cairo_draw_setup(canvas, BBox(point(nx0, ny0), point(nx1, ny1)), cw, ch, margin)
        fill_cell(canvas, ix, iy, cell_size, 0.85, 0.85, 0.85)
        draw_grid(canvas, grid_n, cell_size)
        Cairo.set_source_rgb(canvas, 0.85, 0.15, 0.15)
        draw_points(canvas, pts, p14_pt_r)
        if !isempty(curved_boundary)
            Cairo.set_source_rgba(canvas, 0.2, 0.6, 0.2, 0.3)
            draw_curved_cell(canvas, curved_boundary)
            Cairo.fill(canvas)
            Cairo.set_source_rgb(canvas, 0.0, 0.5, 0.0)
            Cairo.set_line_width(canvas, 0.005)
            draw_curved_cell(canvas, curved_boundary)
            Cairo.stroke(canvas)
        end
        Cairo.set_source_rgb(canvas, 0.0, 0.3, 0.9)
        Cairo.arc(canvas, singleton_pt[1], singleton_pt[2], p14_sing_r, 0, 2pi)
        Cairo.fill(canvas)
        Cairo.restore(canvas)
        draw_page_number(canvas, page_num, cw, ch, margin; show=show_page_numbers)
        description(canvas, "Page $page_num: Grid cell and neighbors with curved polygon")
        page_num += 1
        Cairo.show_page(canvas)

        # Repeat last frame 4 times
        p15_x0 = max(0.0, (ix - 2) * cell_size)
        p15_x1 = min(1.0, (ix + 1) * cell_size)
        p15_y0 = max(0.0, (iy - 2) * cell_size)
        p15_y1 = min(1.0, (iy + 1) * cell_size)
        p15_w = p15_x1 - p15_x0
        p15_max_r = p15_w / 200.0
        p15_pt_r = min(0.003, p15_max_r)
        p15_sing_r = min(0.008, p15_max_r * 1.6)
        for rep in 1:4
            Cairo.save(canvas)
            Cairo.set_source_rgb(canvas, bg_r, bg_g, bg_b)
            Cairo.paint(canvas)
            cairo_draw_setup(canvas, BBox(point(p15_x0, p15_y0), point(p15_x1, p15_y1)), cw, ch, margin)
            fill_cell(canvas, ix, iy, cell_size, 0.85, 0.85, 0.85)
            draw_grid(canvas, grid_n, cell_size)
            Cairo.set_source_rgb(canvas, 0.85, 0.15, 0.15)
            draw_points(canvas, pts, p15_pt_r)
            if length(voronoi_cell) >= 3
                Cairo.set_source_rgba(canvas, 0.8, 0.3, 0.3, 0.25)
                Cairo.new_path(canvas)
                Cairo.move_to(canvas, voronoi_cell[1][1], voronoi_cell[1][2])
                for v in voronoi_cell[2:end]
                    Cairo.line_to(canvas, v[1], v[2])
                end
                Cairo.close_path(canvas)
                Cairo.fill(canvas)
                Cairo.set_source_rgb(canvas, 0.8, 0.1, 0.1)
                Cairo.set_line_width(canvas, 0.004)
                Cairo.new_path(canvas)
                Cairo.move_to(canvas, voronoi_cell[1][1], voronoi_cell[1][2])
                for v in voronoi_cell[2:end]
                    Cairo.line_to(canvas, v[1], v[2])
                end
                Cairo.close_path(canvas)
                Cairo.stroke(canvas)
            end
            if !isempty(curved_boundary) && length(curved_boundary) >= 3
                Cairo.set_source_rgba(canvas, 0.2, 0.6, 0.2, 0.3)
                draw_curved_cell(canvas, curved_boundary)
                Cairo.fill(canvas)
                Cairo.set_source_rgb(canvas, 0.0, 0.5, 0.0)
                Cairo.set_line_width(canvas, 0.005)
                draw_curved_cell(canvas, curved_boundary)
                Cairo.stroke(canvas)
            end
            Cairo.set_source_rgb(canvas, 0.0, 0.3, 0.9)
            Cairo.arc(canvas, singleton_pt[1], singleton_pt[2], p15_sing_r, 0, 2pi)
            Cairo.fill(canvas)
            Cairo.restore(canvas)
            draw_page_number(canvas, page_num, cw, ch, margin; show=show_page_numbers)
            description(canvas, "Page $page_num: Voronoi cell (red) and curved polygon (green) of singleton point — 3x3 neighborhood around cell ($ix, $iy)")
            page_num += 1
            Cairo.show_page(canvas)
        end

        # Zoom out to full grid
        for zi in 1:n_zoom_out
            bb = zoom_out_boxes[zi]
            bl = bottom_left(bb)
            tr = top_right(bb)
            bb_w = tr[1] - bl[1]
            max_r = bb_w / 200.0
            pt_r = min(0.005, max_r)
            sing_r = min(0.008, max_r * 1.6)
            Cairo.save(canvas)
            Cairo.set_source_rgb(canvas, bg_r, bg_g, bg_b)
            Cairo.paint(canvas)
            cairo_draw_setup(canvas, bb, cw, ch, margin)
            draw_grid(canvas, grid_n, cell_size)
            Cairo.set_source_rgb(canvas, 0.85, 0.15, 0.15)
            draw_points(canvas, pts, pt_r)
            for (ci, (ssx, ssy)) in enumerate(all_singletons)
                sp = cell_points[ssx, ssy][1]
                vc = singleton_voronoi_cells[ci]
                if length(vc) >= 3
                    Cairo.set_source_rgba(canvas, 0.2, 0.3, 0.8, 0.15)
                    Cairo.new_path(canvas)
                    Cairo.move_to(canvas, vc[1][1], vc[1][2])
                    for v in vc[2:end]
                        Cairo.line_to(canvas, v[1], v[2])
                    end
                    Cairo.close_path(canvas)
                    Cairo.fill(canvas)
                    Cairo.set_source_rgb(canvas, 0.0, 0.3, 0.9)
                    Cairo.set_line_width(canvas, 1.0)
                    Cairo.new_path(canvas)
                    Cairo.move_to(canvas, vc[1][1], vc[1][2])
                    for v in vc[2:end]
                        Cairo.line_to(canvas, v[1], v[2])
                    end
                    Cairo.close_path(canvas)
                    Cairo.stroke(canvas)
                end
                Cairo.set_source_rgb(canvas, 0.0, 0.3, 0.9)
                Cairo.arc(canvas, sp[1], sp[2], sing_r, 0, 2pi)
                Cairo.fill(canvas)
            end
            Cairo.restore(canvas)
            draw_page_number(canvas, page_num, cw, ch, margin; show=show_page_numbers)
            description(canvas, "Page $page_num: Zoom out — all singleton Voronoi cells (step $zi of $n_zoom_out)")
            page_num += 1
            Cairo.show_page(canvas)
        end

        n_last = page_num - 1
        Random.seed!(123)
        vc_colors = [(rand(), rand(), rand()) for _ in 1:length(all_singletons)]
        for extra in 1:n_last
            for ci in 1:length(all_singletons)
                if rand() < 1/3
                    vc_colors[ci] = (rand(), rand(), rand())
                end
            end
            Cairo.save(canvas)
            Cairo.set_source_rgb(canvas, bg_r, bg_g, bg_b)
            Cairo.paint(canvas)
            cairo_draw_setup(canvas, BBox(point(0, 0), point(1, 1)), cw, ch, margin)
            draw_grid(canvas, grid_n, cell_size)
            Cairo.set_source_rgb(canvas, 0.85, 0.15, 0.15)
            draw_points(canvas, pts, 0.005)
            for (ci, (ssx, ssy)) in enumerate(all_singletons)
                sp = cell_points[ssx, ssy][1]
                vc = singleton_voronoi_cells[ci]
                col = vc_colors[ci]
                if length(vc) >= 3
                    Cairo.set_source_rgb(canvas, col[1], col[2], col[3])
                    Cairo.new_path(canvas)
                    Cairo.move_to(canvas, vc[1][1], vc[1][2])
                    for v in vc[2:end]
                        Cairo.line_to(canvas, v[1], v[2])
                    end
                    Cairo.close_path(canvas)
                    Cairo.fill(canvas)
                    Cairo.set_source_rgb(canvas, 0.0, 0.3, 0.9)
                    Cairo.set_line_width(canvas, 1.0)
                    Cairo.new_path(canvas)
                    Cairo.move_to(canvas, vc[1][1], vc[1][2])
                    for v in vc[2:end]
                        Cairo.line_to(canvas, v[1], v[2])
                    end
                    Cairo.close_path(canvas)
                    Cairo.stroke(canvas)
                end
                Cairo.set_source_rgb(canvas, 0.0, 0.3, 0.9)
                Cairo.arc(canvas, sp[1], sp[2], 0.008, 0, 2pi)
                Cairo.fill(canvas)
            end
            Cairo.restore(canvas)
            draw_page_number(canvas, page_num, cw, ch, margin; show=show_page_numbers)
            description(canvas, "Page $page_num: Random Voronoi cell colors ($extra of $n_last)")
            page_num += 1
            Cairo.show_page(canvas)
        end
    end

    open_canvas(pdf_path, cw, ch) do canvas
        draw_all_pages(canvas, pts, grid_n, cell_size, all_singletons, k, ix, iy,
            singleton_pt, cell_points, cx0, cx1, cy0, cy1, curved_boundary, voronoi_cell,
            singleton_voronoi_cells, n_zoom, zoom_boxes, n_zoom_out, zoom_out_boxes,
            cw, ch, margin; show_page_numbers=show_page_numbers)
    end

    gif_path = joinpath(output_dir, "points_and_bins.gif")
    open_canvas(gif_path, cw, ch; fps=4) do canvas
        draw_all_pages(canvas, pts, grid_n, cell_size, all_singletons, k, ix, iy,
            singleton_pt, cell_points, cx0, cx1, cy0, cy1, curved_boundary, voronoi_cell,
            singleton_voronoi_cells, n_zoom, zoom_boxes, n_zoom_out, zoom_out_boxes,
            cw, ch, margin; show_page_numbers=show_page_numbers)
    end

    println("Output: ", relpath(get_file_path(pdf_path), normpath(joinpath(@__DIR__, ".."))))
    println("Output: ", relpath(get_file_path(gif_path), normpath(joinpath(@__DIR__, ".."))))

    # ==============================================================================
    # Key pages PDF + HTML (many_cells)
    # ==============================================================================
    function _render_many_cells!(canvas, pts, grid_n, cell_size, N, ix, iy, singleton_pt, cx0, cx1, cy0, cy1, curved_boundary, voronoi_cell, all_singletons, cell_points, singleton_voronoi_cells, cw, ch, margin; show_page_numbers=true)
        page_num = 1
        bg_r, bg_g, bg_b = 1.0, 0.995, 0.95

        # Page 1: Empty grid
        Cairo.save(canvas)
        Cairo.set_source_rgb(canvas, bg_r, bg_g, bg_b)
        Cairo.paint(canvas)
        cairo_draw_setup(canvas, BBox(point(0, 0), point(1, 1)), cw, ch, margin)
        draw_grid(canvas, grid_n, cell_size)
        Cairo.restore(canvas)
        draw_page_number(canvas, page_num, cw, ch, margin; show=show_page_numbers)
        description(canvas, "Page $page_num: Empty 10x10 grid")
        page_num += 1
        Cairo.show_page(canvas)

        # Page 2: Grid with points
        Cairo.save(canvas)
        Cairo.set_source_rgb(canvas, bg_r, bg_g, bg_b)
        Cairo.paint(canvas)
        cairo_draw_setup(canvas, BBox(point(0, 0), point(1, 1)), cw, ch, margin)
        fill_cell(canvas, ix, iy, cell_size, 0.85, 0.85, 0.85)
        draw_grid(canvas, grid_n, cell_size)
        Cairo.set_source_rgb(canvas, 0.85, 0.15, 0.15)
        draw_points(canvas, pts, 0.005)
        Cairo.set_source_rgb(canvas, 0.0, 0.3, 0.9)
        Cairo.arc(canvas, singleton_pt[1], singleton_pt[2], 0.008, 0, 2pi)
        Cairo.fill(canvas)
        Cairo.restore(canvas)
        draw_page_number(canvas, page_num, cw, ch, margin; show=show_page_numbers)
        description(canvas, "Page $page_num: $N random points in 10x10 grid (blue = singleton)")
        page_num += 1
        Cairo.show_page(canvas)

        # Curved polygon
        p13_w = cx1 - cx0
        p13_max_r = p13_w / 200.0
        p13_sing_r = min(0.008, p13_max_r * 1.6)
        Cairo.save(canvas)
        Cairo.set_source_rgb(canvas, bg_r, bg_g, bg_b)
        Cairo.paint(canvas)
        cairo_draw_setup(canvas, BBox(point(cx0, cy0), point(cx1, cy1)), cw, ch, margin)
        fill_cell(canvas, ix, iy, cell_size, 0.85, 0.85, 0.85)
        draw_grid(canvas, grid_n, cell_size)
        Cairo.set_source_rgb(canvas, 0.7, 0.7, 0.7)
        Cairo.set_line_width(canvas, 0.002)
        Cairo.rectangle(canvas, cx0, cy0, cell_size, cell_size)
        Cairo.stroke(canvas)
        if !isempty(curved_boundary)
            Cairo.set_source_rgba(canvas, 0.2, 0.6, 0.2, 0.3)
            draw_curved_cell(canvas, curved_boundary)
            Cairo.fill(canvas)
            Cairo.set_source_rgb(canvas, 0.0, 0.5, 0.0)
            Cairo.set_line_width(canvas, 0.005)
            draw_curved_cell(canvas, curved_boundary)
            Cairo.stroke(canvas)
        end
        Cairo.set_source_rgb(canvas, 0.0, 0.3, 0.9)
        Cairo.arc(canvas, singleton_pt[1], singleton_pt[2], p13_sing_r, 0, 2pi)
        Cairo.fill(canvas)
        Cairo.restore(canvas)
        draw_page_number(canvas, page_num, cw, ch, margin; show=show_page_numbers)
        description(canvas, "Page $page_num: Curved polygon (bisector cell against grid cell boundary)")
        page_num += 1
        Cairo.show_page(canvas)

        # Zoom out to show neighbors
        nx0 = max(0.0, cx0 - cell_size)
        nx1 = min(1.0, cx1 + cell_size)
        ny0 = max(0.0, cy0 - cell_size)
        ny1 = min(1.0, cy1 + cell_size)
        p14_w = nx1 - nx0
        p14_max_r = p14_w / 200.0
        p14_pt_r = min(0.005, p14_max_r)
        p14_sing_r = min(0.008, p14_max_r * 1.6)
        Cairo.save(canvas)
        Cairo.set_source_rgb(canvas, bg_r, bg_g, bg_b)
        Cairo.paint(canvas)
        cairo_draw_setup(canvas, BBox(point(nx0, ny0), point(nx1, ny1)), cw, ch, margin)
        fill_cell(canvas, ix, iy, cell_size, 0.85, 0.85, 0.85)
        draw_grid(canvas, grid_n, cell_size)
        Cairo.set_source_rgb(canvas, 0.85, 0.15, 0.15)
        draw_points(canvas, pts, p14_pt_r)
        if !isempty(curved_boundary)
            Cairo.set_source_rgba(canvas, 0.2, 0.6, 0.2, 0.3)
            draw_curved_cell(canvas, curved_boundary)
            Cairo.fill(canvas)
            Cairo.set_source_rgb(canvas, 0.0, 0.5, 0.0)
            Cairo.set_line_width(canvas, 0.005)
            draw_curved_cell(canvas, curved_boundary)
            Cairo.stroke(canvas)
        end
        Cairo.set_source_rgb(canvas, 0.0, 0.3, 0.9)
        Cairo.arc(canvas, singleton_pt[1], singleton_pt[2], p14_sing_r, 0, 2pi)
        Cairo.fill(canvas)
        Cairo.restore(canvas)
        draw_page_number(canvas, page_num, cw, ch, margin; show=show_page_numbers)
        description(canvas, "Page $page_num: Grid cell and neighbors with curved polygon")
        page_num += 1
        Cairo.show_page(canvas)

        # Voronoi cell + curved polygon (one frame)
        p15_x0 = max(0.0, (ix - 2) * cell_size)
        p15_x1 = min(1.0, (ix + 1) * cell_size)
        p15_y0 = max(0.0, (iy - 2) * cell_size)
        p15_y1 = min(1.0, (iy + 1) * cell_size)
        p15_w = p15_x1 - p15_x0
        p15_max_r = p15_w / 200.0
        p15_pt_r = min(0.003, p15_max_r)
        p15_sing_r = min(0.008, p15_max_r * 1.6)
        Cairo.save(canvas)
        Cairo.set_source_rgb(canvas, bg_r, bg_g, bg_b)
        Cairo.paint(canvas)
        cairo_draw_setup(canvas, BBox(point(p15_x0, p15_y0), point(p15_x1, p15_y1)), cw, ch, margin)
        fill_cell(canvas, ix, iy, cell_size, 0.85, 0.85, 0.85)
        draw_grid(canvas, grid_n, cell_size)
        Cairo.set_source_rgb(canvas, 0.85, 0.15, 0.15)
        draw_points(canvas, pts, p15_pt_r)
        if length(voronoi_cell) >= 3
            Cairo.set_source_rgba(canvas, 0.8, 0.3, 0.3, 0.25)
            Cairo.new_path(canvas)
            Cairo.move_to(canvas, voronoi_cell[1][1], voronoi_cell[1][2])
            for v in voronoi_cell[2:end]
                Cairo.line_to(canvas, v[1], v[2])
            end
            Cairo.close_path(canvas)
            Cairo.fill(canvas)
            Cairo.set_source_rgb(canvas, 0.8, 0.1, 0.1)
            Cairo.set_line_width(canvas, 0.004)
            Cairo.new_path(canvas)
            Cairo.move_to(canvas, voronoi_cell[1][1], voronoi_cell[1][2])
            for v in voronoi_cell[2:end]
                Cairo.line_to(canvas, v[1], v[2])
            end
            Cairo.close_path(canvas)
            Cairo.stroke(canvas)
        end
        if !isempty(curved_boundary) && length(curved_boundary) >= 3
            Cairo.set_source_rgba(canvas, 0.2, 0.6, 0.2, 0.3)
            draw_curved_cell(canvas, curved_boundary)
            Cairo.fill(canvas)
            Cairo.set_source_rgb(canvas, 0.0, 0.5, 0.0)
            Cairo.set_line_width(canvas, 0.005)
            draw_curved_cell(canvas, curved_boundary)
            Cairo.stroke(canvas)
        end
        Cairo.set_source_rgb(canvas, 0.0, 0.3, 0.9)
        Cairo.arc(canvas, singleton_pt[1], singleton_pt[2], p15_sing_r, 0, 2pi)
        Cairo.fill(canvas)
        Cairo.restore(canvas)
        draw_page_number(canvas, page_num, cw, ch, margin; show=show_page_numbers)
        description(canvas, "Page $page_num: Voronoi cell (red) and curved polygon (green) of singleton point — 3x3 neighborhood around cell ($ix, $iy)")
        page_num += 1
        Cairo.show_page(canvas)

        # Final frame: all singleton Voronoi cells with random colors
        Random.seed!(123)
        vc_colors = [(rand(), rand(), rand()) for _ in 1:length(all_singletons)]
        Cairo.save(canvas)
        Cairo.set_source_rgb(canvas, bg_r, bg_g, bg_b)
        Cairo.paint(canvas)
        cairo_draw_setup(canvas, BBox(point(0, 0), point(1, 1)), cw, ch, margin)
        draw_grid(canvas, grid_n, cell_size)
        Cairo.set_source_rgb(canvas, 0.85, 0.15, 0.15)
        draw_points(canvas, pts, 0.005)
        for (ci, (ssx, ssy)) in enumerate(all_singletons)
            sp = cell_points[ssx, ssy][1]
            vc = singleton_voronoi_cells[ci]
            col = vc_colors[ci]
            if length(vc) >= 3
                Cairo.set_source_rgb(canvas, col[1], col[2], col[3])
                Cairo.new_path(canvas)
                Cairo.move_to(canvas, vc[1][1], vc[1][2])
                for v in vc[2:end]
                    Cairo.line_to(canvas, v[1], v[2])
                end
                Cairo.close_path(canvas)
                Cairo.fill(canvas)
                Cairo.set_source_rgb(canvas, 0.0, 0.3, 0.9)
                Cairo.set_line_width(canvas, 1.0)
                Cairo.new_path(canvas)
                Cairo.move_to(canvas, vc[1][1], vc[1][2])
                for v in vc[2:end]
                    Cairo.line_to(canvas, v[1], v[2])
                end
                Cairo.close_path(canvas)
                Cairo.stroke(canvas)
            end
            Cairo.set_source_rgb(canvas, 0.0, 0.3, 0.9)
            Cairo.arc(canvas, sp[1], sp[2], 0.008, 0, 2pi)
            Cairo.fill(canvas)
        end
        Cairo.restore(canvas)
        draw_page_number(canvas, page_num, cw, ch, margin; show=show_page_numbers)
        description(canvas, "Page $page_num: All singleton Voronoi cells with random colors")
        page_num += 1
        Cairo.show_page(canvas)
    end

    many_cells_pdf = joinpath(output_dir, "many_cells.pdf")
    open_canvas(many_cells_pdf, cw, ch) do canvas
        _render_many_cells!(canvas, pts, grid_n, cell_size, N, ix, iy, singleton_pt, cx0, cx1, cy0, cy1, curved_boundary, voronoi_cell, all_singletons, cell_points, singleton_voronoi_cells, cw, ch, margin; show_page_numbers=show_page_numbers)
    end
    println("Output: ", relpath(get_file_path(many_cells_pdf), normpath(joinpath(@__DIR__, ".."))))

    many_cells_svg_dir = joinpath(output_dir, "many_cells_html")
    many_cells_html = joinpath(many_cells_svg_dir, "index.html")
    open_canvas(many_cells_html, cw, ch; title="Many Cells — Key Pages") do canvas
        _render_many_cells!(canvas, pts, grid_n, cell_size, N, ix, iy, singleton_pt, cx0, cx1, cy0, cy1, curved_boundary, voronoi_cell, all_singletons, cell_points, singleton_voronoi_cells, cw, ch, margin; show_page_numbers=show_page_numbers)
    end
    println("Output: ", relpath(get_file_path(many_cells_html), normpath(joinpath(@__DIR__, ".."))))
end

main(; show_page_numbers=true)