#!/usr/bin/env julia

using Pkg
Pkg.activate(@__DIR__)

"""
    Example: Prophet Voronoi Diagram (Periodic 3x3 Grid)

Supported Modes:
  Mode 1: Prefix Voronoi sequence (starts from 4 points up to N points, default N = 100).
  Mode 2: Powers of 2 comparison (N = 2^3 ... 2^12, 10 pages) with ≤ 25% ink area calibration.
  Mode 3: (Default) 2-page Prophet comparison:
          - Page 1: Final Voronoi diagram of N = 1000 points, highlighting the prophet site p_k 
                    with the maximum cell area.
          - Page 2: Prefix Voronoi diagram of {p_1, ..., p_k} at the moment the prophet site p_k 
                    arrived, with its cell highlighted in the same blue style.

Usage:
  ./examples/prophet.jl             # Runs default Mode 3 (N = 1000, 2 pages)
  ./examples/prophet.jl 1 [N] [seed] # Mode 1: Prefix Voronoi sequence
  ./examples/prophet.jl 2 [seed]     # Mode 2: Powers of 2 (N = 2^3..2^12, 10 pages)
  ./examples/prophet.jl 3 [N] [seed] # Mode 3: 2-page Prophet comparison (default N = 1000)
"""

using BasicCompGeometry
using LinearAlgebra
using Cairo
using Printf
using Random

"""
    clip_polygon_halfplane(poly, M, n; eps=1e-11)

Clip convex polygon `poly` by the halfplane { x | dot(x - M, n) <= 0 } using the Sutherland-Hodgman algorithm.
`M` is a point on the dividing hyperplane and `n` is the outward normal vector pointing to the region to be clipped.
"""
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
                # Entering inside: compute intersection
                denom = dot(B - A, n)
                t = abs(denom) > 1e-14 ? clamp(dot(M - A, n) / denom, 0.0, 1.0) : 0.5
                push!(out, A + t * (B - A))
            end
            push!(out, B)
        else
            if inA
                # Leaving inside: compute intersection
                denom = dot(B - A, n)
                t = abs(denom) > 1e-14 ? clamp(dot(M - A, n) / denom, 0.0, 1.0) : 0.5
                push!(out, A + t * (B - A))
            end
        end
    end
    return out
end

"""
    create_3x3_replicated_points(pts)

Copy the points in the center square to the 8 surrounding squares in a 3x3 grid.
The first `length(pts)` points in the returned array are the original points in [0, 1]^2 (shift (0, 0)).
"""
function create_3x3_replicated_points(pts::Vector{Point2F})
    replicated = Point2F[]
    shifts = [
        point(0.0, 0.0),
        point(-1.0, 0.0), point(1.0, 0.0),
        point(0.0, -1.0), point(0.0, 1.0),
        point(-1.0, -1.0), point(1.0, -1.0),
        point(-1.0, 1.0), point(1.0, 1.0)
    ]
    for s in shifts
        for p in pts
            push!(replicated, p + s)
        end
    end
    return replicated
end

"""
    compute_voronoi_cells_periodic(pts, outer_box)

Compute the Voronoi cells of all center points `pts` with respect to all 9x replicated points in a 3x3 grid.
Uses spatial neighbor filtering and distance ordering for optimal performance.
"""
function compute_voronoi_cells_periodic(pts::Vector{Point2F}, outer_box::Vector{Point2F})
    N = length(pts)
    all_pts = create_3x3_replicated_points(pts)
    cells = Vector{Vector{Point2F}}(undef, N)
    search_radius = max(0.12, 4.5 / sqrt(N))

    for j in 1:N
        p_j = all_pts[j]
        poly = copy(outer_box)

        candidates = Int[]
        for k in 1:length(all_pts)
            k == j && continue
            p_k = all_pts[k]
            if abs(p_k[1] - p_j[1]) <= search_radius && abs(p_k[2] - p_j[2]) <= search_radius
                push!(candidates, k)
            end
        end
        # Sort closest neighbors first to rapidly shrink the polygon
        sort!(candidates, by = k -> (all_pts[k][1] - p_j[1])^2 + (all_pts[k][2] - p_j[2])^2)

        for k in candidates
            p_k = all_pts[k]
            M = (p_j + p_k) / 2.0
            n = p_k - p_j
            poly = clip_polygon_halfplane(poly, M, n)
            isempty(poly) && break
        end
        cells[j] = poly
    end
    return cells
end

"""
    polygon_area(poly)

Calculate the area of a 2D polygon using the Shoelace formula.
"""
function polygon_area(poly::Vector{Point2F})
    nv = length(poly)
    nv < 3 && return 0.0
    s = 0.0
    for i in 1:nv
        p1 = poly[i]
        p2 = poly[mod1(i + 1, nv)]
        s += p1[1] * p2[2] - p2[1] * p1[2]
    end
    return abs(s) / 2.0
end

"""
    compute_total_unique_edge_length(cells; eps=1e-8)

Calculate the total length of all unique Voronoi edges in the collection of cells.
"""
function compute_total_unique_edge_length(cells::Vector{Vector{Point2F}}; eps::Float64=1e-8)
    edges = Dict{Tuple{Tuple{Int,Int}, Tuple{Int,Int}}, Float64}()
    for cell in cells
        nv = length(cell)
        nv < 2 && continue
        for i in 1:nv
            p1 = cell[i]
            p2 = cell[mod1(i + 1, nv)]
            k1 = (round(Int, p1[1] / eps), round(Int, p1[2] / eps))
            k2 = (round(Int, p2[1] / eps), round(Int, p2[2] / eps))
            edge_key = k1 < k2 ? (k1, k2) : (k2, k1)
            edges[edge_key] = dist(p1, p2)
        end
    end
    return sum(values(edges))
end

"""
    render_voronoi_page!(cr, pts, cells, k_star, p_star, cw, ch, margin, scale_factor,
                         title_text, subtitle_text; highlight_k_star=true)

Render a single page showing the Voronoi diagram with pink frame, gray fills, yellow unit square, and blue prophet site.
"""
function render_voronoi_page!(
    cr::CairoContext,
    pts::Vector{Point2F},
    cells::Vector{Vector{Point2F}},
    k_star::Int,
    p_star::Point2F,
    cw::Int,
    ch::Int,
    margin::Float64,
    scale_factor::Float64,
    title_text::String,
    subtitle_text::String;
    highlight_k_star::Bool = true
)
    N = length(pts)
    L = compute_total_unique_edge_length(cells)

    # Calibrate point radius and boundary width to fill ≤ 25% of unit square (10% points + 15% boundaries)
    r_unit = sqrt(0.10 / (pi * N))
    w_unit = 0.15 / max(L, 1e-4)

    r_pt = r_unit * scale_factor
    w_pt = w_unit * scale_factor
    r_blue_pt = 1.35 * r_pt
    r_blue_halo_pt = 2.0 * r_pt

    Cairo.save(cr)

    # Page Background
    set_source_rgb(cr, 0.98, 0.98, 0.99)
    Cairo.rectangle(cr, 0, 0, cw, ch)
    Cairo.fill(cr)

    # Header Title
    set_source_rgb(cr, 0.12, 0.12, 0.2)
    Cairo.select_font_face(cr, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_BOLD)
    Cairo.set_font_size(cr, 16.5)
    Cairo.move_to(cr, margin, margin * 0.60)
    Cairo.show_text(cr, title_text)

    # Subtitle / Status Note
    Cairo.select_font_face(cr, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
    Cairo.set_font_size(cr, 11.5)
    Cairo.move_to(cr, margin, ch - margin * 0.45)
    set_source_rgb(cr, 0.3, 0.3, 0.35)
    Cairo.show_text(cr, subtitle_text)

    # Setup geometry transform: y-upwards, origin placed so [-0.1, 1.1]^2 is centered on canvas
    Cairo.translate(cr, margin + 0.1 * scale_factor, ch - margin - 0.1 * scale_factor)
    Cairo.scale(cr, scale_factor, -scale_factor)

    # 1. Draw outer square boundary around [-0.1, 1.1]^2 in light pink
    set_source_rgb(cr, 1.0, 0.72, 0.80)
    Cairo.set_line_width(cr, 3.0)
    Cairo.rectangle(cr, -0.1, -0.1, 1.2, 1.2)
    Cairo.stroke(cr)

    # 2. Fill all Voronoi cells of the middle square points in light gray
    set_source_rgba(cr, 0.88, 0.88, 0.91, 0.55)
    for cell in cells
        length(cell) < 3 && continue
        Cairo.new_path(cr)
        Cairo.move_to(cr, cell[1][1], cell[1][2])
        for v in cell[2:end]
            Cairo.line_to(cr, v[1], v[2])
        end
        Cairo.close_path(cr)
        Cairo.fill(cr)
    end

    # 3. Highlight prophet site cell in soft blue (if requested and present)
    if highlight_k_star && k_star <= length(cells)
        largest_cell = cells[k_star]
        if length(largest_cell) >= 3
            set_source_rgba(cr, 0.35, 0.65, 0.95, 0.50)
            Cairo.new_path(cr)
            Cairo.move_to(cr, largest_cell[1][1], largest_cell[1][2])
            for v in largest_cell[2:end]
                Cairo.line_to(cr, v[1], v[2])
            end
            Cairo.close_path(cr)
            Cairo.fill_preserve(cr)

            set_source_rgb(cr, 0.1, 0.4, 0.85)
            Cairo.set_line_width(cr, max(w_pt * 1.5, 1.5))
            Cairo.stroke(cr)
        end
    end

    # 4. Stroke all Voronoi cell boundaries of the middle square points in black
    set_source_rgb(cr, 0.0, 0.0, 0.0)
    Cairo.set_line_width(cr, w_pt)
    for cell in cells
        length(cell) < 2 && continue
        Cairo.new_path(cr)
        Cairo.move_to(cr, cell[1][1], cell[1][2])
        for v in cell[2:end]
            Cairo.line_to(cr, v[1], v[2])
        end
        Cairo.close_path(cr)
        Cairo.stroke(cr)
    end

    # 5. Draw the boundary of the unit square [0, 1]^2 in light yellow (after diagram is drawn)
    set_source_rgb(cr, 0.98, 0.86, 0.20)
    Cairo.set_line_width(cr, max(w_pt * 1.5, 2.0))
    Cairo.rectangle(cr, 0.0, 0.0, 1.0, 1.0)
    Cairo.stroke(cr)

    # 6. Draw only the points in the middle square with calibrated radius
    for j in 1:N
        (highlight_k_star && j == k_star) && continue
        pj = pts[j]
        set_source_rgb(cr, 0.85, 0.15, 0.15)
        Cairo.new_path(cr)
        Cairo.arc(cr, pj[1], pj[2], r_unit, 0.0, 2pi)
        Cairo.fill(cr)
    end

    # Draw prophet point p* in blue if present
    if highlight_k_star && k_star <= N
        if r_pt >= 1.5
            set_source_rgba(cr, 0.1, 0.4, 0.9, 0.30)
            Cairo.new_path(cr)
            Cairo.arc(cr, p_star[1], p_star[2], (r_blue_halo_pt / scale_factor), 0.0, 2pi)
            Cairo.fill(cr)
        end

        set_source_rgb(cr, 0.0, 0.3, 0.9)
        Cairo.new_path(cr)
        Cairo.arc(cr, p_star[1], p_star[2], (r_blue_pt / scale_factor), 0.0, 2pi)
        Cairo.fill_preserve(cr)

        set_source_rgb(cr, 0.0, 0.15, 0.6)
        Cairo.set_line_width(cr, max(w_pt, 1.0))
        Cairo.stroke(cr)
    end

    Cairo.restore(cr)
    Cairo.show_page(cr)
end

# ==============================================================================
# MODE 1: Prefix Voronoi Sequence
# ==============================================================================
function generate_mode1_pdf(
    filename::String;
    N::Int = 100,
    seed::Union{Int, Nothing} = nothing,
    cw::Int = 800,
    ch::Int = 800,
    margin::Float64 = 50.0
)
    if !isnothing(seed)
        Random.seed!(seed)
    end

    if N < 4
        throw(ArgumentError("N must be at least 4 for Mode 1 (requested N = $N)"))
    end

    pts = [point(rand(), rand()) for _ in 1:N]
    outer_box = [point(-0.1, -0.1), point(1.1, -0.1), point(1.1, 1.1), point(-0.1, 1.1)]

    println("Precomputing final Voronoi diagram for N = $N points...")
    final_cells = compute_voronoi_cells_periodic(pts, outer_box)
    final_areas = [polygon_area(c) for c in final_cells]
    k_star = argmax(final_areas)
    p_star = pts[k_star]
    max_area = final_areas[k_star]

    @printf("Prophet site identified: p_%d = (%.4f, %.4f) with final cell area = %.6f (%.2f%% of torus)\n",
            k_star, p_star[1], p_star[2], max_area, max_area * 100)

    mkpath(dirname(filename))
    surface = CairoPDFSurface(filename, cw, ch)
    cr = CairoContext(surface)

    usable_w = cw - 2 * margin
    usable_h = ch - 2 * margin
    scale_factor = min(usable_w, usable_h) / 1.2
    total_pages = N - 3

    println("Mode 1: Rendering $total_pages prefix pages (from 4 to $N points) into $filename...")
    for i in 4:N
        page_num = i - 3
        prefix_pts = pts[1:i]
        prefix_cells = compute_voronoi_cells_periodic(prefix_pts, outer_box)
        has_prophet = (i >= k_star)

        title = if i == N
            @sprintf("Page %d / %d: Final Voronoi Diagram (%d Points, Largest Cell: p_%d, Area: %.4f)",
                     page_num, total_pages, N, k_star, max_area)
        else
            @sprintf("Page %d / %d: Voronoi Diagram of Prefix p_1 ... p_%d (3x3 Grid)",
                     page_num, total_pages, i)
        end

        subtitle = if i < k_star
            @sprintf("Prophet target p_%d arrives at step %d", k_star, k_star)
        elseif i == k_star
            @sprintf("★ Prophet target p_%d arrived at this step (drawn in blue, current area: %.4f)",
                     k_star, polygon_area(prefix_cells[k_star]))
        else
            @sprintf("Prophet target p_%d is present (drawn in blue, current area: %.4f)",
                     k_star, polygon_area(prefix_cells[k_star]))
        end

        render_voronoi_page!(cr, prefix_pts, prefix_cells, k_star, p_star, cw, ch, margin, scale_factor,
                             title, subtitle; highlight_k_star = has_prophet)
    end

    Cairo.finish(surface)
    println("Successfully generated Mode 1 PDF ($total_pages pages): $filename")
end

# ==============================================================================
# MODE 2: Powers of 2 Calibrated Diagrams
# ==============================================================================
function generate_mode2_pdf(
    filename::String;
    powers::UnitRange{Int} = 3:12,
    seed::Union{Int, Nothing} = nothing,
    cw::Int = 800,
    ch::Int = 800,
    margin::Float64 = 50.0
)
    if !isnothing(seed)
        Random.seed!(seed)
    end

    mkpath(dirname(filename))
    surface = CairoPDFSurface(filename, cw, ch)
    cr = CairoContext(surface)

    usable_w = cw - 2 * margin
    usable_h = ch - 2 * margin
    scale_factor = min(usable_w, usable_h) / 1.2
    outer_box = [point(-0.1, -0.1), point(1.1, -0.1), point(1.1, 1.1), point(-0.1, 1.1)]

    total_pages = length(powers)
    println("Mode 2: Rendering $total_pages powers-of-2 pages (N = 2^$(first(powers))..2^$(last(powers))) into $filename...")

    for (page_idx, p) in enumerate(powers)
        N = 2^p
        pts = [point(rand(), rand()) for _ in 1:N]
        cells = compute_voronoi_cells_periodic(pts, outer_box)
        areas = [polygon_area(c) for c in cells]

        k_star = argmax(areas)
        p_star = pts[k_star]
        max_area = areas[k_star]
        L = compute_total_unique_edge_length(cells)

        r_unit = sqrt(0.10 / (pi * N))
        w_unit = 0.15 / L
        r_pt = r_unit * scale_factor
        w_pt = w_unit * scale_factor
        total_area_pct = (N * pi * r_unit^2 + L * w_unit) * 100.0

        @printf("  Page %2d / %2d: N = %4d (2^%2d) | L = %5.1f | r = %4.1f pt, w = %4.2f pt | Prophet: p_%d (Area: %.4f)\n",
                page_idx, total_pages, N, p, L, r_pt, w_pt, k_star, max_area)

        title = @sprintf("Page %d / %d: Final Voronoi Diagram (N = 2^%d = %d Points, 3x3 Grid)",
                         page_idx, total_pages, p, N)
        subtitle = @sprintf("Calibrated Styling (Ink Area = %.1f%% ≤ 25%%): Point Radius = %.1f pt, Boundary Width = %.2f pt | Largest Cell: p_%d (Area = %.4f)",
                            total_area_pct, r_pt, w_pt, k_star, max_area)

        render_voronoi_page!(cr, pts, cells, k_star, p_star, cw, ch, margin, scale_factor,
                             title, subtitle; highlight_k_star = true)
    end

    Cairo.finish(surface)
    println("Successfully generated Mode 2 PDF ($total_pages pages): $filename")
end

# ==============================================================================
# MODE 3: 2-Page Prophet Comparison (Default Mode)
# ==============================================================================
function generate_mode3_pdf(
    filename::String;
    N::Int = 1000,
    seed::Union{Int, Nothing} = nothing,
    cw::Int = 800,
    ch::Int = 800,
    margin::Float64 = 50.0
)
    if !isnothing(seed)
        Random.seed!(seed)
    end

    if N < 4
        throw(ArgumentError("N must be at least 4 for Mode 3 (requested N = $N)"))
    end

    pts = [point(rand(), rand()) for _ in 1:N]
    outer_box = [point(-0.1, -0.1), point(1.1, -0.1), point(1.1, 1.1), point(-0.1, 1.1)]

    # 1. Compute final Voronoi diagram for N points
    println("Mode 3: Precomputing final Voronoi diagram for N = $N points...")
    cells_N = compute_voronoi_cells_periodic(pts, outer_box)
    areas_N = [polygon_area(c) for c in cells_N]
    k_star = argmax(areas_N)
    p_star = pts[k_star]
    max_area_final = areas_N[k_star]

    @printf("Prophet site identified: p_%d = (%.4f, %.4f) with final cell area = %.6f (%.2f%% of torus)\n",
            k_star, p_star[1], p_star[2], max_area_final, max_area_final * 100)

    # 2. Compute prefix Voronoi diagram for the first k_star points
    println("Mode 3: Computing prefix Voronoi diagram for k = $k_star points (when prophet site p_$k_star arrived)...")
    prefix_pts = pts[1:k_star]
    cells_k = compute_voronoi_cells_periodic(prefix_pts, outer_box)
    area_at_arrival = polygon_area(cells_k[k_star])
    @printf("Cell area of prophet site p_%d at arrival (k = %d): %.6f (%.2f%% of torus)\n",
            k_star, k_star, area_at_arrival, area_at_arrival * 100)

    mkpath(dirname(filename))
    surface = CairoPDFSurface(filename, cw, ch)
    cr = CairoContext(surface)

    usable_w = cw - 2 * margin
    usable_h = ch - 2 * margin
    scale_factor = min(usable_w, usable_h) / 1.2

    # --- Page 1: Final Voronoi diagram of N points ---
    println("Rendering Page 1 / 2: Final Voronoi diagram of N = $N points...")
    title_p1 = @sprintf("Page 1 / 2: Final Voronoi Diagram of N = %d Points (3x3 Grid)", N)
    subtitle_p1 = @sprintf("Prophet Site p_%d (blue) has the Maximum Final Cell Area: %.6f (%.2f%% of torus)",
                           k_star, max_area_final, max_area_final * 100)
    render_voronoi_page!(cr, pts, cells_N, k_star, p_star, cw, ch, margin, scale_factor,
                         title_p1, subtitle_p1; highlight_k_star = true)

    # --- Page 2: Prefix Voronoi diagram of the first k points ---
    println("Rendering Page 2 / 2: Prefix Voronoi diagram of k = $k_star points...")
    title_p2 = @sprintf("Page 2 / 2: Voronoi Diagram of Prefix p_1 ... p_%d at Arrival Step", k_star)
    subtitle_p2 = @sprintf("Prophet Site p_%d Cell Area at Arrival: %.6f (%.2f%%) vs Final Area at Step %d: %.6f (%.2f%%)",
                           k_star, area_at_arrival, area_at_arrival * 100, N, max_area_final, max_area_final * 100)
    render_voronoi_page!(cr, prefix_pts, cells_k, k_star, p_star, cw, ch, margin, scale_factor,
                         title_p2, subtitle_p2; highlight_k_star = true)

    Cairo.finish(surface)
    println("Successfully generated Mode 3 PDF (2 pages): $filename")
end

# ==============================================================================
# CLI Entrypoint
# ==============================================================================
function main(args=ARGS)
    mode_str = length(args) >= 1 ? lowercase(args[1]) : "3"

    outdir = normpath(joinpath(@__DIR__, "..", "output"))
    mkpath(outdir)
    pdf_path = joinpath(outdir, "prophet.pdf")

    println("="^70)
    if mode_str in ["1", "mode1", "prefix"]
        N = length(args) >= 2 ? parse(Int, args[2]) : 100
        seed = length(args) >= 3 ? parse(Int, args[3]) : nothing
        println("Prophet Voronoi Diagram - Mode 1 (Prefix Sequence, N = $N)")
        println("="^70)
        generate_mode1_pdf(pdf_path; N=N, seed=seed)
    elseif mode_str in ["2", "mode2", "powers"]
        seed = length(args) >= 2 ? parse(Int, args[2]) : nothing
        println("Prophet Voronoi Diagram - Mode 2 (Powers of 2: N = 2^3..2^12)")
        println("="^70)
        generate_mode2_pdf(pdf_path; powers=3:12, seed=seed)
    elseif mode_str in ["3", "mode3", "prophet"] || isempty(args)
        N = length(args) >= 2 ? parse(Int, args[2]) : 1000
        seed = length(args) >= 3 ? parse(Int, args[3]) : nothing
        println("Prophet Voronoi Diagram - Mode 3 [Default] (2-Page Prophet Comparison, N = $N)")
        println("="^70)
        generate_mode3_pdf(pdf_path; N=N, seed=seed)
    else
        println("Unknown mode: '$mode_str'. Available modes: 1 (prefix), 2 (powers), 3 (default prophet comparison).")
        println("Falling back to default Mode 3 (N = 1000)...")
        println("="^70)
        generate_mode3_pdf(pdf_path; N=1000)
    end
    println("="^70)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
