#!/usr/bin/env julia

using Pkg
Pkg.activate(@__DIR__)

"""
    Example: Prophet Voronoi Diagram (Periodic 3x3 Grid)

Supported Modes:
  Mode 1: Prefix Voronoi sequence (starts from 4 points up to N points, default N = 100).
  Mode 2: Powers of 2 comparison (N = 2^3 ... 2^12, 10 pages) with ≤ 25% ink area calibration.
  Mode 3: (Default) 7-page Prophet comparison:
          - Page 1: Final Voronoi diagram of N = 1000 points (Cairo, no text).
          - Page 2: Text description of Page 1 (LuaLaTeX).
          - Page 3: Prefix Voronoi diagram of {p_1, ..., p_k} at prophet arrival (Cairo, no text).
          - Page 4: Text description of Page 3 (LuaLaTeX).
          - Page 5: Mathematical formulations and prophet inequality background (LuaLaTeX).
          - Page 6: Final Voronoi diagram without sites, highlighting in red all cells with area ≥ 1/2 max area (Cairo).
          - Page 7: Final Voronoi diagram without sites, highlighting in red all cells with area ≥ 1/4 max area (Cairo).

Usage:
  ./examples/prophet.jl             # Runs default Mode 3 (N = 1000, 7 pages)
  ./examples/prophet.jl 1 [N] [seed] # Mode 1: Prefix Voronoi sequence
  ./examples/prophet.jl 2 [seed]     # Mode 2: Powers of 2 (N = 2^3..2^12, 10 pages)
  ./examples/prophet.jl 3 [N] [seed] # Mode 3: 7-page Prophet comparison (default N = 1000)
"""

using BasicCompGeometry
using LinearAlgebra
using Cairo
using Printf
using Random
using LaTeXStrings

# Register session LaTeX macros for figures and documentation pages
add_latex_packages!("amsmath", "amssymb", "amsfonts", "microtype")
add_latex_macros!(raw"""
\newcommand{\R}{\mathbb{R}}
\newcommand{\E}{\mathbb{E}}
\newcommand{\eps}{\varepsilon}
\DeclareMathOperator{\Vor}{Vor}
\DeclareMathOperator{\Area}{Area}
""")

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
                         title_text="", subtitle_text=""; highlight_k_star=true, show_text=false)

Render a single page showing the Voronoi diagram with pink frame, gray fills, yellow unit square, and blue prophet site.
When `show_text` is `false`, the title header and subtitle footer are skipped.
"""
function render_voronoi_page!(
    cr,
    pts::Vector{Point2F},
    cells::Vector{Vector{Point2F}},
    k_star::Int,
    p_star::Point2F,
    cw::Int,
    ch::Int,
    margin::Float64,
    scale_factor::Float64,
    title_text::String = "",
    subtitle_text::String = "";
    highlight_k_star::Bool = true,
    show_text::Bool = false,
    show_points::Bool = true,
    highlight_cells::Vector{Int} = Int[],
    highlight_color::NTuple{4, Float64} = (0.92, 0.22, 0.22, 0.55),
    highlight_stroke_color::NTuple{3, Float64} = (0.80, 0.10, 0.10)
)
    N = length(pts)
    L = compute_total_unique_edge_length(cells)

    # Calibrate point radius and boundary width to fill ≤ 25% of unit square (10% points + 15% boundaries)
    r_unit = sqrt(0.10 / (pi * N))
    w_unit = 0.15 / max(L, 1e-4)

    r_pt = r_unit * scale_factor
    w_pt = w_unit * scale_factor

    # Point size (diameter) should never exceed 2% of the width of the image
    max_radius_pt = 0.01 * cw
    if r_pt > max_radius_pt
        r_pt = max_radius_pt
        r_unit = r_pt / scale_factor
    end

    r_blue_pt = min(1.35 * r_pt, max_radius_pt * 1.35)
    r_blue_halo_pt = min(2.0 * r_pt, max_radius_pt * 2.0)

    Cairo.save(cr)

    # Page Background
    set_source_rgb(cr, 0.98, 0.98, 0.99)
    Cairo.rectangle(cr, 0, 0, cw, ch)
    Cairo.fill(cr)

    # Header Title
    if show_text && !isempty(title_text)
        set_source_rgb(cr, 0.12, 0.12, 0.2)
        Cairo.select_font_face(cr, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_BOLD)
        Cairo.set_font_size(cr, 16.5)
        Cairo.move_to(cr, margin, margin * 0.60)
        Cairo.show_text(cr, title_text)
    end

    # Subtitle / Status Note
    if show_text && !isempty(subtitle_text)
        Cairo.select_font_face(cr, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
        Cairo.set_font_size(cr, 11.5)
        Cairo.move_to(cr, margin, ch - margin * 0.45)
        set_source_rgb(cr, 0.3, 0.3, 0.35)
        Cairo.show_text(cr, subtitle_text)
    end

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

    # 3. Highlight specific cells (if provided)
    if !isempty(highlight_cells)
        set_source_rgba(cr, highlight_color[1], highlight_color[2], highlight_color[3], highlight_color[4])
        for idx in highlight_cells
            (idx < 1 || idx > length(cells)) && continue
            cell = cells[idx]
            length(cell) < 3 && continue
            Cairo.new_path(cr)
            Cairo.move_to(cr, cell[1][1], cell[1][2])
            for v in cell[2:end]
                Cairo.line_to(cr, v[1], v[2])
            end
            Cairo.close_path(cr)
            Cairo.fill_preserve(cr)

            set_source_rgb(cr, highlight_stroke_color[1], highlight_stroke_color[2], highlight_stroke_color[3])
            Cairo.set_line_width(cr, max(w_pt * 1.2, 1.2))
            Cairo.stroke(cr)
        end
    end

    # 4. Highlight prophet cell in blue (if requested)
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

    # 5. Stroke all Voronoi cell boundaries of the middle square points in black
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

    # 6. Draw points in the middle square (if show_points is true)
    if show_points
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
    end

    Cairo.restore(cr)
end

# ==============================================================================
# MODE 1: Prefix Voronoi Sequence
# ==============================================================================
function _render_mode1_content!(canvas, pts, outer_box, k_star, p_star, max_area, cw, ch, margin, total_pages, N)
    usable_w = cw - 2 * margin
    usable_h = ch - 2 * margin
    scale_factor = min(usable_w, usable_h) / 1.2

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

        render_voronoi_page!(canvas, prefix_pts, prefix_cells, k_star, p_star, cw, ch, margin, scale_factor,
                             title, subtitle; highlight_k_star = has_prophet, show_text = true)
        description(canvas, subtitle)
        if i < N
            Cairo.show_page(canvas)
        end
    end
end

function generate_mode1(
    pdf_filename::String,
    svg_dir::String;
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

    mkpath(dirname(pdf_filename))
    total_pages = N - 3

    # 1. Multi-page PDF via Canvas
    println("Mode 1: Rendering $total_pages prefix pages into PDF via Canvas: $pdf_filename...")
    open_canvas(pdf_filename, cw, ch) do canvas
        _render_mode1_content!(canvas, pts, outer_box, k_star, p_star, max_area, cw, ch, margin, total_pages, N)
    end
    println("Successfully generated Mode 1 PDF ($total_pages pages): $pdf_filename")

    # 2. HTML Presentation via Canvas
    html_target = joinpath(svg_dir, "index.html")
    println("Mode 1: Rendering $total_pages HTML presentation slides via Canvas: $html_target...")
    c_html = open_canvas(html_target, cw, ch; title="Prophet Voronoi Sequence (Mode 1)") do canvas
        _render_mode1_content!(canvas, pts, outer_box, k_star, p_star, max_area, cw, ch, margin, total_pages, N)
    end
    println("Generated $total_pages SVG slides in: $svg_dir")
    println("HTML Presentation: ", get_file_path(c_html))
end

# ==============================================================================
# MODE 2: Powers of 2 Calibrated Diagrams
# ==============================================================================
function _render_mode2_content!(canvas, powers, seed, cw, ch, margin)
    usable_w = cw - 2 * margin
    usable_h = ch - 2 * margin
    scale_factor = min(usable_w, usable_h) / 1.2
    outer_box = [point(-0.1, -0.1), point(1.1, -0.1), point(1.1, 1.1), point(-0.1, 1.1)]
    total_pages = length(powers)

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

        max_radius_pt = 0.01 * cw
        if r_pt > max_radius_pt
            r_pt = max_radius_pt
            r_unit = r_pt / scale_factor
        end

        total_area_pct = (N * pi * r_unit^2 + L * w_unit) * 100.0

        @printf("  Page %2d / %2d: N = %4d (2^%2d) | L = %5.1f | r = %4.1f pt, w = %4.2f pt | Prophet: p_%d (Area: %.4f)\n",
                page_idx, total_pages, N, p, L, r_pt, w_pt, k_star, max_area)

        title = @sprintf("Page %d / %d: Final Voronoi Diagram (N = 2^%d = %d Points, 3x3 Grid)",
                         page_idx, total_pages, p, N)
        subtitle = @sprintf("Calibrated Styling (Ink Area = %.1f%% ≤ 25%%): Point Radius = %.1f pt, Boundary Width = %.2f pt | Largest Cell: p_%d (Area = %.4f)",
                            total_area_pct, r_pt, w_pt, k_star, max_area)

        render_voronoi_page!(canvas, pts, cells, k_star, p_star, cw, ch, margin, scale_factor,
                             title, subtitle; highlight_k_star = true, show_text = true)
        description(canvas, subtitle)
        if page_idx < total_pages
            Cairo.show_page(canvas)
        end
    end
end

function generate_mode2(
    pdf_filename::String,
    svg_dir::String;
    powers::UnitRange{Int} = 3:12,
    seed::Union{Int, Nothing} = nothing,
    cw::Int = 800,
    ch::Int = 800,
    margin::Float64 = 50.0
)
    if !isnothing(seed)
        Random.seed!(seed)
    end

    mkpath(dirname(pdf_filename))
    total_pages = length(powers)

    # 1. Multi-page PDF via Canvas
    println("Mode 2: Rendering $total_pages powers-of-2 pages into PDF via Canvas: $pdf_filename...")
    open_canvas(pdf_filename, cw, ch) do canvas
        _render_mode2_content!(canvas, powers, seed, cw, ch, margin)
    end
    println("Successfully generated Mode 2 PDF ($total_pages pages): $pdf_filename")

    # 2. HTML Presentation via Canvas
    html_target = joinpath(svg_dir, "index.html")
    println("Mode 2: Rendering $total_pages HTML presentation slides via Canvas: $html_target...")
    c_html = open_canvas(html_target, cw, ch; title="Prophet Voronoi Powers of 2 (Mode 2)") do canvas
        _render_mode2_content!(canvas, powers, seed, cw, ch, margin)
    end
    println("Generated $total_pages SVG slides in: $svg_dir")
    println("HTML Presentation: ", get_file_path(c_html))
end

# ==============================================================================
# LaTeX Text Generation Helpers
# ==============================================================================

"""
    build_text_page2(N, k_star, max_area)

Generate LaTeX string for Page 2 describing the Final Voronoi Diagram on Page 1.
"""
function build_text_page2(N::Int, k_star::Int, max_area::Float64)
    expected_area = 1.0 / N
    ratio = max_area / expected_area
    pct_max = max_area * 100.0
    pct_exp = expected_area * 100.0

    return """
\\section*{Page 1: Final Voronoi Diagram on the Flat Torus}

\\noindent
\\textbf{Configuration:} \$N = $N\$ points sampled independently and uniformly at random from the unit square \$[0, 1]^2\$:
\\[
p_1, p_2, \\dots, p_N \\stackrel{\\text{i.i.d.}}{\\sim} \\mathrm{Uniform}([0, 1]^2).
\\]
Periodic boundary conditions are enforced on the 2D flat torus \$\\mathbb{T}^2 \\cong [0, 1)^2\$ by replicating the unit square across a \$3 \\times 3\$ grid of adjacent domains.

\\subsection*{The Prophet Site \$p_{k^*}\$}

\\noindent
An offline \\emph{prophet} with complete hindsight knowledge of all \$N\$ site coordinates identifies the site index \$k^*\$ whose final Voronoi cell achieves maximum area:
\\[
k^* = \\arg\\max_{1 \\le i \\le N} \\operatorname{Area}\\bigl(V_N(p_i)\\bigr).
\\]

\\begin{itemize}
  \\item \\textbf{Prophet Site Index:} \$k^* = $k_star\$
  \\item \\textbf{Final Cell Area:} \$\\operatorname{Area}(V_N(p_{k^*})) = $(round(max_area, digits=6))\$ ($(round(pct_max, digits=4))\\% of total torus area)
  \\item \\textbf{Expected Cell Area:} \$\\mathbb{E}[\\operatorname{Area}(V_N)] = 1/N = $(round(expected_area, digits=6))\$ ($(round(pct_exp, digits=4))\\%)
  \\item \\textbf{Ratio over Mean:} \$\\operatorname{Area}(V_N(p_{k^*})) / (1/N) = $(round(ratio, digits=2))\\times\$
\\end{itemize}

\\subsection*{Visual Encoding}
\\begin{itemize}
  \\item \\textbf{Golden Boundary:} Outlines the fundamental unit domain \$[0, 1]^2\$.
  \\item \\textbf{Pink Margin:} Outlines the rendered \$[-0.1, 1.1]^2\$ bounding region for periodic boundary continuity.
  \\item \\textbf{Gray Polygons:} The \$N\$ periodic Voronoi cells \$V_N(p_1), \\dots, V_N(p_N)\$.
  \\item \\textbf{Blue Highlight:} The prophet cell \$V_N(p_{k^*})\$ with highlighted site marker.
\\end{itemize}
"""
end

"""
    build_text_page4(k_star, area_at_arrival, N, max_area_final)

Generate LaTeX string for Page 4 describing the Prefix Voronoi Diagram on Page 3.
"""
function build_text_page4(k_star::Int, area_at_arrival::Float64, N::Int, max_area_final::Float64, num_larger::Int=0)
    pct_arrival = area_at_arrival * 100.0
    pct_final = max_area_final * 100.0
    retention_pct = (max_area_final / area_at_arrival) * 100.0
    shrink_pct = 100.0 - retention_pct
    k_prev = k_star - 1

    return """
\\section*{Pages 3--4: Prefix Voronoi Diagram at Arrival Step \$k^*\$ }

\\noindent
\\textbf{Prefix Configuration:} At step \$k = $k_star\$, the prophet site \$p_{k^*}\$ arrives into the stream of sites \$\\{p_1, \\dots, p_{$k_prev}\\}\$. The diagrams on Pages 3 and 4 illustrate the state of the Voronoi partition on the torus \$\\mathbb{T}^2\$ using only the first \$k^*\$ points:
\\[
\\mathcal{P}_{k^*} = \\{p_1, p_2, \\dots, p_{k^*}\\}.
\\]

\\subsection*{Cell Evolution and Shrinkage}

\\noindent
When site \$p_{k^*}\$ is first inserted, it only competes with the preceding \$k^* - 1\$ points. As subsequent points \$p_{k^*+1}, \\dots, p_N\$ arrive, each new site clips existing cells via bisecting halfplanes:
\\[
V_m(p_{k^*}) = V_{m-1}(p_{k^*}) \\cap H(p_{k^*}, p_m), \\quad \\text{for } m = k^*+1, \\dots, N.
\\]

\\begin{itemize}
  \\item \\textbf{Arrival Step:} \$k^* = $k_star\$ of \$N = $N\$
  \\item \\textbf{Cell Area at Arrival:} \$\\operatorname{Area}(V_{k^*}(p_{k^*})) = $(round(area_at_arrival, digits=6))\$ ($(round(pct_arrival, digits=4))\\% of torus)
  \\item \\textbf{Final Cell Area:} \$\\operatorname{Area}(V_N(p_{k^*})) = $(round(max_area_final, digits=6))\$ ($(round(pct_final, digits=4))\\% of torus)
  \\item \\textbf{Area Retention Ratio:} \$\\rho = \\frac{\\operatorname{Area}(V_N(p_{k^*}))}{\\operatorname{Area}(V_{k^*}(p_{k^*}))} = $(round(retention_pct, digits=2))\\%\$
  \\item \\textbf{Page 4 (Green Highlight):} Displays the $(num_larger) cells in the prefix configuration that are strictly larger in area than the prophet cell \$p_{k^*}\$ at arrival.
\\end{itemize}

\\subsection*{Key Geometric Insights}
\\begin{itemize}
  \\item The cell area sequence \$m \\mapsto \\operatorname{Area}(V_m(p_{k^*}))\$ is monotonically non-increasing.
  \\item Despite losing $(round(shrink_pct, digits=2))\\% of its initial area to later arrivals, \$p_{k^*}\$ retains sufficient volume to emerge as the maximum-area cell across the entire final configuration.
\\end{itemize}
"""
end

"""
    build_text_page5()

Generate LaTeX string for Page 5 containing mathematical formulations and prophet inequality foundations.
"""
function build_text_page5()
    return """
\\section*{Theoretical Background \\& Mathematical Foundations}

\\subsection*{1. Classical Prophet Inequality}
\\noindent
Let \$X_1, X_2, \\dots, X_N\$ be independent non-negative random variables drawn from known distributions \$D_1, \\dots, D_N\$. 
An offline prophet with complete foresight achieves reward \$\\mathbb{E}[\\max_{1 \\le i \\le N} X_i]\$. 
A sequential gambler selecting an online stopping time \$\\tau\$ with threshold strategies satisfies the classic Krengel-Sucheston theorem (1977):
\\[
\\mathbb{E}\\left[ \\max_{1 \\le i \\le N} X_i \\right] \\le 2 \\cdot \\mathbb{E}[X_\\tau].
\\]

\\subsection*{2. Poisson-Voronoi Cell Area Distribution}
\\noindent
For a homogeneous Poisson point process in \$\\mathbb{R}^2\$ with intensity \$\\lambda\$, the mean cell area is \$\\mathbb{E}[A] = 1/\\lambda\$.
The normalized cell area \$S = A / \\mathbb{E}[A]\$ closely follows a generalized Gamma distribution:
\\[
f_S(s) = \\frac{c \\, b^a}{\\Gamma(a)} s^{a c - 1} \\exp\\left(-b s^c\\right), \\quad a \\approx 3.5, \\; b \\approx 3.5, \\; c \\approx 1.07.
\\]

\\subsection*{3. Extremal Statistics of Voronoi Cells}
\\noindent
For \$N\$ uniform random points in the 2D torus \$[0, 1)^2\$, the maximum cell area satisfies:
\\[
\\max_{1 \\le i \\le N} \\operatorname{Area}(V_N(p_i)) = \\Theta\\left( \\frac{\\ln N}{N} \\right) \\quad \\text{with high probability as } N \\to \\infty.
\\]

\\subsection*{4. Geometric Online Selection}
\\noindent
In geometric prophet settings, sites arrive sequentially, and the decision-maker must select or irrevocably evaluate geometric functionals (e.g., Voronoi territory, nearest neighbor distance, Delaunay edge length) in an online fashion under unknown future arrivals.
"""
end

# ==============================================================================
# MODE 3: 8-Page Prophet Comparison (Default Mode)
# ==============================================================================
function _render_mode3_content!(canvas, pts, cells_N, prefix_pts, cells_k, k_star, p_star,
                                max_area_final, area_at_arrival, large_cells_half, large_cells_quarter,
                                N, cw, ch, margin, scale_factor)
    areas_k = [polygon_area(c) for c in cells_k]
    larger_cells_k = findall(a -> a > area_at_arrival + 1e-9, areas_k)
    num_larger = length(larger_cells_k)

    # Page 1: Diagram 1
    println("  Rendering Page 1 / 8: Final Voronoi Diagram...")
    render_voronoi_page!(canvas, pts, cells_N, k_star, p_star, cw, ch, margin, scale_factor;
                         highlight_k_star = true, show_text = false)
    description(canvas, "Page 1: Final Voronoi diagram on the flat torus for N = $N points. Prophet site p_$k_star achieves maximum cell area $(round(max_area_final, digits=4)).")
    Cairo.show_page(canvas)

    # Page 2: Text 1
    println("  Rendering Page 2 / 8: Text documentation for Page 1...")
    cairo_draw_latex_page(canvas, build_text_page2(N, k_star, max_area_final))
    description(canvas, "Page 2: Mathematical formulation and statistics of the final Voronoi configuration.")
    Cairo.show_page(canvas)

    # Page 3: Diagram 2
    println("  Rendering Page 3 / 8: Prefix Voronoi Diagram at arrival step k = $k_star...")
    render_voronoi_page!(canvas, prefix_pts, cells_k, k_star, p_star, cw, ch, margin, scale_factor;
                         highlight_k_star = true, show_text = false)
    description(canvas, "Page 3: Prefix Voronoi diagram of the first k = $k_star points when prophet site p_$k_star arrived.")
    Cairo.show_page(canvas)

    # Page 4: Diagram 3 (NEW: Prefix cells larger in area than blue cell in green)
    println("  Rendering Page 4 / 8: Prefix Voronoi Diagram (cells larger than prophet cell in green)...")
    render_voronoi_page!(canvas, prefix_pts, cells_k, k_star, p_star, cw, ch, margin, scale_factor;
                         highlight_k_star = true, show_text = false, show_points = true,
                         highlight_cells = larger_cells_k,
                         highlight_color = (0.20, 0.75, 0.30, 0.55),
                         highlight_stroke_color = (0.10, 0.60, 0.20))
    description(canvas, "Page 4: Prefix Voronoi diagram at step k = $k_star highlighting $(num_larger) cells with area larger than the prophet cell (green).")
    Cairo.show_page(canvas)

    # Page 5: Text 2
    println("  Rendering Page 5 / 8: Text documentation for Pages 3 & 4...")
    cairo_draw_latex_page(canvas, build_text_page4(k_star, area_at_arrival, N, max_area_final, num_larger))
    description(canvas, "Page 5: Analysis of cell area evolution, shrinkage, and prefix area rankings.")
    Cairo.show_page(canvas)

    # Page 6: Text 3
    println("  Rendering Page 6 / 8: Theoretical Background & Mathematical Foundations...")
    cairo_draw_latex_page(canvas, build_text_page5())
    description(canvas, "Page 6: Classical prophet inequality theorem and Poisson-Voronoi extremal bounds.")
    Cairo.show_page(canvas)

    # Page 7: Diagram 4
    println("  Rendering Page 7 / 8: Final Voronoi diagram (cells >= 1/2 max area in red)...")
    render_voronoi_page!(canvas, pts, cells_N, k_star, p_star, cw, ch, margin, scale_factor;
                         highlight_k_star = false, show_text = false, show_points = false,
                         highlight_cells = large_cells_half,
                         highlight_color = (0.92, 0.22, 0.22, 0.55))
    description(canvas, "Page 7: Final Voronoi diagram highlighting $(length(large_cells_half)) cells with area ≥ 1/2 of maximum area (red).")
    Cairo.show_page(canvas)

    # Page 8: Diagram 5
    println("  Rendering Page 8 / 8: Final Voronoi diagram (cells >= 1/4 max area in red)...")
    render_voronoi_page!(canvas, pts, cells_N, k_star, p_star, cw, ch, margin, scale_factor;
                         highlight_k_star = false, show_text = false, show_points = false,
                         highlight_cells = large_cells_quarter,
                         highlight_color = (0.92, 0.22, 0.22, 0.55))
    description(canvas, "Page 8: Final Voronoi diagram highlighting $(length(large_cells_quarter)) cells with area ≥ 1/4 of maximum area (red).")
end

function generate_mode3(
    pdf_filename::String,
    svg_dir::String;
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

    # 3. Identify cells by area threshold
    threshold_half = 0.5 * max_area_final
    large_cells_half = findall(a -> a >= threshold_half, areas_N)
    @printf("Mode 3: Found %d cells with area >= 1/2 max area (>= %.6f)\n", length(large_cells_half), threshold_half)

    threshold_quarter = 0.25 * max_area_final
    large_cells_quarter = findall(a -> a >= threshold_quarter, areas_N)
    @printf("Mode 3: Found %d cells with area >= 1/4 max area (>= %.6f)\n", length(large_cells_quarter), threshold_quarter)

    mkpath(dirname(pdf_filename))

    usable_w = cw - 2 * margin
    usable_h = ch - 2 * margin
    scale_factor = min(usable_w, usable_h) / 1.2

    # 1. Multi-page PDF via Canvas
    println("Mode 3: Rendering 8 pages into PDF via Canvas: $pdf_filename...")
    open_canvas(pdf_filename, cw, ch) do canvas
        _render_mode3_content!(canvas, pts, cells_N, prefix_pts, cells_k, k_star, p_star,
                               max_area_final, area_at_arrival, large_cells_half, large_cells_quarter,
                               N, cw, ch, margin, scale_factor)
    end
    println("Successfully generated Mode 3 PDF (8 pages): $pdf_filename")

    # 2. HTML Presentation via Canvas
    html_target = joinpath(svg_dir, "index.html")
    println("Mode 3: Rendering 8 HTML presentation slides via Canvas: $html_target...")
    c_html = open_canvas(html_target, cw, ch; title="Prophet Voronoi Presentation (Mode 3)") do canvas
        _render_mode3_content!(canvas, pts, cells_N, prefix_pts, cells_k, k_star, p_star,
                               max_area_final, area_at_arrival, large_cells_half, large_cells_quarter,
                               N, cw, ch, margin, scale_factor)
    end
    println("Generated 8 SVG slides in: $svg_dir")
    println("HTML Presentation: ", get_file_path(c_html))
end

# Backward compatibility aliases
const generate_mode1_pdf = (f; kwargs...) -> generate_mode1(f, normpath(joinpath(dirname(f), "svg")); kwargs...)
const generate_mode2_pdf = (f; kwargs...) -> generate_mode2(f, normpath(joinpath(dirname(f), "svg")); kwargs...)
const generate_mode3_pdf = (f; kwargs...) -> generate_mode3(f, normpath(joinpath(dirname(f), "svg")); kwargs...)

# ==============================================================================
# CLI Entrypoint
# ==============================================================================
function main(args=ARGS)
    # Check for help
    if any(a -> a in ("-h", "--help", "help"), args)
        println("""
Usage: prophet.jl [mode] [N]
       prophet.jl [N]

Arguments:
  mode   Mode number or name:
           3 or 'prophet' : 8-page prophet comparison (default)
           1 or 'prefix'  : Prefix Voronoi sequence
           2 or 'powers'  : Powers of 2 sequence (N = 2^3..2^12)
  N      Number of points (default: 100)

Examples:
  julia --project=examples examples/prophet.jl
  julia --project=examples examples/prophet.jl 50
  julia --project=examples examples/prophet.jl 1 50
  julia --project=examples examples/prophet.jl 2
""")
        return
    end

    # Extract optional --seed / -s flag if specified
    seed_arg = nothing
    filtered_args = String[]
    i = 1
    while i <= length(args)
        arg = args[i]
        if (arg == "--seed" || arg == "-s") && i < length(args)
            seed_arg = tryparse(Int, args[i+1])
            i += 2
        elseif startswith(arg, "--seed=")
            seed_arg = tryparse(Int, split(arg, '=', limit=2)[2])
            i += 1
        else
            push!(filtered_args, arg)
            i += 1
        end
    end

    outdir = normpath(joinpath(@__DIR__, "..", "output"))
    mkpath(outdir)
    pdf_path = joinpath(outdir, "prophet.pdf")
    svg_dir = joinpath(outdir, "svg")

    first_arg = length(filtered_args) >= 1 ? filtered_args[1] : ""
    mode_str = "3"
    N_arg = nothing

    if isempty(filtered_args)
        mode_str = "3"
    elseif tryparse(Int, first_arg) !== nothing && parse(Int, first_arg) > 3
        # e.g. prophet.jl 200
        mode_str = "3"
        N_arg = parse(Int, first_arg)
        if length(filtered_args) >= 2 && isnothing(seed_arg)
            seed_arg = tryparse(Int, filtered_args[2])
        end
    else
        mode_str = lowercase(first_arg)
        if length(filtered_args) >= 2
            if mode_str in ("2", "mode2", "powers") && isnothing(seed_arg)
                seed_arg = tryparse(Int, filtered_args[2])
            else
                N_arg = tryparse(Int, filtered_args[2])
                if length(filtered_args) >= 3 && isnothing(seed_arg)
                    seed_arg = tryparse(Int, filtered_args[3])
                end
            end
        end
    end

    println("="^70)
    if mode_str in ["1", "mode1", "prefix"]
        N = !isnothing(N_arg) ? N_arg : 100
        println("Prophet Voronoi Diagram - Mode 1 (Prefix Sequence, N = $N)")
        println("="^70)
        generate_mode1(pdf_path, svg_dir; N=N, seed=seed_arg)
    elseif mode_str in ["2", "mode2", "powers"]
        println("Prophet Voronoi Diagram - Mode 2 (Powers of 2: N = 2^3..2^12)")
        println("="^70)
        generate_mode2(pdf_path, svg_dir; powers=3:12, seed=seed_arg)
    elseif mode_str in ["3", "mode3", "prophet"]
        N = !isnothing(N_arg) ? N_arg : 100
        println("Prophet Voronoi Diagram - Mode 3 [Default] (8-Page Prophet Comparison, N = $N)")
        println("="^70)
        generate_mode3(pdf_path, svg_dir; N=N, seed=seed_arg)
    else
        N = !isnothing(N_arg) ? N_arg : 100
        println("Unknown mode: '$mode_str'. Available modes: 1 (prefix), 2 (powers), 3 (default prophet comparison).")
        println("Falling back to default Mode 3 (N = $N)...")
        println("="^70)
        generate_mode3(pdf_path, svg_dir; N=N, seed=seed_arg)
    end
    println("="^70)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
