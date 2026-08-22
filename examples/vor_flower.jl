#!/usr/bin/env julia

using Pkg
Pkg.activate(@__DIR__)

"""
    Example: Voronoi Flower (vor_flower.jl)

Generates a 10-page interleaved PDF visualization (`output/flower.pdf`) exploring concentric scaled
regular polygon sites and their Voronoi diagrams. Each figure page is immediately followed by a concise summary page:
- **Page 1**: Figure 1 — Voronoi diagram of site set \$T = \\bigcup_{i=0}^{10} S_i \\cup \\{\\text{origin}\\}\$.
- **Page 2**: Summary of Figure 1 (point sets, colors, roles).
- **Page 3**: Figure 2 — Voronoi diagram of \$T\$ with origin cell (amber) and its vertices \$U\$ in red.
- **Page 4**: Summary of Figure 2.
- **Page 5**: Figure 3 — Voronoi diagram of \$T \\cup U\$, showing the closer new vertices in blue.
- **Page 6**: Summary of Figure 3.
- **Page 7**: Figure 4 — Voronoi diagram of \$T \\cup U\$, with new Voronoi cells of \$U\$ filled in light green.
- **Page 8**: Summary of Figure 4.
- **Page 9**: Figure 5 — Original Voronoi of \$T\$ with empty disks at \$V\$ colored by distance to nearest site.
- **Page 10**: Summary of Figure 5.

Usage:
    ./examples/vor_flower.jl
"""

using BasicCompGeometry
using LinearAlgebra
using Cairo
using Printf

const PROJROOT = normpath(joinpath(@__DIR__, ".."))
using LaTeXStrings

# Register LaTeX macros & packages
add_latex_packages!("amsmath", "amssymb", "amsfonts", "microtype")
add_latex_macros!(raw"""
\newcommand{\R}{\mathbb{R}}
\newcommand{\eps}{\varepsilon}
\DeclareMathOperator{\Vor}{Vor}
\DeclareMathOperator{\dist}{dist}
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

"""
    compute_voronoi_cells(pts, bbox_poly)

Compute the Voronoi cell for each point in `pts`, restricted to `bbox_poly`.
"""
function compute_voronoi_cells(pts::Vector{Point2F}, bbox_poly::Vector{Point2F})
    N = length(pts)
    cells = Vector{Vector{Point2F}}(undef, N)
    for i in 1:N
        p_i = pts[i]
        poly = copy(bbox_poly)
        other_indices = [j for j in 1:N if j != i]
        sort!(other_indices, by = j -> dist_sq(p_i, pts[j]))
        for j in other_indices
            p_j = pts[j]
            M = (p_i + p_j) / 2.0
            n = p_j - p_i
            poly = clip_polygon_halfplane(poly, M, n)
            isempty(poly) && break
        end
        cells[i] = poly
    end
    return cells
end

"""
    collect_unique_vertices(cells; eps=1e-6)

Collect all unique vertices from a list of Voronoi polygon cells.
"""
function collect_unique_vertices(cells::Vector{Vector{Point2F}}; eps::Float64=1e-6)
    vertices = Point2F[]
    for cell in cells
        for v in cell
            if !any(u -> dist(u, v) < eps, vertices)
                push!(vertices, v)
            end
        end
    end
    return vertices
end

"""
    colormap_turbo(t)

Turbo / viridis inspired smooth colormap for a normalized parameter `t` in [0, 1].
Returns (r, g, b) with values in [0, 1].
"""
function colormap_turbo(t::Float64)
    t_clamped = clamp(t, 0.0, 1.0)
    # Smooth RGB multi-segment colormap from vibrant blue to teal, yellow, orange, and red
    r = clamp(1.5 * (t_clamped - 0.2) * (t_clamped > 0.2 ? 1.0 : 0.0) + 0.1 * (1.0 - t_clamped), 0.0, 1.0)
    g = clamp(sin(pi * t_clamped), 0.0, 1.0)
    b = clamp(1.2 * (0.8 - t_clamped) * (t_clamped < 0.8 ? 1.0 : 0.0) + 0.15, 0.0, 1.0)
    return (r, g, b)
end

"""
    draw_voronoi_diagram(cr, cells, sites, cw, ch, r_view;
                         fill_color=(0.96, 0.96, 0.98),
                         edge_color=(0.35, 0.35, 0.40),
                         edge_width=1.0,
                         site_color=(0.15, 0.15, 0.20),
                         site_radius=2.5)

Helper to draw the basic Voronoi cells and sites in the standard coordinate system.
"""
function draw_voronoi_diagram(cr, cells, sites, cw, ch, r_view;
                              fill_color=(0.96, 0.96, 0.98),
                              edge_color=(0.35, 0.35, 0.40),
                              edge_width=1.0,
                              site_color=(0.15, 0.15, 0.20),
                              site_radius=2.5)
    # Fill cells
    set_source_rgb(cr, fill_color[1], fill_color[2], fill_color[3])
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

    # Stroke edges
    set_source_rgb(cr, edge_color[1], edge_color[2], edge_color[3])
    Cairo.set_line_width(cr, edge_width)
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

    # Draw sites
    if !isempty(sites) && site_radius > 0
        set_source_rgb(cr, site_color[1], site_color[2], site_color[3])
        for p in sites
            Cairo.new_path(cr)
            Cairo.arc(cr, p[1], p[2], site_radius, 0.0, 2pi)
            Cairo.fill(cr)
        end
    end
end

"""
    generate_flower_pdf(output_pdf)

Generate the 6-page PDF document for the Voronoi Flower visualization.
"""
function generate_flower_pdf(output_pdf::String="output/flower.pdf")
    mkpath(dirname(abspath(output_pdf)))

    # 1. Construct Site Set T
    k = 10
    angles = [2 * pi * j / k for j in 0:k-1]
    S0 = [point(cos(a), sin(a)) for a in angles]

    T = Point2F[point(0.0, 0.0)]
    for i in 0:10
        scale = (1.1)^i
        for p in S0
            push!(T, point(scale * p[1], scale * p[2]))
        end
    end

    # Bounding box for Voronoi computation
    R_bound = 3.5
    bbox = Point2F[
        point(-R_bound, -R_bound),
        point(R_bound, -R_bound),
        point(R_bound, R_bound),
        point(-R_bound, R_bound)
    ]

    # Compute Voronoi cells of T
    cells_T = compute_voronoi_cells(T, bbox)
    verts_T = collect_unique_vertices(cells_T)
    origin_cell = cells_T[1]
    U = collect_unique_vertices([origin_cell])

    # Compute Voronoi diagram of T union U
    T_prime = vcat(T, U)
    cells_T_prime = compute_voronoi_cells(T_prime, bbox)
    all_v_prime = collect_unique_vertices(cells_T_prime)
    origin_cell_prime = cells_T_prime[1]

    # Find truly NEW vertices of Vor(T U U) that did not exist in Vor(T)
    new_vertices = Point2F[]
    for v in all_v_prime
        (abs(v[1]) >= R_bound - 1e-3 || abs(v[2]) >= R_bound - 1e-3) && continue
        if !any(u -> dist(u, v) < 1e-5, verts_T)
            push!(new_vertices, v)
        end
    end

    # Identify set V: new Voronoi vertices of Vor(T U U) not adjacent to the origin cell,
    # keeping the closer distance group to the origin (norm ≈ 0.724)
    V_all_non_origin = Point2F[]
    for v in new_vertices
        if !any(u -> dist(u, v) < 1e-5, origin_cell_prime)
            push!(V_all_non_origin, v)
        end
    end
    V = filter(v -> norm(v) < 0.80, V_all_non_origin)

    # For Page 4: new vertices, keeping only the closer distance group
    # Note: the new vertices not on origin cell have norm ≈ 0.724 (closer) vs 0.851 (farther)
    new_vertices_closer = filter(v -> norm(v) < 0.80, new_vertices)

    # Find nearest site distances for points in V (with respect to T_prime)
    dists_V = Float64[]
    for v in V
        min_d = minimum(dist(v, s) for s in T_prime)
        push!(dists_V, min_d)
    end
    min_dist_V = minimum(dists_V)
    max_dist_V = maximum(dists_V)

    # Canvas setup
    cw = 800.0
    ch = 800.0
    margin = 50.0
    view_radius = 2.85  # visible radius in world coordinates
    scale_factor = (cw - 2 * margin) / (2 * view_radius)

    # Temporary files for diagram pages and text page
    tmp_diagrams_pdf = joinpath(mktempdir(), "diagrams.pdf")
    tmp_text_pdf = joinpath(mktempdir(), "text_page.pdf")

    # --- Draw 5 Diagram Pages ---
    surf = Cairo.CairoPDFSurface(tmp_diagrams_pdf, cw, ch)
    cr = Cairo.CairoContext(surf)

    # Transform function for each diagram page
    function setup_page_transform(page_title::String="")
        # Background
        set_source_rgb(cr, 0.99, 0.99, 1.0)
        Cairo.rectangle(cr, 0, 0, cw, ch)
        Cairo.fill(cr)

        # Title
        if !isempty(page_title)
            set_source_rgb(cr, 0.12, 0.15, 0.25)
            Cairo.select_font_face(cr, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_BOLD)
            Cairo.set_font_size(cr, 18.0)
            Cairo.move_to(cr, margin, margin * 0.65)
            Cairo.show_text(cr, page_title)
        end

        # World Coordinate Setup: Center at (cw/2, ch/2), Y-upwards
        Cairo.translate(cr, cw / 2.0, ch / 2.0)
        Cairo.scale(cr, scale_factor, -scale_factor)
    end

    site_rad_world = 3.2 / scale_factor
    vtx_rad_world = 4.2 / scale_factor

    # Consistent Point Set Colors across all figures:
    # 1. Site set T: Dark Navy Blue (0.15, 0.20, 0.35)
    # 2. Vertex set U (sites of U when augmented): Bright Crimson Red (0.88, 0.10, 0.10)
    # 3. Vertex set V / new Voronoi vertices: Vibrant Royal Blue (0.10, 0.40, 0.92)
    COLOR_T = (0.15, 0.20, 0.35)
    COLOR_U = (0.88, 0.10, 0.10)
    COLOR_V = (0.10, 0.40, 0.92)

    # ================= PAGE 1 =================
    # Voronoi diagram of T
    Cairo.save(cr)
    setup_page_transform("Page 1: Voronoi Diagram of Set T")
    draw_voronoi_diagram(cr, cells_T, T, cw, ch, view_radius;
                         fill_color=(0.94, 0.95, 0.98),
                         edge_color=(0.40, 0.45, 0.55),
                         edge_width=1.0 / scale_factor,
                         site_color=COLOR_T,
                         site_radius=site_rad_world)
    Cairo.restore(cr)
    Cairo.show_page(cr)

    # ================= PAGE 2 =================
    # Voronoi diagram of T, highlighting origin cell and red vertices U
    Cairo.save(cr)
    setup_page_transform("Page 2: Origin Cell & Voronoi Vertices U (Red)")
    draw_voronoi_diagram(cr, cells_T, T, cw, ch, view_radius;
                         fill_color=(0.95, 0.95, 0.97),
                         edge_color=(0.50, 0.52, 0.60),
                         edge_width=0.9 / scale_factor,
                         site_color=COLOR_T,
                         site_radius=site_rad_world)

    # Highlight origin cell in soft gold/amber
    set_source_rgba(cr, 1.0, 0.85, 0.30, 0.65)
    Cairo.new_path(cr)
    Cairo.move_to(cr, origin_cell[1][1], origin_cell[1][2])
    for v in origin_cell[2:end]
        Cairo.line_to(cr, v[1], v[2])
    end
    Cairo.close_path(cr)
    Cairo.fill(cr)

    # Stroke origin cell boundary in dark amber
    set_source_rgb(cr, 0.85, 0.50, 0.05)
    Cairo.set_line_width(cr, 2.0 / scale_factor)
    Cairo.new_path(cr)
    Cairo.move_to(cr, origin_cell[1][1], origin_cell[1][2])
    for v in origin_cell[2:end]
        Cairo.line_to(cr, v[1], v[2])
    end
    Cairo.close_path(cr)
    Cairo.stroke(cr)

    # Draw vertices U in bright red with halo
    for u in U
        # Red halo
        set_source_rgba(cr, COLOR_U[1], COLOR_U[2], COLOR_U[3], 0.35)
        Cairo.new_path(cr)
        Cairo.arc(cr, u[1], u[2], vtx_rad_world * 1.8, 0.0, 2pi)
        Cairo.fill(cr)

        # Solid red point
        set_source_rgb(cr, COLOR_U[1], COLOR_U[2], COLOR_U[3])
        Cairo.new_path(cr)
        Cairo.arc(cr, u[1], u[2], vtx_rad_world, 0.0, 2pi)
        Cairo.fill(cr)

        # White border
        set_source_rgb(cr, 1.0, 1.0, 1.0)
        Cairo.set_line_width(cr, 1.0 / scale_factor)
        Cairo.stroke(cr)
    end
    Cairo.restore(cr)
    Cairo.show_page(cr)

    # ================= PAGE 4 (Diagram 3) =================
    # Voronoi diagram of T U U, new vertices (closer ones, V) in royal blue
    Cairo.save(cr)
    setup_page_transform("Page 4: Voronoi Diagram of T ∪ U (New Vertices V in Blue)")
    # Draw Voronoi cells and sites of T (dark navy)
    draw_voronoi_diagram(cr, cells_T_prime, T, cw, ch, view_radius;
                         fill_color=(0.95, 0.96, 0.98),
                         edge_color=(0.42, 0.46, 0.56),
                         edge_width=0.9 / scale_factor,
                         site_color=COLOR_T,
                         site_radius=site_rad_world)

    # Draw sites of U in consistent red
    for u in U
        set_source_rgb(cr, COLOR_U[1], COLOR_U[2], COLOR_U[3])
        Cairo.new_path(cr)
        Cairo.arc(cr, u[1], u[2], site_rad_world * 1.1, 0.0, 2pi)
        Cairo.fill(cr)
        set_source_rgb(cr, 1.0, 1.0, 1.0)
        Cairo.set_line_width(cr, 0.7 / scale_factor)
        Cairo.stroke(cr)
    end

    # Highlight origin cell in Vor(T U U)
    set_source_rgba(cr, 0.40, 0.80, 0.90, 0.45)
    Cairo.new_path(cr)
    Cairo.move_to(cr, origin_cell_prime[1][1], origin_cell_prime[1][2])
    for v in origin_cell_prime[2:end]
        Cairo.line_to(cr, v[1], v[2])
    end
    Cairo.close_path(cr)
    Cairo.fill(cr)

    # Draw new Voronoi vertices (closer ones, V) in royal blue
    for v in new_vertices_closer
        # Blue halo
        set_source_rgba(cr, COLOR_V[1], COLOR_V[2], COLOR_V[3], 0.30)
        Cairo.new_path(cr)
        Cairo.arc(cr, v[1], v[2], vtx_rad_world * 1.6, 0.0, 2pi)
        Cairo.fill(cr)

        # Bright blue vertex
        set_source_rgb(cr, COLOR_V[1], COLOR_V[2], COLOR_V[3])
        Cairo.new_path(cr)
        Cairo.arc(cr, v[1], v[2], vtx_rad_world, 0.0, 2pi)
        Cairo.fill(cr)

        # White outline
        set_source_rgb(cr, 1.0, 1.0, 1.0)
        Cairo.set_line_width(cr, 0.8 / scale_factor)
        Cairo.stroke(cr)
    end
    Cairo.restore(cr)
    Cairo.show_page(cr)

    # ================= PAGE 5 (Diagram 4) =================
    # Voronoi diagram of T U U with new cells (cells of U) filled in light green
    Cairo.save(cr)
    setup_page_transform("Page 5: Voronoi of T ∪ U (New Cells of U in Light Green)")
    draw_voronoi_diagram(cr, cells_T_prime, T, cw, ch, view_radius;
                         fill_color=(0.95, 0.96, 0.98),
                         edge_color=(0.42, 0.46, 0.56),
                         edge_width=0.9 / scale_factor,
                         site_color=COLOR_T,
                         site_radius=site_rad_world)

    # Fill new cells of U (indices length(T)+1 to length(T_prime)) in light green
    new_cell_indices = (length(T) + 1):length(T_prime)
    for idx in new_cell_indices
        cell = cells_T_prime[idx]
        length(cell) < 3 && continue
        set_source_rgba(cr, 0.55, 0.88, 0.55, 0.70)
        Cairo.new_path(cr)
        Cairo.move_to(cr, cell[1][1], cell[1][2])
        for v in cell[2:end]
            Cairo.line_to(cr, v[1], v[2])
        end
        Cairo.close_path(cr)
        Cairo.fill(cr)

        # Stroke green boundary
        set_source_rgb(cr, 0.20, 0.65, 0.25)
        Cairo.set_line_width(cr, 1.6 / scale_factor)
        Cairo.new_path(cr)
        Cairo.move_to(cr, cell[1][1], cell[1][2])
        for v in cell[2:end]
            Cairo.line_to(cr, v[1], v[2])
        end
        Cairo.close_path(cr)
        Cairo.stroke(cr)
    end

    # Draw sites of U in consistent red
    for u in U
        set_source_rgb(cr, COLOR_U[1], COLOR_U[2], COLOR_U[3])
        Cairo.new_path(cr)
        Cairo.arc(cr, u[1], u[2], site_rad_world * 1.2, 0.0, 2pi)
        Cairo.fill(cr)

        set_source_rgb(cr, 1.0, 1.0, 1.0)
        Cairo.set_line_width(cr, 0.8 / scale_factor)
        Cairo.stroke(cr)
    end
    Cairo.restore(cr)
    Cairo.show_page(cr)

    # ================= PAGE 6 (Diagram 5) =================
    # Original Voronoi of T + points V (blue) + semi-transparent disks colored by nearest site distance
    Cairo.save(cr)
    setup_page_transform("Page 6: Original Voronoi & Disks at V Colored by Nearest Distance")
    draw_voronoi_diagram(cr, cells_T, T, cw, ch, view_radius;
                         fill_color=(0.97, 0.97, 0.98),
                         edge_color=(0.55, 0.58, 0.65),
                         edge_width=0.8 / scale_factor,
                         site_color=COLOR_T,
                         site_radius=site_rad_world * 0.8)

    # Draw semi-transparent disks for each vertex in V
    dist_span = max(max_dist_V - min_dist_V, 1e-6)
    for (i, v) in enumerate(V)
        d = dists_V[i]
        t_norm = (d - min_dist_V) / dist_span
        r_col, g_col, b_col = colormap_turbo(t_norm)

        # Semi-transparent disk
        set_source_rgba(cr, r_col, g_col, b_col, 0.32)
        Cairo.new_path(cr)
        Cairo.arc(cr, v[1], v[2], d, 0.0, 2pi)
        Cairo.fill(cr)

        # Disk boundary stroke
        set_source_rgba(cr, r_col * 0.85, g_col * 0.85, b_col * 0.85, 0.60)
        Cairo.set_line_width(cr, 1.0 / scale_factor)
        Cairo.stroke(cr)
    end

    # Overlay points in V with consistent royal blue center and white halo
    for (i, v) in enumerate(V)
        set_source_rgb(cr, COLOR_V[1], COLOR_V[2], COLOR_V[3])
        Cairo.new_path(cr)
        Cairo.arc(cr, v[1], v[2], vtx_rad_world * 0.9, 0.0, 2pi)
        Cairo.fill(cr)

        set_source_rgb(cr, 1.0, 1.0, 1.0)
        Cairo.set_line_width(cr, 0.7 / scale_factor)
        Cairo.stroke(cr)
    end

    # Colorbar legend in the bottom right corner (in world coordinates)
    cb_x = view_radius * 0.25
    cb_y = -view_radius * 0.88
    cb_w = view_radius * 0.65
    cb_h = view_radius * 0.06
    n_bars = 50
    for b in 0:n_bars-1
        bx = cb_x + (b / n_bars) * cb_w
        bw = cb_w / n_bars + 1e-4
        bt = b / (n_bars - 1)
        r_c, g_c, b_c = colormap_turbo(bt)
        set_source_rgb(cr, r_c, g_c, b_c)
        Cairo.rectangle(cr, bx, cb_y, bw, cb_h)
        Cairo.fill(cr)
    end
    set_source_rgb(cr, 0.2, 0.2, 0.3)
    Cairo.set_line_width(cr, 0.8 / scale_factor)
    Cairo.rectangle(cr, cb_x, cb_y, cb_w, cb_h)
    Cairo.stroke(cr)

    # Colorbar text labels
    Cairo.select_font_face(cr, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_BOLD)
    Cairo.set_font_size(cr, 9.5 / scale_factor)
    Cairo.move_to(cr, cb_x, cb_y - 0.07)
    Cairo.restore(cr)

    # Draw colorbar labels in screen coordinates
    cb_screen_x = cw / 2.0 + cb_x * scale_factor
    cb_screen_y = ch / 2.0 - cb_y * scale_factor
    cb_screen_w = cb_w * scale_factor

    # Small background badge for colorbar label
    set_source_rgba(cr, 1.0, 1.0, 1.0, 0.85)
    Cairo.rectangle(cr, cb_screen_x - 6.0, cb_screen_y - 20.0, cb_screen_w + 12.0, 36.0)
    Cairo.fill(cr)
    set_source_rgba(cr, 0.7, 0.7, 0.8, 0.6)
    Cairo.set_line_width(cr, 1.0)
    Cairo.rectangle(cr, cb_screen_x - 6.0, cb_screen_y - 20.0, cb_screen_w + 12.0, 36.0)
    Cairo.stroke(cr)

    # Redraw colorbar inside badge
    for b in 0:n_bars-1
        bx = cb_screen_x + (b / n_bars) * cb_screen_w
        bw = cb_screen_w / n_bars + 0.5
        bt = b / (n_bars - 1)
        r_c, g_c, b_c = colormap_turbo(bt)
        set_source_rgb(cr, r_c, g_c, b_c)
        Cairo.rectangle(cr, bx, cb_screen_y, bw, 8.0)
        Cairo.fill(cr)
    end
    set_source_rgb(cr, 0.2, 0.2, 0.3)
    Cairo.set_line_width(cr, 0.8)
    Cairo.rectangle(cr, cb_screen_x, cb_screen_y, cb_screen_w, 8.0)
    Cairo.stroke(cr)

    set_source_rgb(cr, 0.15, 0.18, 0.25)
    Cairo.select_font_face(cr, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_BOLD)
    Cairo.set_font_size(cr, 10.5)
    Cairo.move_to(cr, cb_screen_x, cb_screen_y - 5.0)
    Cairo.show_text(cr, @sprintf("Nearest Site Dist: %.3f (min) -> %.3f (max)", min_dist_V, max_dist_V))

    Cairo.show_page(cr)
    Cairo.finish(surf)

    # --- Generate Interleaved Text Summary Pages with LaTeX ---
    text_page1 = raw"""
\vspace*{30pt}
\section*{Figure 1: Voronoi Diagram of Set $T$}
\vspace{15pt}

\subsection*{Point Sets and Visual Elements}
\begin{itemize}
    \setlength{\itemsep}{12pt}
    \item \textbf{Site Set $T$ (Dark Navy Blue Dots):}
    The input collection of sites consisting of the origin $(0, 0)$ together with $11$ concentric rings of regular $10$-gons, each scaled by a factor of $1.1$.
    \item \textbf{Voronoi Cells (Light Gray-Blue Fills):}
    The Voronoi cell for each site in $T$, forming a symmetric, concentric flower-like pattern radiating outward from the origin.
    \item \textbf{Voronoi Edges (Gray Lines):}
    The bisecting boundaries separating adjacent Voronoi cells.
\end{itemize}
"""

    text_page2 = raw"""
\vspace*{30pt}
\section*{Figure 2: Origin Cell and Voronoi Vertices $U$}
\vspace{15pt}

\subsection*{Point Sets and Visual Elements}
\begin{itemize}
    \setlength{\itemsep}{12pt}
    \item \textbf{Origin Cell (Amber Fill and Boundary):}
    The central Voronoi cell belonging to the origin site $(0, 0)$.
    \item \textbf{Vertex Set $U$ (Bright Red Dots with Halos):}
    The $10$ Voronoi vertices of the cell of the origin. Each red vertex is equidistant to the origin and two adjacent sites on the innermost ring $S_0$.
    \item \textbf{Site Set $T$ (Dark Navy Blue Dots):}
    The background sites of the original set $T$.
\end{itemize}
"""

    text_page3 = raw"""
\vspace*{30pt}
\section*{Figure 3: Voronoi Diagram of $T \cup U$ with Closer New Vertices $V$}
\vspace{15pt}

\subsection*{Point Sets and Visual Elements}
\begin{itemize}
    \setlength{\itemsep}{12pt}
    \item \textbf{New Vertices $V$ (Royal Blue Dots with Halos):}
    The new Voronoi vertices formed after inserting the vertex set $U$ as new sites, restricted to the closer ring of vertices around the origin.
    \item \textbf{Original Sites $T$ (Dark Navy Blue Dots):}
    The original input points forming the concentric rings and origin.
    \item \textbf{Added Sites $U$ (Bright Red Dots):}
    The $10$ site locations corresponding to the Voronoi vertices $U$.
    \item \textbf{New Origin Cell (Cyan Fill):}
    The updated, reduced Voronoi cell of the origin in $\Vor(T \cup U)$.
\end{itemize}
"""

    text_page4 = raw"""
\vspace*{30pt}
\section*{Figure 4: Augmented Voronoi Diagram with New Cells Highlighted}
\vspace{15pt}

\subsection*{Point Sets and Visual Elements}
\begin{itemize}
    \setlength{\itemsep}{12pt}
    \item \textbf{New Voronoi Cells (Light Green Fills and Dark Green Borders):}
    The $10$ new Voronoi cells corresponding to the newly added sites $U$.
    \item \textbf{Added Sites $U$ (Bright Red Dots):}
    The $10$ site locations corresponding to the Voronoi vertices $U$.
    \item \textbf{Original Sites $T$ (Dark Navy Blue Dots):}
    The original input points forming the concentric rings and origin.
    \item \textbf{Voronoi Edges (Gray Lines):}
    The updated Voronoi tessellation of $T \cup U$.
\end{itemize}
"""

    text_page5 = raw"""
\vspace*{30pt}
\section*{Figure 5: Original Voronoi and Empty Disks at $V$}
\vspace{15pt}

\subsection*{Point Sets and Visual Elements}
\begin{itemize}
    \setlength{\itemsep}{12pt}
    \item \textbf{Vertex Set $V$ (Royal Blue Center Dots):}
    The set of new Voronoi vertices in $\Vor(T \cup U)$ not adjacent to the origin.
    \item \textbf{Original Sites $T$ (Dark Navy Blue Dots):}
    The original input points forming the concentric rings and origin.
    \item \textbf{Empty Disks (Semi-Transparent Colored Disks):}
    For each vertex $v \in V$, a disk centered at $v$ with radius equal to its distance to the nearest site in $T \cup U$.
    \item \textbf{Disk Colors (Blue to Red Spectrum):}
    The color of each disk designates the distance from its center to the nearest site.
    \item \textbf{Background Voronoi Diagram (Gray Lines):}
    The original Voronoi diagram $\Vor(T)$ overlaid for geometric reference.
\end{itemize}
"""

    latex_to_pdf([text_page1, text_page2, text_page3, text_page4, text_page5],
                 tmp_text_pdf; paperwidth=800, paperheight=800, margin=55)

    # --- Merge into Interleaved 10-Page PDF (Figure followed by Summary) ---
    pdf_merge(output_pdf,
              (tmp_diagrams_pdf, 1), (tmp_text_pdf, 1),
              (tmp_diagrams_pdf, 2), (tmp_text_pdf, 2),
              (tmp_diagrams_pdf, 3), (tmp_text_pdf, 3),
              (tmp_diagrams_pdf, 4), (tmp_text_pdf, 4),
              (tmp_diagrams_pdf, 5), (tmp_text_pdf, 5))

    println("Successfully generated: ", relpath(output_pdf, PROJROOT))
end

# CLI entry point
if abspath(PROGRAM_FILE) == @__FILE__
    out = length(ARGS) >= 1 ? ARGS[1] : "output/flower.pdf"
    generate_flower_pdf(out)
end

