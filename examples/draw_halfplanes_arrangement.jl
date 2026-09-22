#!/usr/bin/env julia

using QuickEnv

# Example: draw_halfplanes_arrangement.jl
# Reads a saved set of halfplanes (or generates random ones), computes their arrangement
# and default view, computes the faces inside the view via vertical trapezoid chains,
# and outputs a 2-page publication-quality PDF:
#   Page 1: Polygonal faces filled with distinct random colors, with boundary lines.
#   Page 2: Arrangement depth heatmap with inward directional whiskers along boundary lines.

using BasicCompGeometry
using Cairo
import BasicCompGeometry: clip, width, height
using LinearAlgebra
using Printf
using Random

function hsv_to_rgb(h::Real, s::Real, v::Real)
    h_norm = mod(h, 360.0) / 60.0
    c = v * s
    x = c * (1.0 - abs(mod(h_norm, 2.0) - 1.0))
    m = v - c
    r, g, b = 0.0, 0.0, 0.0
    if 0.0 <= h_norm < 1.0
        r, g, b = c, x, 0.0
    elseif 1.0 <= h_norm < 2.0
        r, g, b = x, c, 0.0
    elseif 2.0 <= h_norm < 3.0
        r, g, b = 0.0, c, x
    elseif 3.0 <= h_norm < 4.0
        r, g, b = 0.0, x, c
    elseif 4.0 <= h_norm < 5.0
        r, g, b = x, 0.0, c
    else
        r, g, b = c, 0.0, x
    end
    return (r + m, g + m, b + m)
end

function distinct_face_color(i::Int)
    golden_ratio = 0.618033988749895
    h = mod(i * golden_ratio, 1.0) * 360.0
    s = 0.45 + 0.20 * (mod(i, 3) / 2.0)
    v = 0.90 + 0.08 * mod(i, 2)
    return hsv_to_rgb(h, s, v)
end

function viridis_color(t::Real)
    t_clamped = clamp(Float64(t), 0.0, 1.0)
    pts = [
        (0.00, 0.267, 0.005, 0.329),
        (0.25, 0.231, 0.322, 0.545),
        (0.50, 0.129, 0.569, 0.551),
        (0.75, 0.369, 0.789, 0.383),
        (1.00, 0.993, 0.906, 0.144)
    ]
    for i in 1:(length(pts)-1)
        t0, r0, g0, b0 = pts[i]
        t1, r1, g1, b1 = pts[i+1]
        if t_clamped <= t1
            u = (t_clamped - t0) / (t1 - t0)
            return (r0 + u * (r1 - r0), g0 + u * (g1 - g0), b0 + u * (b1 - b0))
        end
    end
    return (0.993, 0.906, 0.144)
end

function parse_args()
    input_file = nothing
    output_file = joinpath(@__DIR__, "..", "output", "halfplanes_arrangement.pdf")
    n_rand = 12
    cw = 800.0
    ch = 800.0
    factor = 1.08

    args = ARGS
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--input" || arg == "-i"
            i += 1
            input_file = args[i]
        elseif arg == "--output" || arg == "-o"
            i += 1
            output_file = args[i]
        elseif arg == "--n" || arg == "-n"
            i += 1
            n_rand = parse(Int, args[i])
        elseif arg == "--width"
            i += 1
            cw = parse(Float64, args[i])
        elseif arg == "--height"
            i += 1
            ch = parse(Float64, args[i])
        elseif arg == "--factor"
            i += 1
            factor = parse(Float64, args[i])
        elseif !startswith(arg, "-") && input_file === nothing
            if isfile(arg)
                input_file = arg
            elseif tryparse(Int, arg) !== nothing
                n_rand = parse(Int, arg)
            else
                input_file = arg
            end
        end
        i += 1
    end
    return (input_file = input_file, output_file = output_file, n_rand = n_rand,
            cw = cw, ch = ch, factor = factor)
end

function main()
    opts = parse_args()
    mkpath(dirname(abspath(opts.output_file)))

    # Step 1: Load or generate halfplanes
    hps = Halfplane2F[]
    source_desc = ""
    if opts.input_file !== nothing && isfile(opts.input_file)
        hps = read_halfplanes(opts.input_file)
        source_desc = "Loaded $(length(hps)) halfplanes from $(opts.input_file)"
    else
        default_candidate = joinpath(@__DIR__, "..", "output", "pruned_halfplanes_k5.txt")
        if isfile(default_candidate)
            hps = read_halfplanes(default_candidate)
            source_desc = "Loaded $(length(hps)) halfplanes from $(default_candidate)"
        else
            Random.seed!(42)
            hps = [rand_halfplane() for _ in 1:opts.n_rand]
            source_desc = "Generated $(length(hps)) random halfplanes (seed=42)"
        end
    end

    println("=========================================================")
    println("  Arrangement Visualizer: Faces & Depth Heatmap")
    println("=========================================================")
    println("  Source       : ", source_desc)
    println("  Halfplanes   : ", length(hps))

    # Step 2: Compute default view and vertical decomposition
    view = default_view(hps; factor = opts.factor)
    println("  Default View : ", view)

    decomp = vertical_decomposition(hps, view)
    traps = decomp.trapezoids
    println("  Trapezoids   : ", length(traps))

    # Step 3: Compute trapezoid chains and polygonal faces
    chains = trapezoid_chains(decomp)
    face_polys = faces(decomp)
    num_faces = length(face_polys)
    println("  Chains/Faces : ", num_faces)

    # Step 4: Compute centroid and depth of each face
    centroids = Point2F[sum(f.pnts) / length(f.pnts) for f in face_polys]
    face_depths = [depth(hps, c) for c in centroids]
    d_min, d_max = isempty(face_depths) ? (0, 0) : extrema(face_depths)
    println("  Depth range  : min = $d_min, max = $d_max")

    margin = 55.0
    cw = opts.cw
    ch = opts.ch

    open_canvas(opts.output_file, cw, ch) do canvas
        # -------------------------------------------------------------
        # PAGE 1: Arrangement Faces with Distinct Colors & Lines
        # -------------------------------------------------------------
        cairo_draw_setup(canvas, view, cw, ch, margin)

        # 1. Fill each face with a distinct pleasant color
        for i in 1:num_faces
            r, g, b = distinct_face_color(i)
            Cairo.set_source_rgb(canvas, r, g, b)
            cairo_draw_polygon(canvas, face_polys[i]; fill=true, stroke=true, line_width=1.0)
        end

        # 2. Draw input halfplane boundary lines across the view
        Cairo.set_source_rgb(canvas, 0.1, 0.1, 0.1)
        cairo_set_line_width(canvas, 2.0)
        for h in hps
            seg = clip(h.boundary, view)
            if seg !== nothing
                Cairo.move_to(canvas, seg.p[1], seg.p[2])
                Cairo.line_to(canvas, seg.q[1], seg.q[2])
            end
        end
        Cairo.stroke(canvas)

        # 3. Draw view bounding box outline
        Cairo.set_source_rgb(canvas, 0.15, 0.15, 0.15)
        cairo_set_line_width(canvas, 2.5)
        Cairo.rectangle(canvas, view.mini[1], view.mini[2], width(view), height(view))
        Cairo.stroke(canvas)

        # 4. Header title in device pixel coordinates
        Cairo.save(canvas)
        Cairo.reset_transform(canvas)
        Cairo.select_font_face(canvas, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_BOLD)
        Cairo.set_font_size(canvas, 16.0)
        Cairo.set_source_rgb(canvas, 0.15, 0.15, 0.2)
        Cairo.move_to(canvas, 40.0, 32.0)
        Cairo.show_text(canvas, @sprintf("Arrangement Faces: %d faces from %d trapezoids (%d halfplanes)", num_faces, length(traps), length(hps)))
        Cairo.restore(canvas)

        description(canvas, "Page 1: $(num_faces) polygonal faces of $(length(hps)) halfplanes, formed from $(length(chains)) trapezoid chains.")

        Cairo.show_page(canvas)

        # -------------------------------------------------------------
        # PAGE 2: Depth Heatmap with Whiskers
        # -------------------------------------------------------------
        cairo_draw_setup(canvas, view, cw, ch, margin)

        # 1. Fill each face according to its arrangement depth
        for i in 1:num_faces
            d = face_depths[i]
            t = d_max == d_min ? 0.5 : (d - d_min) / (d_max - d_min)
            r, g, b = viridis_color(t)
            Cairo.set_source_rgb(canvas, r, g, b)
            cairo_draw_polygon(canvas, face_polys[i]; fill=true, stroke=true, line_width=0.8)
        end

        # 2. Draw halfplane boundary lines with inward whiskers
        cairo_draw_halfplanes(
            canvas, hps, view;
            line_width = 1.8,
            tick_len = 8.0,
            tick_spacing = 24.0,
            color = (0.1, 0.1, 0.1)
        )

        # 3. Draw view bounding box outline
        Cairo.set_source_rgb(canvas, 0.15, 0.15, 0.15)
        cairo_set_line_width(canvas, 2.5)
        Cairo.rectangle(canvas, view.mini[1], view.mini[2], width(view), height(view))
        Cairo.stroke(canvas)

        # 4. Header title and colorbar in device pixel coordinates
        Cairo.save(canvas)
        Cairo.reset_transform(canvas)
        Cairo.select_font_face(canvas, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_BOLD)
        Cairo.set_font_size(canvas, 16.0)
        Cairo.set_source_rgb(canvas, 0.15, 0.15, 0.2)
        Cairo.move_to(canvas, 40.0, 32.0)
        Cairo.show_text(canvas, @sprintf("Arrangement Depth Heatmap (depths %d .. %d) with Inward Whiskers", d_min, d_max))

        # Colorbar
        cb_x = cw - 320.0
        cb_y = ch - 26.0
        cb_w = 260.0
        cb_h = 12.0
        steps = 60
        for s in 0:(steps-1)
            t_col = s / (steps - 1)
            r, g, b = viridis_color(t_col)
            Cairo.set_source_rgb(canvas, r, g, b)
            Cairo.rectangle(canvas, cb_x + s * (cb_w / steps), cb_y, cb_w / steps + 0.5, cb_h)
            Cairo.fill(canvas)
        end
        Cairo.set_source_rgb(canvas, 0.2, 0.2, 0.2)
        Cairo.set_line_width(canvas, 1.0)
        Cairo.rectangle(canvas, cb_x, cb_y, cb_w, cb_h)
        Cairo.stroke(canvas)

        # Colorbar labels
        Cairo.set_font_size(canvas, 11.0)
        Cairo.set_source_rgb(canvas, 0.3, 0.3, 0.3)
        Cairo.move_to(canvas, cb_x - 70.0, cb_y + 10.0)
        Cairo.show_text(canvas, @sprintf("min = %d", d_min))
        Cairo.move_to(canvas, cb_x + cb_w + 8.0, cb_y + 10.0)
        Cairo.show_text(canvas, @sprintf("max = %d", d_max))
        Cairo.restore(canvas)

        description(canvas, "Page 2: Arrangement depth heatmap (depths $d_min..$d_max) with inward halfplane whiskers.")
    end

    println("  Output PDF   : ", opts.output_file)
    println("=========================================================")
    println("  Done.")
end

main()
