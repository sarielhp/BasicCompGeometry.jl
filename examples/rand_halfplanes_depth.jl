#!/usr/bin/env julia

using QuickEnv

#
# examples/rand_halfplanes_depth.jl
#
# Randomly picks n (default 100) random halfplanes, computes their default view,
# and then computes the vertical decomposition inside the view.
# Next, for each vertical trapezoid, computes its center (average of its vertices),
# and evaluates its depth (number of halfplanes containing the center).
# The total depth of the arrangement is the minimum depth of any of these centers.
# Outputs this number along with arrangement statistics.
#

using BasicCompGeometry
using LinearAlgebra
using Printf
using Random
using Statistics

function run_rand_halfplanes_depth(n::Int = 100; seed::Int = 42, do_plot::Bool = false, plot_path::String = "output/rand_halfplanes_depth.pdf")
    rng = Random.Xoshiro(seed)

    println("="^70)
    @printf("Program: rand_halfplanes_depth\n")
    @printf("Random Halfplanes Arrangement & Vertical Decomposition Depth\n")
    println("="^70)
    @printf("Number of halfplanes (n):       %d\n", n)
    @printf("Random seed:                    %d\n\n", seed)

    # 1. Randomly pick n halfplanes
    println("[1/4] Generating random halfplanes...")
    t_gen = @elapsed hps = [rand_halfplane(rng) for _ = 1:n]
    @printf("      Generated %d halfplanes in %.3f s\n\n", n, t_gen)

    # 2. Compute default view
    println("[2/4] Computing default view rectangle...")
    t_view = @elapsed view = default_view(hps; factor = 1.1)
    view_w = width(view)
    view_h = height(view)
    view_area = view_w * view_h
    @printf("      View box:                 [%.4f, %.4f] × [%.4f, %.4f]\n",
            view.mini[1], view.maxi[1], view.mini[2], view.maxi[2])
    @printf("      Dimensions:               width = %.4f, height = %.4f, area = %.4f\n\n",
            view_w, view_h, view_area)

    # 3. Compute vertical decomposition inside the view
    println("[3/4] Computing vertical decomposition inside view...")
    t_decomp = @elapsed decomp = vertical_decomposition(hps, view)
    num_traps = length(decomp)
    num_verts = length(vertices(decomp))
    total_area = area(decomp)
    rel_area_err = abs(total_area - view_area) / view_area
    @printf("      Decomposition time:       %.3f s\n", t_decomp)
    @printf("      Vertices in view:         %d\n", num_verts)
    @printf("      Vertical trapezoids:      %d\n", num_traps)
    @printf("      Area conservation check:  sum(areas) = %.6f (relative error: %.2e)\n\n",
            total_area, rel_area_err)

    # 4. Compute center of each vertical trapezoid and its depth
    println("[4/4] Computing trapezoid centers and evaluating depths...")
    centers = [centroid(trap.corners) for trap in decomp]

    t_depth = @elapsed depths = [depth(hps, c) for c in centers]
    @printf("      Depth query time:         %.3f s (%.1f µs/trapezoid)\n\n",
            t_depth, (t_depth / num_traps) * 1e6)

    # 5. Summary statistics and total depth
    min_depth = minimum(depths)
    max_depth = maximum(depths)
    mean_depth = mean(depths)
    med_depth = median(depths)

    println("="^70)
    @printf("ARRANGEMENT DEPTH SUMMARY\n")
    println("="^70)
    @printf("  • Total depth of arrangement (min depth): %d\n", min_depth)
    @printf("  • Maximum depth among centers:            %d\n", max_depth)
    @printf("  • Median depth:                           %.1f\n", med_depth)
    @printf("  • Mean depth:                             %.2f\n", mean_depth)
    println("="^70)

    # Clear direct output of the requested number
    @printf("\nTotal depth of the arrangement: %d\n", min_depth)

    # Optional plot
    if do_plot
        render_plot(decomp, depths, min_depth, max_depth, plot_path)
    end

    return min_depth
end

function render_plot(decomp, depths, min_d, max_d, outfile)
    if Base.find_package("Cairo") === nothing
        println("Note: Cairo not available in environment; skipping plot rendering.")
        return
    end

    # Use Cairo via extension
    cairo_mod = Base.get_extension(BasicCompGeometry, :CairoExt)
    if cairo_mod === nothing
        # Try loading Cairo if possible
        try
            @eval using Cairo
        catch
            println("Note: Could not activate CairoExt; skipping plot.")
            return
        end
    end

    mkpath(dirname(outfile))
    view = decomp.view
    cw = 900
    ch = 900
    c = Canvas(outfile, cw, ch; title = "Random Halfplanes Arrangement Depth")
    cairo_draw_setup(c, view, cw, ch, 25)

    span = max(1, max_d - min_d)
    for (i, trap) in enumerate(decomp.trapezoids)
        d = depths[i]
        frac = clamp((d - min_d) / span, 0.0, 1.0)
        # Heatmap colormap: blue (low depth) -> green -> red (high depth)
        r = frac < 0.5 ? 2.0 * frac * 0.7 : 0.7 + 0.3 * (frac - 0.5) * 2.0
        g = frac < 0.5 ? 0.2 + 1.2 * frac : 0.8 * (1.0 - (frac - 0.5) * 2.0)
        b = frac < 0.5 ? 0.9 * (1.0 - 2.0 * frac) : 0.1
        Cairo.set_source_rgba(c, r, g, b, 0.55)

        poly = trap.corners
        Cairo.move_to(c, poly[1][1], poly[1][2])
        for j = 2:4
            Cairo.line_to(c, poly[j][1], poly[j][2])
        end
        Cairo.close_path(c)
        Cairo.fill(c)

        # Trapezoid border
        Cairo.set_source_rgba(c, 0.2, 0.2, 0.2, 0.3)
        cairo_set_line_width(c, 0.5)
        for j = 1:4
            j_next = j == 4 ? 1 : j + 1
            Cairo.move_to(c, poly[j][1], poly[j][2])
            Cairo.line_to(c, poly[j_next][1], poly[j_next][2])
            Cairo.stroke(c)
        end
    end

    # Stroke halfplane boundary lines
    Cairo.set_source_rgba(c, 0.05, 0.05, 0.1, 0.75)
    cairo_set_line_width(c, 1.2)
    for l in decomp.lines
        if abs(l.u[1]) > 1e-12
            x1 = view.mini[1]
            y1 = l.p[2] + (l.u[2] / l.u[1]) * (x1 - l.p[1])
            x2 = view.maxi[1]
            y2 = l.p[2] + (l.u[2] / l.u[1]) * (x2 - l.p[1])
            Cairo.move_to(c, x1, y1)
            Cairo.line_to(c, x2, y2)
            Cairo.stroke(c)
        end
    end

    description(c, @sprintf("Vertical decomposition of %d halfplanes inside default view. Total arrangement depth = %d.",
                            length(decomp.lines), min_d))
    Cairo.finish(c)
    println("Saved visualization to: ", outfile)
end

function main()
    # CLI args: [n] [seed] [--plot]
    n = 100
    seed = 42
    do_plot = false

    args = String[]
    for a in ARGS
        if a == "--plot" || a == "-p"
            do_plot = true
        else
            push!(args, a)
        end
    end

    if length(args) >= 1
        n = parse(Int, args[1])
    end
    if length(args) >= 2
        seed = parse(Int, args[2])
    end

    run_rand_halfplanes_depth(n; seed = seed, do_plot = do_plot)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
