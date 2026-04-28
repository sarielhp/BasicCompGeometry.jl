#! /usr/bin/env julia
#
# examples/sphere_diameter.jl
#
# Compute the diameter of a random point set on the sphere of various dimensions.
# This example demonstrates the use of BoundingBoxTrees (BBT), Well-Separated 
# Pairs Decomposition (WSPD), and diameter approximation algorithms.
#

using BasicCompGeometry
using Printf
using DataFrames
using PrettyTables
using LinearAlgebra

# Visualization requires Cairo and Colors
const HAS_VIS = try
    import Cairo, Colors
    true
catch e
    false
end

if !HAS_VIS
    @warn "Cairo and/or Colors not found. Visualization disabled. Install them to enable PDF output."
end

"""
    BBTree_draw(tree, filename)

A visualization utility for BoundingBoxTrees. 
Draws the tree levels into a multi-page PDF.
"""
function BBTree_draw(tree::BBT.Tree{2,T,S}, filename::String) where {T,S}
    if !HAS_VIS
        println("BBTree_draw requires Cairo and Colors packages.")
        return
    end

    # Import necessary symbols from Cairo and Colors
    # Using 'using' inside function requires a specific trick or just 'import'
    # Here we just use the module prefix to be safe and clear.
    C = Cairo
    Clr = Colors

    # Utility to draw a single bounding box
    function bbox_draw(context, bb, color)
        C.set_source(context, color)
        bl = bottom_left(bb)
        w = width(bb, 1)
        h = height(bb, 2)
        C.rectangle(context, bl[1], bl[2], w, h)
        C.fill(context)
    end

    # Recursive utility to draw nodes at specific levels
    function node_draw(context, node, level, range)
        node.f_leaf && return
        if level ∈ range
            yellow_transparent = Clr.coloralpha(Clr.parse(Clr.Colorant, "yellow"), 0.1)
            C.set_source(context, yellow_transparent)

            # Expand slightly for better visibility
            bb = node.bb + (diam(node.bb) * 0.01)
            bl = bottom_left(bb)
            w = width(bb, 1)
            h = height(bb, 2)

            C.rectangle(context, bl[1], bl[2], w, h)
            C.fill_preserve(context)

            C.set_source(context, Clr.parse(Clr.Colorant, "blue"))
            C.set_line_width(context, 1.0)
            C.stroke(context)
        end
        if !isnothing(node.left)
            node_draw(context, node.left, level + 1, range)
        end
        if !isnothing(node.right)
            node_draw(context, node.right, level + 1, range)
        end
    end

    # Main drawing logic
    # Expand root BB for margins
    bb_root = tree.root.bb + (diam(tree.root.bb) * 0.2)

    # We use the max coordinates for surface size
    tr = top_right(bb_root)
    surface = C.CairoPDFSurface(filename, tr[1], tr[2])
    context = C.CairoContext(surface)

    # Page 1: Overview
    bbox_draw(context, bb_root, Clr.parse(Clr.Colorant, "lightblue"))
    d = BBT.depth(tree.root)
    node_draw(context, tree.root, 0, 0:d)
    C.show_page(context)

    # Subsequent pages: Level-by-level visualization
    for i = (d-1):-1:0
        bbox_draw(context, bb_root, Clr.parse(Clr.Colorant, "lightblue"))
        node_draw(context, tree.root, 0, i:(i+1))
        C.show_page(context)
    end

    C.finish(surface)
    println("Created tree visualization: $filename")
end

function test_diameter(P, N::Int, D::Int)
    # If 2D, scale and shift for visualization purposes
    if D == 2
        translate!(P, point(200.0, 200.0))
        scale!(P, 500.0)
    end

    t_approx = @timed approx_diam = approx_diameter(P, 0.1)
    t_exact = @timed exact_diam = exact_diameter(P)

    @printf(
        "D: %d  N: %8d  Approx: %10.6f  Exact: %10.6f    T_approx: %10.4fs  T_exact: %10.4fs\n",
        D,
        N,
        approx_diam,
        exact_diam,
        t_approx.time,
        t_exact.time
    )

    return DataFrame(
        dimension = D,
        N = N,
        approx_runtime = @sprintf("%.6f", t_approx.time),
        exact_runtime = @sprintf("%.6f", t_exact.time),
        approx_ratio = @sprintf("%.4f", approx_diam/exact_diam)
    )
end

function diam_test_sphere(D, iters, filename)
    # Warmup
    P_warm = Polygon_random_sphere(D, Float64, 20)
    approx_diameter(P_warm, 0.1)
    exact_diameter(P_warm)

    df = DataFrame()
    for i = 1:iters
        N = 2^i
        P = Polygon_random_sphere(D, Float64, N)
        new_row = test_diameter(P, N, D)
        append!(df, new_row)
    end

    mkpath(dirname(filename))
    open(filename, "w") do fl
        # Rename columns to match desired headers if using DataFrame
        df_display = rename(
            df,
            :dimension => "Dim",
            :N => "N",
            :approx_runtime => "RT Approx",
            :exact_runtime => "RT Exact",
            :approx_ratio => "Ratio",
        )
        pretty_table(fl, df_display, alignment = :r, backend = :markdown)
    end
    println("Results saved to: $filename")
end

function main()
    println("--- BasicCompGeometry Example: Sphere Diameter ---")

    # 1. BBT Visualization Example (2D)
    if HAS_VIS
        println("\nGenerating BBT visualization...")
        P_viz = Polygon_random_sphere(2, Float64, 200)
        translate!(P_viz, point(400.0, 400.0))
        scale!(P_viz, 300.0)
        tree = BBT.Tree_init(P_viz)
        BBT.Tree_fully_expand(tree)
        BBTree_draw(tree, "sphere_bbt.pdf")
    end

    # 2. Performance Comparison
    println("\nRunning diameter tests (Approx vs Exact)...")
    # Using smaller iters for the example run
    diam_test_sphere(2, 12, "results/example_2d.md")
    diam_test_sphere(3, 10, "results/example_3d.md")

    println("\nExample completed successfully.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
