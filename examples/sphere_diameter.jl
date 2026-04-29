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
using Cairo
using Colors

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
    println("\nGenerating BBT visualization...")
    P_viz = Polygon_random_sphere(2, Float64, 200)
    translate!(P_viz, point(400.0, 400.0))
    scale!(P_viz, 300.0)
    tree = BBT.Tree_init(P_viz)
    BBT.Tree_fully_expand(tree)
    mkpath("output")
    BBT.Tree_draw(tree, joinpath("output", "sphere_bbt.pdf"))

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
