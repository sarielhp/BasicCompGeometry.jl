#! /bin/env julial

# Comparison of two approx_diam_subspace implementations against exact naive version.
# 1. WSPD-based: approx_diameter_subspace (from src/Diameter.jl)
# 2. FST-based: MVBB.approx_diam_subspace (from src/MVBB.jl)
# 3. Naive: exact_diameter_subspace (from src/Diameter.jl)

using BasicCompGeometry
using BasicCompGeometry.MVBB
using LinearAlgebra
using Printf
using Random

"""
    generate_random_points(n, rng)
"""
function generate_random_points(n, rng)
    pts = Point{3, Float64}[]
    for _ in 1:n
        push!(pts, point(rand(rng)*1.0, rand(rng)*2.0, rand(rng)*4.0))
    end
    return PntSeq(pts)
end

function (@main)(args)
    num_tests = 100
    n = 100
    ε = 0.01
    
    @printf("Verifying approx_diam_subspace against Naive Exact (ε=%.3f, n=%d)\n\n", ε, n)

    # Prepare datasets and directions
    rng_data = Random.Xoshiro(123)
    datasets = []
    directions = []
    for i in 1:num_tests
        push!(datasets, generate_random_points(n, rng_data))
        push!(directions, normalize(Point{3, Float64}(randn(rng_data, 3))))
    end

    # Results
    res_exact = zeros(num_tests)
    res_wspd = zeros(num_tests)
    res_fst = zeros(num_tests)

    # 1. Exact
    for i in 1:num_tests
        res_exact[i] = exact_diameter_subspace(datasets[i], directions[i])
    end

    # 2. WSPD
    for i in 1:num_tests
        res_wspd[i] = approx_diameter_subspace(datasets[i], ε, directions[i])
    end

    # 3. FST
    for i in 1:num_tests
        res = MVBB.approx_diam_subspace(datasets[i], ε, directions[i])
        res_fst[i] = res.distance
    end

    # Analysis
    ratio_wspd = res_exact ./ res_wspd
    ratio_fst = res_exact ./ res_fst
    
    @printf("| %-10s | %-12s | %-12s | %-12s | %-12s |\n", 
            "Metric", "Min Ratio", "Max Ratio", "Avg Ratio", "Max % Error")
    @printf("|------------|--------------|--------------|--------------|--------------|\n")
    
    avg_wspd = sum(ratio_wspd)/num_tests
    max_err_wspd = (maximum(ratio_wspd) - 1.0) * 100.0
    @printf("| %-10s | %-12.6f | %-12.6f | %-12.6f | %-11.4f%% |\n",
            "WSPD", minimum(ratio_wspd), maximum(ratio_wspd), avg_wspd, max_err_wspd)
            
    avg_fst = sum(ratio_fst)/num_tests
    max_err_fst = (maximum(ratio_fst) - 1.0) * 100.0
    @printf("| %-10s | %-12.6f | %-12.6f | %-12.6f | %-11.4f%% |\n",
            "FST", minimum(ratio_fst), maximum(ratio_fst), avg_fst, max_err_fst)

    # Theoretical Bound Check: Exact / Approx should be in [1, 1+ε]
    if maximum(ratio_wspd) > (1.0 + ε) + 1e-12 || minimum(ratio_wspd) < 1.0 - 1e-12
        @printf("\nWARNING: WSPD ratio (Exact/WSPD) outside theoretical bounds [1, 1+ε]!\n")
        @printf("Max Ratio: %.12f, Min Ratio: %.12f\n", maximum(ratio_wspd), minimum(ratio_wspd))
    end
    if maximum(ratio_fst) > (1.0 + ε) + 1e-12 || minimum(ratio_fst) < 1.0 - 1e-12
        @printf("\nWARNING: FST ratio (Exact/FST) outside theoretical bounds [1, 1+ε]!\n")
        @printf("Max Ratio: %.12f, Min Ratio: %.12f\n", maximum(ratio_fst), minimum(ratio_fst))
    end

    return 0
end
