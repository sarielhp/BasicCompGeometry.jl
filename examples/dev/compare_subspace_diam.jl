#! /bin/env julial

# Comparison of two approx_diam_subspace implementations:
# 1. WSPD-based: approx_diameter_subspace (from src/Diameter.jl)
# 2. FST-based: MVBB.approx_diam_subspace (from src/MVBB.jl)
#
# This script performs 100 tests with n=10000 points each.
# Algorithms are run sequentially (all tests for alg A, then all for alg B)
# to minimize interference and ordering effects.

using BasicCompGeometry
using BasicCompGeometry.MVBB
using LinearAlgebra
using Printf
using Random

"""
    generate_random_points(n, rng)

Generate n random 3D points in [0,1]x[0,2]x[0,4].
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
    
    @printf("Comparing two versions of approx_diam_subspace (Ordered Execution)\n")
    @printf("Settings: n=%d points, %d tests, ε=%.3f\n\n", n, num_tests, ε)

    # 1. JIT Warm-up (without measuring time)
    @printf("Warming up JIT compiler...\n")
    rng_warm = Random.Xoshiro(42)
    warmup_ps = generate_random_points(100, rng_warm)
    warmup_dir = normalize(Point{3, Float64}(randn(rng_warm, 3)))
    
    _ = approx_diameter_subspace(warmup_ps, ε, warmup_dir)
    _ = MVBB.approx_diam_subspace(warmup_ps, ε, warmup_dir)
    @printf("Warm-up complete.\n\n")

    # Prepare datasets and directions beforehand for consistency
    datasets = []
    directions = []
    rng_data = Random.Xoshiro(123)
    for i in 1:num_tests
        push!(datasets, generate_random_points(n, rng_data))
        push!(directions, normalize(Point{3, Float64}(randn(rng_data, 3))))
    end

    # 2. Run WSPD-based tests
    @printf("Running %d tests for WSPD-based version...\n", num_tests)
    results_wspd = zeros(num_tests)
    t_start_wspd = time()
    for i in 1:num_tests
        results_wspd[i] = approx_diameter_subspace(datasets[i], ε, directions[i])
        if i % 20 == 0
            @printf("  Completed %d/%d...\n", i, num_tests)
        end
    end
    t_total_wspd = time() - t_start_wspd

    # 3. Run FST-based (MVBB) tests
    @printf("Running %d tests for FST-based version...\n", num_tests)
    results_fst = zeros(num_tests)
    t_start_fst = time()
    for i in 1:num_tests
        res = MVBB.approx_diam_subspace(datasets[i], ε, directions[i])
        results_fst[i] = res.distance
        if i % 20 == 0
            @printf("  Completed %d/%d...\n", i, num_tests)
        end
    end
    t_total_fst = time() - t_start_fst

    # 4. Analysis
    ratios = results_wspd ./ results_fst
    min_ratio = minimum(ratios)
    max_ratio = maximum(ratios)
    max_diff = maximum(abs.(results_wspd .- results_fst))

    @printf("\nRESULTS:\n")
    @printf("--------------------------------------------------\n")
    @printf("Maximum Difference:          %.10e\n", max_diff)
    @printf("Min Ratio (WSPD/FST):        %.10f\n", min_ratio)
    @printf("Max Ratio (WSPD/FST):        %.10f\n", max_ratio)
    @printf("Total Time (WSPD):           %.4f s\n", t_total_wspd)
    @printf("Total Time (FST):            %.4f s\n", t_total_fst)
    @printf("Avg Time per test (WSPD):    %.6f s\n", t_total_wspd / num_tests)
    @printf("Avg Time per test (FST):     %.6f s\n", t_total_fst / num_tests)
    @printf("Speedup (WSPD/FST):          %.2fx\n", t_total_fst / t_total_wspd)
    @printf("--------------------------------------------------\n")

    return 0
end
