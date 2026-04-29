#! /bin/env julial

# Multi-scale verification of approx_diam_subspace implementations:
# 1. WSPD-based (Diameter.jl)
# 2. FST-based (MVBB.jl)
# Compared against Naive O(n^2) Exact (Diameter.jl)

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

function run_verify(n, num_tests, ε)
    rng_data = Random.Xoshiro(123)
    datasets = []
    directions = []
    for i in 1:num_tests
        push!(datasets, generate_random_points(n, rng_data))
        push!(directions, normalize(Point{3, Float64}(randn(rng_data, 3))))
    end

    ratios_wspd = zeros(num_tests)
    ratios_fst = zeros(num_tests)
    t_exact = 0.0
    t_wspd = 0.0
    t_fst = 0.0

    for i in 1:num_tests
        # Exact
        t0 = time()
        exact = exact_diameter_subspace(datasets[i], directions[i])
        t_exact += (time() - t0)

        # WSPD
        t0 = time()
        wspd = approx_diameter_subspace(datasets[i], ε, directions[i])
        t_wspd += (time() - t0)

        # FST
        t0 = time()
        res = MVBB.approx_diam_subspace(datasets[i], ε, directions[i])
        fst = res.distance
        t_fst += (time() - t0)

        ratios_wspd[i] = exact / wspd
        ratios_fst[i] = exact / fst
    end

    return (
        n = n,
        avg_t_exact = t_exact / num_tests,
        avg_t_wspd = t_wspd / num_tests,
        avg_t_fst = t_fst / num_tests,
        max_ratio_wspd = maximum(ratios_wspd),
        max_ratio_fst = maximum(ratios_fst)
    )
end

function (@main)(args)
    ε = 0.01
    
    # Warm-up
    warm_rng = Random.Xoshiro(42)
    warm_ps = generate_random_points(100, warm_rng)
    warm_dir = normalize(Point{3, Float64}(randn(warm_rng, 3)))
    _ = approx_diameter_subspace(warm_ps, ε, warm_dir)
    _ = MVBB.approx_diam_subspace(warm_ps, ε, warm_dir)
    _ = exact_diameter_subspace(warm_ps, warm_dir)

    # (n, num_tests)
    configs = [(100, 100), (1000, 50), (10000, 5)]

    @printf("Multi-Scale Verification (ε=%.3f)\n", ε)
    @printf("| %-7s | %-12s | %-12s | %-12s | %-10s | %-10s |\n", 
            "n", "Avg Naive(s)", "Avg WSPD(s)", "Avg FST(s)", "Max R WSPD", "Max R FST")
    @printf("|---------|--------------|--------------|--------------|------------|------------|\n")

    for (n, nt) in configs
        res = run_verify(n, nt, ε)
        @printf("| %-7d | %-12.6f | %-12.6f | %-12.6f | %-10.4f | %-10.4f |\n",
                res.n, res.avg_t_exact, res.avg_t_wspd, res.avg_t_fst, res.max_ratio_wspd, res.max_ratio_fst)
    end

    return 0
end
