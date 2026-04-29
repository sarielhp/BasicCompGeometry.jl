#! /bin/env julial

# High-scale comparison of subspace diameter implementations:
# 1. WSPD-based (Standard)
# 2. FST-based (MVBB)
# 3. Naive Exact (only for n <= 10000)

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

function run_benchmark(n, num_tests, ε, run_exact=true)
    rng_data = Random.Xoshiro(123)
    datasets = []
    directions = []
    for i in 1:num_tests
        push!(datasets, generate_random_points(n, rng_data))
        push!(directions, normalize(Point{3, Float64}(randn(rng_data, 3))))
    end

    # 1. Exact Naive
    results_exact = zeros(num_tests)
    t_exact = 0.0
    if run_exact
        t_start_exact = time()
        for i in 1:num_tests
            results_exact[i] = exact_diameter_subspace(datasets[i], directions[i])
        end
        t_exact = (time() - t_start_exact) / num_tests
    end

    # 2. WSPD Standard
    results_wspd = zeros(num_tests)
    t_start_wspd = time()
    for i in 1:num_tests
        results_wspd[i] = approx_diameter_subspace(datasets[i], ε, directions[i])
    end
    t_wspd = (time() - t_start_wspd) / num_tests

    # 3. FST
    results_fst = zeros(num_tests)
    t_start_fst = time()
    for i in 1:num_tests
        res = MVBB.approx_diam_subspace(datasets[i], ε, directions[i])
        results_fst[i] = res.distance
    end
    t_fst = (time() - t_start_fst) / num_tests

    # Ratios (Exact / Approx) or (FST / WSPD if no exact)
    max_r_wspd = run_exact ? maximum(results_exact ./ results_wspd) : 0.0
    max_r_fst = run_exact ? maximum(results_exact ./ results_fst) : 0.0
    max_diff_wf = maximum(abs.(results_wspd .- results_fst))
    
    return (
        n = n,
        ε = ε,
        t_exact = t_exact,
        t_wspd = t_wspd,
        t_fst = t_fst,
        max_r_wspd = max_r_wspd,
        max_r_fst = max_r_fst,
        max_diff_wf = max_diff_wf
    )
end

function (@main)(args)
    num_tests = 100
    
    @printf("Benchmark: High-Scale Subspace Diameter Comparison\n")
    @printf("Settings: num_tests=%d\n\n", num_tests)
    
    # Warm-up
    warm_rng = Random.Xoshiro(42)
    warm_ps = generate_random_points(100, warm_rng)
    warm_dir = normalize(Point{3, Float64}(randn(warm_rng, 3)))
    _ = exact_diameter_subspace(warm_ps, warm_dir)
    _ = approx_diameter_subspace(warm_ps, 0.01, warm_dir)
    _ = MVBB.approx_diam_subspace(warm_ps, 0.01, warm_dir)

    # Table 1: Standard sizes
    @printf("Table 1: Standard Sizes (ε=0.010)\n")
    sizes = [100, 1000, 10000]
    results1 = []
    for n in sizes
        @printf("  Running n=%d...\n", n)
        push!(results1, run_benchmark(n, num_tests, 0.010, true))
    end

    @printf("\n| %-7s | %-12s | %-12s | %-12s | %-8s | %-8s | %-8s |\n", 
            "n", "Avg Naive(s)", "Avg WSPD(s)", "Avg FST(s)", "Max R W", "Max R F", "Speedup")
    @printf("|---------|--------------|--------------|--------------|----------|----------|----------|\n")
    for r in results1
        @printf("| %-7d | %-12.6f | %-12.6f | %-12.6f | %-8.4f | %-8.4f | %-8.2fx |\n",
                r.n, r.t_exact, r.t_wspd, r.t_fst, r.max_r_wspd, r.max_r_fst, r.t_wspd / r.t_fst)
    end

    # Table 2: High Scale (n=200,000)
    @printf("\nTable 2: High Scale (n=200,000, No Naive)\n")
    eps_values = [0.0, 0.1, 1.0]
    results2 = []
    for ε in eps_values
        @printf("  Running ε=%.1f...\n", ε)
        push!(results2, run_benchmark(200000, num_tests, ε, false))
    end

    @printf("\n| %-7s | %-8s | %-12s | %-12s | %-12s | %-8s |\n", 
            "n", "ε", "Avg WSPD(s)", "Avg FST(s)", "Max Diff WF", "Speedup")
    @printf("|---------|----------|--------------|--------------|--------------|----------|\n")
    for r in results2
        @printf("| %-7d | %-8.1f | %-12.6f | %-12.6f | %-12.4e | %-8.2fx |\n",
                r.n, r.ε, r.t_wspd, r.t_fst, r.max_diff_wf, r.t_wspd / r.t_fst)
    end

    return 0
end
