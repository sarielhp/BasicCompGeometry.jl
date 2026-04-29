#! /bin/env julial

using BasicCompGeometry
using BasicCompGeometry.MVBB
using LinearAlgebra
using Printf
using Random

"""
    run_comparison_test(n, ε, rng)

Run a single comparison test and return (d_wspd, d_fst, t_wspd, t_fst).
"""
function run_comparison_test(n, ε, rng)
    # Generate points
    pts = Point{3, Float64}[]
    for _ in 1:n
        push!(pts, point(rand(rng)*1.0, rand(rng)*2.0, rand(rng)*4.0))
    end
    ps = PntSeq(pts)

    # Random unit direction
    dir = normalize(Point{3, Float64}(randn(rng, 3)))

    # Measure WSPD version
    t_wspd = @elapsed d_wspd = approx_diameter_subspace(ps, ε, dir)

    # Measure FST version
    t_fst = @elapsed res_fst = MVBB.approx_diam_subspace(ps, ε, dir)
    d_fst = res_fst.distance

    return d_wspd, d_fst, t_wspd, t_fst
end

function (@main)(args)
    num_tests = 1000
    n = 1000
    ε = 0.1
    rng = Random.Xoshiro(123) # Reproducible tests
    
    @printf("Starting comparison: %d tests, n=%d, ε=%.2f\n", num_tests, n, ε)
    @printf("Points in [0,1]x[0,2]x[0,4], random directions.\n\n")

    # 1. Warm-up (JIT compilation)
    @printf("Warming up JIT...\n")
    _ = run_comparison_test(100, ε, rng)
    @printf("Warm-up complete.\n\n")

    # 2. Main Loop
    total_t_wspd = 0.0
    total_t_fst = 0.0
    max_diff = 0.0
    
    @printf("Running %d tests...\n", num_tests)
    for i in 1:num_tests
        d1, d2, t1, t2 = run_comparison_test(n, ε, rng)
        
        total_t_wspd += t1
        total_t_fst += t2
        
        diff = abs(d1 - d2)
        if diff > max_diff
            max_diff = diff
        end

        if i % 100 == 0
            @printf("  Progress: %d/%d...\n", i, num_tests)
        end
    end

    @printf("\nRESULTS:\n")
    @printf("--------------------------------------------------\n")
    @printf("Maximum diameter difference: %.10e\n", max_diff)
    @printf("Total time (WSPD-based):     %.4f seconds\n", total_t_wspd)
    @printf("Total time (FST-based):      %.4f seconds\n", total_t_fst)
    @printf("Average time per test (WSPD): %.6f seconds\n", total_t_wspd / num_tests)
    @printf("Average time per test (FST):  %.6f seconds\n", total_t_fst / num_tests)
    @printf("Overall Speedup:             %.2fx\n", total_t_fst / total_t_wspd)
    @printf("--------------------------------------------------\n")

    if max_diff < 1e-12
        @printf("SUCCESS: Results are consistent.\n")
    else
        @printf("WARNING: Results differ by %e\n", max_diff)
    end
    
    return 0
end
