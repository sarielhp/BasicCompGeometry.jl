#! /bin/env julial

"""
    Example: Nearest Neighbor Search Comparison (Accuracy and Performance)

This example compares:
1. Exact Nearest Neighbor (Linear Scan)
2. Approximate Nearest Neighbor (approx_nn)
3. "Silly" Nearest Neighbor (silly_nn - simple descent)
4. Hybrid Nearest Neighbor (hybrid_nn - silly_nn followed by approx_nn)

It reports both the approximation quality (ratios) and the total running time.
"""

using BasicCompGeometry
using BasicCompGeometry.BBT
using Printf
using Random

function (@main)(args)
    if length(args) < 3
        println("Usage: julia --project=. examples/ann.jl <n> <d> <q> [c]")
        println("  n: number of points to store in BBT")
        println("  d: dimensions")
        println("  q: number of query points")
        println("  c: approximation parameter for approx_nn (default: 1.0)")
        return 1
    end

    n = parse(Int, args[1])
    d = parse(Int, args[2])
    num_queries = parse(Int, args[3])
    c = length(args) >= 4 ? parse(Float64, args[4]) : 1.0

    @printf("Generating %d points in %d dimensions...\n", n, d)
    data_points = [rand_point(d) for _ in 1:n]
    ps = PntSeq(data_points)

    @printf("Building BBT tree...\n")
    t_build = @timed tree = Tree_init(ps)
    Tree_fully_expand(tree)
    @printf("Tree built in %.4f seconds. Depth: %d\n", t_build.time, depth(tree.root))

    # Generate query points in advance
    queries = [rand_point(d) for _ in 1:num_queries]

    # Warm-up to ensure JIT compilation doesn't skew timing
    @printf("Warming up JIT...\n")
    wp = rand_point(d)
    exact_naive_scan(ps, wp)
    approx_nn(tree, wp, c)
    silly_nn(tree, wp)
    hybrid_nn(tree, wp, c)

    @printf("\nRunning %d queries (c = %.2f)...\n", num_queries, c)
    
    # 1. Linear Scan
    @printf("Running Linear Scan (Exact)...\n")
    exact_results = Vector{Float64}(undef, num_queries)
    t_exact = @elapsed for i in 1:num_queries
        d_res, _, _ = exact_naive_scan(ps, queries[i])
        exact_results[i] = d_res
    end

    # 2. approx_nn
    @printf("Running approx_nn (ANN)...\n")
    ann_results = Vector{Float64}(undef, num_queries)
    t_ann = @elapsed for i in 1:num_queries
        d_res, _, _ = approx_nn(tree, queries[i], c)
        ann_results[i] = d_res
    end

    # 3. silly_nn
    @printf("Running silly_nn (Silly)...\n")
    silly_results = Vector{Float64}(undef, num_queries)
    t_silly = @elapsed for i in 1:num_queries
        d_res, _, _ = silly_nn(tree, queries[i])
        silly_results[i] = d_res
    end

    # 4. hybrid_nn
    @printf("Running hybrid_nn (Hybrid)...\n")
    hybrid_results = Vector{Float64}(undef, num_queries)
    t_hybrid = @elapsed for i in 1:num_queries
        d_res, _, _ = hybrid_nn(tree, queries[i], c)
        hybrid_results[i] = d_res
    end

    # Calculate Averages
    sum_ann_ratio = 0.0
    sum_silly_ratio = 0.0
    sum_hybrid_ratio = 0.0
    
    for i in 1:num_queries
        ann_ratio = exact_results[i] > 0 ? ann_results[i] / exact_results[i] : 1.0
        silly_ratio = exact_results[i] > 0 ? silly_results[i] / exact_results[i] : 1.0
        hybrid_ratio = exact_results[i] > 0 ? hybrid_results[i] / exact_results[i] : 1.0
        
        sum_ann_ratio += ann_ratio
        sum_silly_ratio += silly_ratio
        sum_hybrid_ratio += hybrid_ratio
    end

    @printf("\nQuality Comparison (Averages over %d queries):\n", num_queries)
    @printf("%-20s : %-10.4f\n", "ANN/Exact Ratio", sum_ann_ratio / num_queries)
    @printf("%-20s : %-10.4f\n", "Silly/Exact Ratio", sum_silly_ratio / num_queries)
    @printf("%-20s : %-10.4f\n", "Hybrid/Exact Ratio", sum_hybrid_ratio / num_queries)

    @printf("\nPerformance Comparison (Total Time for %d queries):\n", num_queries)
    @printf("%-15s : %10.6f seconds\n", "Linear Scan", t_exact)
    @printf("%-15s : %10.6f seconds (Speedup: %.2fx)\n", "approx_nn", t_ann, t_exact / t_ann)
    @printf("%-15s : %10.6f seconds (Speedup: %.2fx)\n", "silly_nn", t_silly, t_exact / t_silly)
    @printf("%-15s : %10.6f seconds (Speedup: %.2fx)\n", "hybrid_nn", t_hybrid, t_exact / t_hybrid)

    return 0
end
