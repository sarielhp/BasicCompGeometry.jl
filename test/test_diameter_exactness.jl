using Test
using BasicCompGeometry
using BasicCompGeometry.MVBB
using LinearAlgebra
using Random
using Printf

@testset "Diameter Exactness Comparison (n=10000, ε=0)" begin
    # Generate 10000 points uniformly on the 3D unit sphere
    n = 10000
    rng = Random.Xoshiro(42) # Fixed seed for reproducibility
    pts = Point{3, Float64}[]
    for _ in 1:n
        p = Point{3, Float64}(randn(rng, 3))
        push!(pts, normalize(p))
    end
    ps = PntSeq(pts)

    # 1. Exact Diameter (Naive O(n^2))
    d_exact = exact_diameter(ps)

    # 2. WSPD-based Diameter with ε=0
    # Note: approx_diameter uses ε/2 as the separation parameter.
    d_wspd = approx_diameter(ps, 0.0)

    # 3. Fair Split Tree Traversal with ε=0
    res_fst = MVBB.approx_diam(ps, 0.0)
    d_fst = res_fst.distance

    @printf("Exact: %.10f\n", d_exact)
    @printf("WSPD:  %.10f\n", d_wspd)
    @printf("FST:   %.10f\n", d_fst)

    # They should all be exactly equal
    @test d_wspd == d_exact
    @test d_fst == d_exact
end
