using Test
using Random
using BasicCompGeometry
using BasicCompGeometry.LongestConvexSubset

@testset "LongestConvexSubset" begin
    @testset "Simple square" begin
        # 4 points forming a square and 1 interior point
        pnts = [
            point(0.0, 0.0),
            point(1.0, 0.0),
            point(1.0, 1.0),
            point(0.0, 1.0),
            point(0.5, 0.5),
        ]
        poly = Polygon(pnts)

        subset = compute_largest_convex_subset(poly)
        @test length(subset) == 4
        # The interior point should be excluded
        @test !(point(0.5, 0.5) in subset)
    end

    @testset "Chains" begin
        # Points: (0,0), (1, 1), (2,0), (1, -1)
        pnts = [point(0.0, 0.0), point(1.0, 1.0), point(2.0, 0.0), point(1.0, -1.0)]
        poly = Polygon(pnts)

        # Upper chain should be (0,0), (1,1), (2,0)
        u_chain = longest_convex_chain(poly)
        @test length(u_chain) == 3
        # Check that it contains the expected points (order might depend on sorting)
        @test point(0.0, 0.0) in u_chain
        @test point(1.0, 1.0) in u_chain
        @test point(2.0, 0.0) in u_chain

        # Lower chain should be (0,0), (1,-1), (2,0)
        l_chain = longest_concave_chain(poly)
        @test length(l_chain) == 3
        @test point(0.0, 0.0) in l_chain
        @test point(1.0, -1.0) in l_chain
        @test point(2.0, 0.0) in l_chain
    end

    @testset "Integer coordinates" begin
        # 4 points forming a square and 1 interior point using Integers
        pnts = [point(0, 0), point(10, 0), point(10, 10), point(0, 10), point(5, 5)]
        poly = Polygon(pnts)

        subset = compute_largest_convex_subset(poly)
        @test length(subset) == 4
        @test eltype(subset[1]) <: Integer
        @test !(point(5, 5) in subset)
    end

    @testset "Concrete context types and context reuse" begin
        rng = Xoshiro(7)
        ps = PntSeq([point(rand(rng), rand(rng)) for _ = 1:60])

        lct = LongestConvexSubset.lcs_init(ps)
        ro = LongestConvexSubset.RotateOrder(ps, ps[1])
        @test isconcretetype(fieldtype(typeof(lct), :PS))
        @test isconcretetype(fieldtype(typeof(ro), :PS))

        # Hoisted context initialization should produce identical results to per-target init
        ps25 = PntSeq([point(rand(rng), rand(rng)) for _ = 1:25])
        subset25 = compute_largest_convex_subset(ps25)
        @test length(subset25) >= 3
        work_ps = deepcopy(ps25)
        LongestConvexSubset.sort_points!(work_ps)
        ps_ref = PntSeq{2,Float64}([point(p[1], -p[2]) for p in work_ps])
        best_sz = -1
        best_sol = Int[]
        for target = 1:length(ps25)
            lu = LongestConvexSubset.lcs_init(work_ps)
            LongestConvexSubset.compute_dp_table!(lu, target)
            ll = LongestConvexSubset.lcs_init(ps_ref)
            LongestConvexSubset.compute_dp_table!(ll, target)
            bs = LongestConvexSubset._best_start(lu, ll)
            bs == -1 && continue
            comb = LongestConvexSubset._combined_loop(lu, ll, bs, target)
            if length(comb) > best_sz
                best_sz = length(comb)
                best_sol = comb
            end
        end
        @test subset25.pnts == [work_ps[i] for i in best_sol]

        # Reusing lct_upper/lct_lower avoids O(n^3) context allocations (~192 MB -> ~175 KB)
        compute_largest_convex_subset(ps)
        alloc = @allocated compute_largest_convex_subset(ps)
        @test alloc < 2_000_000
    end
end
