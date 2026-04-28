using Test
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
end
