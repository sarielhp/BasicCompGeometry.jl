using Test
using BasicCompGeometry
using BasicCompGeometry.MetricSpace
using LinearAlgebra

@testset "MetricSpace" begin
    @testset "PointsSpace" begin
        pnts = [point(0.0, 0.0), point(3.0, 4.0)]
        ps = PointsSpace(pnts)
        @test size(ps) == 2
        @test dist(ps, 1, 2) ≈ 5.0
        @test dist(ps, 1, 1) == 0.0
        @test dist_real(ps, 1, point(1.0, 0.0)) ≈ 1.0
    end

    @testset "MPointsSpace" begin
        m = [0.0 3.0; 0.0 4.0]
        mps = MPointsSpace(m)
        @test size(mps) == 2
        @test dist(mps, 1, 2) ≈ 5.0
        @test dist_real(mps, 1, [1.0, 0.0]) ≈ 1.0
    end

    @testset "PermutMetric" begin
        pnts = [point(0.0, 0.0), point(3.0, 4.0)]
        ps = PointsSpace(pnts)
        pm = PermutMetric(ps)
        @test size(pm) == 2
        @test dist(pm, 1, 2) ≈ 5.0

        MetricSpace.swap!(pm, 1, 2)
        @test original(pm, 1) == 2
        @test original(pm, 2) == 1
        @test dist(pm, 1, 2) ≈ 5.0
    end

    @testset "SpherePMetric" begin
        pnts = [point(0.0, 0.0), point(1.0, 0.0), point(0.0, 1.0)]
        ps = PointsSpace(pnts)
        # Base point is pnts[1] (0,0)
        # pnts[2] is (1,0), dist to base is 1.0
        # pnts[3] is (0,1), dist to base is 1.0
        # dist(pnts[2], pnts[3]) is sqrt(2)
        # angle should be π/2

        sm = SpherePMetric(ps, 1)
        push!(sm, 2, 1.0)
        push!(sm, 3, 1.0)

        @test size(sm) == 2
        @test dist(sm, 1, 2) ≈ π/2
    end

    @testset "Greedy Permutation" begin
        pnts = [point(0.0, 0.0), point(1.0, 0.0), point(0.5, 0.5), point(10.0, 10.0)]
        ps = PointsSpace(pnts)

        I, D, N = greedy_permutation_naive(ps, 4)
        @test I[1] == 1
        @test I[2] == 4 # furthest from 1 is (10,10)
        @test D[1] ≈ dist(pnts[1], pnts[4])
    end

    @testset "Polygon as AbsFMS" begin
        pnts = [point(0.0, 0.0), point(3.0, 4.0), point(0.0, 4.0)]
        poly = Polygon(pnts)
        @test poly isa AbsFMS
        @test size(poly) == 3
        @test dist(poly, 1, 2) ≈ 5.0

        # Test greedy permutation directly on Polygon
        I, D, N = greedy_permutation_naive(poly, 2)
        @test I[1] == 1
        @test I[2] == 2 # (3,4) is further from (0,0) than (0,4)
    end
end

@testset "Metric adapters and angular edge cases" begin
    points = [point(0.0, 0.0), point(3.0, 4.0), point(6.0, 0.0)]
    space = PointsSpace(points)
    @test metric(space, 1, 2) == 5.0
    @test original(space, 3) == 3
    @test dist(MPointsSpace(Matrix(PntSeq(points)), 2), 1, 1) == 0.0
    @test_throws DimensionMismatch MetricSpace.euclidean_distance([1.0], [1.0, 2.0])
    @test_throws AssertionError MPointsSpace(zeros(2, 2), 3)

    indexes = [3, 1, 2]
    permuted = PermutMetric(space, indexes)
    indexes[1] = 1
    @test original(permuted, 1) == 3 # The constructor owns its mapping.
    @test dist(permuted, 1, 1) == 0.0
    @test dist_real(permuted, 1, point(6.0, 2.0)) == 2.0
    MetricSpace.swap!(permuted, 2, 2)
    @test permuted.I == [3, 1, 2]

    angular = SpherePMetric(space, 1)
    push!(angular, 2, 5.0)
    push!(angular, 3, 6.0)
    theta = dist(angular, 1, 2)
    @test dist(angular, 1, 1) == 0.0
    MetricSpace.swap!(angular, 1, 2)
    @test angular.I == [3, 2]
    @test angular.D == [6.0, 5.0]
    @test dist(angular, 1, 2) ≈ theta
    MetricSpace.swap!(angular, 1, 1)
    @test angular.I == [3, 2]

    angle_of = MetricSpace.calculate_angle_facing_a
    @test angle_of(0.0, 1.0, 1.0) == 0.0
    @test angle_of(1.0, 0.0, 1.0) == π / 2
    @test angle_of(2.0, 1.0, 1.0) ≈ π
    @test_throws ArgumentError angle_of(-1.0, 1.0, 1.0)
end
