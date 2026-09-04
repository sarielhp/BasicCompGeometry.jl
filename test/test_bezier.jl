using Test
using BasicCompGeometry

@testset "CubicBezier and CubicSpline Primitives" begin
    # 1. Cubic Bézier evaluation
    p0 = Point(0.0, 0.0)
    p1 = Point(0.0, 2.0)
    p2 = Point(2.0, 2.0)
    p3 = Point(2.0, 0.0)
    b = CubicBezier(p0, p1, p2, p3)

    @test b(0.0) ≈ p0
    @test b(1.0) ≈ p3
    mid = b(0.5)
    @test mid[1] ≈ 1.0
    @test mid[2] ≈ 1.5

    # 2. Derivative
    d0 = derivative(b, 0.0)
    @test d0[1] ≈ 0.0
    @test d0[2] ≈ 6.0 # 3 * (p1 - p0)

    # 3. Subdivision
    left, right = subdivide(b, 0.5)
    @test left == Base.split(b, 0.5)[1]
    @test left.p0 ≈ p0
    @test left.p3 ≈ mid
    @test right.p0 ≈ mid
    @test right.p3 ≈ p3
    @test left(0.5) ≈ b(0.25)
    @test right(0.5) ≈ b(0.75)

    # 4. Exact Bounding Box
    bb = bbox(b)
    @test bottom_left(bb)[1] ≈ 0.0 atol=1e-10
    @test top_right(bb)[1] ≈ 2.0 atol=1e-10
    @test bottom_left(bb)[2] ≈ 0.0 atol=1e-10
    @test top_right(bb)[2] ≈ 1.5 atol=1e-10 # Maximum y is at mid (t=0.5)

    # All sampled points must lie strictly inside bbox
    for t in 0.0:0.05:1.0
        pt = b(t)
        @test is_inside(pt, bb)
    end

    # 5. Flattening
    flat = flatten(b; tol=1e-3)
    @test cardin(flat) >= 4
    @test flat[1] ≈ p0
    @test flat[cardin(flat)] ≈ p3

    # 6. Catmull-Rom Interpolation
    knots = PntSeq([Point(0.0, 0.0), Point(2.0, 3.0), Point(4.0, 1.0), Point(6.0, 4.0)])
    cr_spline = interpolate_catmull_rom(knots)
    @test length(cr_spline) == 3

    # Must exactly interpolate all knot points
    @test cr_spline(0.0) ≈ knots[1]
    @test cr_spline(1.0) ≈ knots[2]
    @test cr_spline(2.0) ≈ knots[3]
    @test cr_spline(3.0) ≈ knots[4]

    # 7. Natural Cubic Spline Interpolation
    nat_spline = interpolate_natural_spline(knots)
    @test length(nat_spline) == 3

    # Must exactly interpolate all knot points
    @test nat_spline(0.0) ≈ knots[1]
    @test nat_spline(1.0) ≈ knots[2]
    @test nat_spline(2.0) ≈ knots[3]
    @test nat_spline(3.0) ≈ knots[4]

    # Spline bbox contains all knots
    s_bb = bbox(nat_spline)
    for p in knots
        @test is_inside(p, s_bb)
    end

    # Spline flattening
    s_flat = flatten(nat_spline; tol=1e-2)
    @test cardin(s_flat) > cardin(knots)
    @test s_flat[1] ≈ knots[1]
    @test s_flat[cardin(s_flat)] ≈ knots[4]
end
