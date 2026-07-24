using Test
using BasicCompGeometry

@testset "Spheres & Circles" begin
    # Constructors and type aliases
    s3d = Sphere(point(1.0, 2.0, 3.0), 5.0)
    @test s3d.center == point(1.0, 2.0, 3.0)
    @test s3d.radius == 5.0
    @test typeof(s3d) == Sphere{3,Float64}

    c2d = Circle(point(1.0, 2.0), 3.0)
    @test c2d.center == point(1.0, 2.0)
    @test c2d.radius == 3.0
    @test typeof(c2d) == Circle{Float64}
    @test typeof(c2d) == Sphere{2,Float64}

    c_tuple = Circle((1, 2), 3)
    @test c_tuple.center == point(1, 2)
    @test c_tuple.radius == 3

    # is_inside & dist
    @test is_inside(point(1.0, 2.0), c2d) == true
    @test is_inside(point(4.0, 2.0), c2d) == true
    @test is_inside(point(5.0, 2.0), c2d) == false
    @test dist(point(5.0, 2.0), c2d) ≈ 1.0
    @test dist(point(1.0, 2.0), c2d) == 0.0

    # Point inversion
    ref = Circle(point(0.0, 0.0), 1.0)
    p = point(2.0, 0.0)
    @test invert(p, ref) ≈ point(0.5, 0.0)
    @test invert(invert(p, ref), ref) ≈ p

    # Line inversion: Line NOT passing through origin -> Circle
    l_sec = Line(point(2.0, 0.0), point(0.0, 1.0)) # x = 2
    inv_l = invert(l_sec, ref)
    @test inv_l isa Circle{Float64}
    @test inv_l.center ≈ point(0.25, 0.0)
    @test inv_l.radius ≈ 0.25

    # Invert back: Circle passing through origin -> Line
    l_back = invert(inv_l, ref)
    @test l_back isa Line{2,Float64}
    # Check that l_back is equivalent to x = 2
    @test BasicCompGeometry.distance(point(2.0, 0.0), l_back) ≈ 0.0 atol = 1e-10
    @test BasicCompGeometry.distance(point(2.0, 5.0), l_back) ≈ 0.0 atol = 1e-10

    # Line inversion: Line passing through origin -> Line
    l_orig = Line(point(0.0, 0.0), point(1.0, 1.0))
    inv_l_orig = invert(l_orig, ref)
    @test inv_l_orig isa Line{2,Float64}
    @test BasicCompGeometry.distance(point(1.0, 1.0), inv_l_orig) ≈ 0.0 atol = 1e-10

    # Circle inversion: Circle NOT passing through origin -> Circle
    c_off = Circle(point(3.0, 0.0), 1.0)
    inv_c_off = invert(c_off, ref)
    @test inv_c_off isa Circle{Float64}
    @test inv_c_off.center ≈ point(0.375, 0.0)
    @test inv_c_off.radius ≈ 0.125

    # Invert back -> Circle
    c_back = invert(inv_c_off, ref)
    @test c_back isa Circle{Float64}
    @test c_back.center ≈ c_off.center atol = 1e-10
    @test c_back.radius ≈ c_off.radius atol = 1e-10

    # Line <-> Plane conversions in 2D
    l_test = Line(point(1.0, 0.0), point(0.0, 1.0)) # x = 1, direction (0,1)
    pl_test = Plane(l_test)
    @test pl_test.n ≈ point(1.0, 0.0) # normal perpendicular to (0,1)
    l_converted = Line(pl_test)
    @test l_converted.p == l_test.p
    @test l_converted.u ≈ l_test.u

    # polar(Line, p) & polar_line(p)
    pt = point(2.0, 0.0)
    pl_line = polar(Line, pt)
    @test pl_line isa Line{2, Float64}
    @test polar_line(pt) == pl_line
    @test BasicCompGeometry.distance(point(0.5, 0.0), pl_line) ≈ 0.0 atol = 1e-10
    @test BasicCompGeometry.distance(point(0.5, 3.0), pl_line) ≈ 0.0 atol = 1e-10

    # Plane distance
    @test BasicCompGeometry.distance(point(3.0, 0.0), pl_test) ≈ 2.0

    # CircleArc & Segment Inversion
    arc1 = CircleArc(point(0.25, 0.0), 0.25, 0.0, pi / 2)
    @test arc1.center == point(0.25, 0.0)
    @test arc1.radius == 0.25
    @test typeof(arc1) == CircleArc2F

    seg = Segment(point(2.0, 0.0), point(2.0, 2.0))
    inv_seg = invert(seg, ref)
    @test inv_seg isa CircleArc2F
    @test inv_seg.center ≈ point(0.25, 0.0)
    @test inv_seg.radius ≈ 0.25
    @test inv_seg.theta1 ≈ 0.0 atol = 1e-10
    @test inv_seg.theta2 ≈ pi / 2 atol = 1e-10
end
