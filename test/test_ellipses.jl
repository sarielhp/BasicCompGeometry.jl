using Test
using BasicCompGeometry
import BasicCompGeometry: cardin, center, end_point, r_major, r_minor, start_point
using BasicCompGeometry: area

@testset "Ellipse and EllipticArc Primitives" begin
    # 1. Basic Ellipse creation and properties
    c = Point(10.0, 20.0)
    e = Ellipse(c, 5.0, 3.0, 0.0)
    @test center(e) == c
    @test r_major(e) == 5.0
    @test r_minor(e) == 3.0
    @test angle(e) == 0.0
    @test area(e) ≈ π * 15.0

    # 2. Point containment
    @test is_inside(c, e)
    @test is_inside(Point(15.0, 20.0), e)
    @test is_inside(Point(10.0, 23.0), e)
    @test !is_inside(Point(15.1, 20.0), e)
    @test !is_inside(Point(10.0, 23.1), e)
    @test !is_inside(Point(0.0, 0.0), e)

    # 3. Parametric evaluation
    @test point_on(e, 0.0) ≈ Point(15.0, 20.0)
    @test point_on(e, π / 2) ≈ Point(10.0, 23.0)
    @test point_on(e, π) ≈ Point(5.0, 20.0)
    @test point_on(e, 3π / 2) ≈ Point(10.0, 17.0)

    # 4. Axis-aligned BBox
    bb = bbox(e)
    @test bottom_left(bb) ≈ Point(5.0, 17.0)
    @test top_right(bb) ≈ Point(15.0, 23.0)

    # 5. Rotated Ellipse BBox
    e_rot = Ellipse(Point(0.0, 0.0), 4.0, 2.0, π / 4)
    bb_rot = bbox(e_rot)
    # semi-axes rotated by 45 deg: half_w = sqrt(16 * 0.5 + 4 * 0.5) = sqrt(10)
    expected_half = sqrt(10.0)
    @test bottom_left(bb_rot)[1] ≈ -expected_half atol=1e-10
    @test top_right(bb_rot)[1] ≈ expected_half atol=1e-10
    @test bottom_left(bb_rot)[2] ≈ -expected_half atol=1e-10
    @test top_right(bb_rot)[2] ≈ expected_half atol=1e-10

    # 6. Uniform sampling
    pts = sample_uniformly(e, 12)
    @test cardin(pts) == 12
    for p in pts
        @test is_inside(p, e)
    end

    # 7. EllipticArc properties
    arc = EllipticArc(e, 0.0, π / 2)
    @test start_point(arc) ≈ Point(15.0, 20.0)
    @test end_point(arc) ≈ Point(10.0, 23.0)
    @test point_on(arc, 0.5) ≈ point_on(e, π / 4)

    # 8. EllipticArc BBox
    bb_arc = bbox(arc)
    @test bottom_left(bb_arc)[1] ≈ 10.0 atol=1e-10
    @test top_right(bb_arc)[1] ≈ 15.0 atol=1e-10
    @test bottom_left(bb_arc)[2] ≈ 20.0 atol=1e-10
    @test top_right(bb_arc)[2] ≈ 23.0 atol=1e-10

    # 9. Uniform sampling on arc
    arc_pts = sample_uniformly(arc, 5)
    @test cardin(arc_pts) == 5
    @test arc_pts[1] ≈ start_point(arc)
    @test arc_pts[5] ≈ end_point(arc)
end

@testset "Ellipse convenience constructors and wrapped arcs" begin
    e = Ellipse((1, 2.0), 4, 2, π / 4)
    @test e isa Ellipse2F
    @test BasicCompGeometry.rotation_angle(e) ≈ π / 4
    @test BBox(e) ≈
          BBox(point(1 - sqrt(10), 2 - sqrt(10)), point(1 + sqrt(10), 2 + sqrt(10)))
    @test occursin("Ellipse(", sprint(show, e))
    arc = EllipticArc(point(1.0, 2.0), 4.0, 2.0, π / 4, 3π / 2, π)
    @test center(arc) == point(1.0, 2.0)
    @test r_major(arc) == 4.0
    @test r_minor(arc) == 2.0
    @test angle(arc) ≈ π / 4
    @test BasicCompGeometry.alpha1(arc) ≈ 3π / 2
    @test BasicCompGeometry.alpha2(arc) ≈ π
    @test point_on(arc, 1.0) ≈ end_point(arc)
    @test occursin("EllipticArc(", sprint(show, arc))
    box = BBox(arc)
    @test all(p -> is_inside(p, box + 1e-10), sample_uniformly(arc, 101))
    singleton = EllipticArc(e, 0.0, 0.0)
    @test bottom_left(bbox(singleton)) == top_right(bbox(singleton))
end
