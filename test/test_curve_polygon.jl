using Test
using BasicCompGeometry
using LinearAlgebra

@testset "CurvePolygon2D and point_on tests" begin
    # Test CircleArc geom_length and at
    arc = CircleArc(point(0.0, 0.0), 1.0, 0.0, pi / 2)
    @test geom_length(arc) ≈ pi / 2
    @test at(arc, 0.0) ≈ point(1.0, 0.0)
    @test at(arc, 1.0) ≈ point(0.0, 1.0)
    @test at(arc, 0.5) ≈ point(cos(pi / 4), sin(pi / 4))

    # Test Segment geom_length and at
    seg = Segment(point(0.0, 1.0), point(0.0, 0.0))
    @test geom_length(seg) == 1.0
    @test at(seg, 0.5) ≈ point(0.0, 0.5)

    # Test CurvePolygon2D
    seg2 = Segment(point(0.0, 0.0), point(1.0, 0.0))
    cp = CurvePolygon2D([seg2, arc, seg])
    @test cp isa CurvePolygon2D{Float64}
    @test geom_length(cp) ≈ 1.0 + pi / 2 + 1.0

    # Test point_on for CurvePolygon2D
    @test point_on(cp, 0.0) ≈ point(0.0, 0.0)
    @test point_on(cp, 1.0) ≈ point(0.0, 0.0) # end of seg
    
    total_l = geom_length(cp)
    t_arc_mid = (1.0 + (pi / 4)) / total_l
    @test point_on(cp, t_arc_mid) ≈ point(cos(pi / 4), sin(pi / 4))

    # Test point_on for regular polygon (PntSeq / AbsPntSeq)
    pnts = [point(0.0, 0.0), point(2.0, 0.0), point(2.0, 2.0), point(0.0, 2.0)]
    poly = PntSeq(pnts)
    # Closed parameterization: total length 8.0
    @test point_on(poly, 0.0, closed=true) ≈ point(0.0, 0.0)
    @test point_on(poly, 0.25, closed=true) ≈ point(2.0, 0.0)
    @test point_on(poly, 0.5, closed=true) ≈ point(2.0, 2.0)
    @test point_on(poly, 0.75, closed=true) ≈ point(0.0, 2.0)
    # Test direction and tangent_line
    v = direction(pi / 4)
    @test v ≈ point(cos(pi / 4), sin(pi / 4))

    line_t, q_t = tangent_line(poly, 0.0)
    @test q_t ≈ point(2.0, 0.0) # max x vertex

    u_t = before_tangent_to_polygon(cp, poly, 0.0)
    v_t = after_tangent_to_polygon(cp, poly, 0.0)
    @test 0.0 <= u_t <= 1.0
    @test 0.0 <= v_t <= 1.0
end
