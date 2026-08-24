using Test
using BasicCompGeometry
using LinearAlgebra

@testset "Parabola" begin
    @testset "Construction" begin
        p = point(0.3, 0.4)
        line = Line(point(0.0, 0.0), point(0.0, 1.0))
        par = Parabola(p, line)
        @test par isa Parabola{Float64}
        @test par isa AbsCurve2D{Float64}
        @test par isa Curve2D{Float64}
        @test vertex(par) ≈ point(0.15, 0.4)
        @test p_coeff(par) ≈ 0.15
        @test axis_direction(par) ≈ point(1.0, 0.0)
    end

    @testset "Point on parabola" begin
        p = point(0.3, 0.4)
        line = Line(point(0.0, 0.0), point(0.0, 1.0))
        par = Parabola(p, line)
        for t in [-2.0, -1.0, 0.0, 1.0, 2.0]
            x = at(par, t)
            d_p = dist(x, p)
            d_line = x.x
            @test abs(d_p - d_line) < 1e-10
        end
    end

    @testset "Intersect line with parabola" begin
        p = point(0.3, 0.4)
        line = Line(point(0.0, 0.0), point(0.0, 1.0))
        par = Parabola(p, line)
        bis = Line(point(0.0, 0.0), point(1.0, 1.0))
        pts = intersect_line_curve(bis, par)
        @test length(pts) == 2
        for (pt, _) in pts
            @test abs(dist(pt, p) - pt.x) < 1e-10
        end

        # Line parallel to axis of parabola (x = 0.5)
        vert_line = Line(point(0.5, 0.0), point(1.0, 0.0))
        pts_vert = intersect_line_curve(vert_line, par)
        @test length(pts_vert) >= 1
    end

    @testset "geom_length and eltype" begin
        p = point(0.3, 0.4)
        line = Line(point(0.0, 0.0), point(0.0, 1.0))
        par = Parabola(p, line)
        @test isinf(geom_length(par))
        @test eltype(par) == Float64
    end

    @testset "Bisector function" begin
        line = Line(point(0.0, 1.0), point(1.0, 0.0))
        pt = point(0.0, 0.0)
        p = bisector(line, pt)
        @test p isa Parabola{Float64}
        @test p isa Curve2D{Float64}
    end

    @testset "ParabolaArc" begin
        para = Parabola(point(0.0, 1.0), Line(point(0.0, -1.0), point(1.0, 0.0)))
        arc = ParabolaArc(para, point(-2.0, 1.0), point(2.0, 1.0))
        @test arc isa ParabolaArc{Float64}
        @test arc isa AbsCurve2D{Float64}
        @test arc isa Curve2D{Float64}
        @test eltype(arc) == Float64
        @test at(arc, 0.0) ≈ point(-2.0, 1.0)
        @test at(arc, 1.0) ≈ point(2.0, 1.0)
        @test at(arc, 0.5) ≈ point(0.0, 0.0)
        @test geom_length(arc) > 4.0

        h_line = Line(point(0.0, 0.5), point(1.0, 0.0))
        pts = intersect_line_curve(h_line, arc)
        @test length(pts) == 2
        for (pt, s) in pts
            @test pt.y ≈ 0.5
            @test -2.0 <= pt.x <= 2.0
            @test 0.0 <= s <= 1.0
        end
    end
end
