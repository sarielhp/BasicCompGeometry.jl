using Test
using BasicCompGeometry
using LinearAlgebra

@testset "Hyperbola" begin
    @testset "Construction" begin
        h = Hyperbola(point(0.0, 0.0), 5.0, 3.0, 0.0)
        @test h isa Hyperbola{Float64}
        @test h isa AbsCurve2D{Float64}
        @test center(h) ≈ point(0.0, 0.0)
        @test a_coeff(h) ≈ 5.0
        @test b_coeff(h) ≈ 3.0
        @test c_coeff(h) ≈ sqrt(34.0)
        @test rotation_angle(h) ≈ 0.0
    end

    @testset "From foci" begin
        f1 = point(-5.0, 0.0)
        f2 = point(5.0, 0.0)
        h = Hyperbola(f1, f2, 6.0)
        @test center(h) ≈ point(0.0, 0.0)
        @test a_coeff(h) ≈ 3.0
        @test b_coeff(h) ≈ 4.0
        @test c_coeff(h) ≈ 5.0
    end

    @testset "From line and point (bisector)" begin
        line = Line(point(0.0, 1.0), point(1.0, 0.0))
        pt = point(0.0, 0.0)
        p = Parabola(pt, line)
        @test p isa Parabola{Float64}
        @test p isa AbsCurve2D{Float64}
        @test p.focus == pt
        @test p.directrix == line
    end

    @testset "Point on hyperbola" begin
        h = Hyperbola(point(0.0, 0.0), 5.0, 3.0, 0.0)
        p0 = at(h, 0.0; branch=1)
        @test p0 ≈ point(5.0, 0.0)
        p0b = at(h, 0.0; branch=2)
        @test p0b ≈ point(-5.0, 0.0)
        p1 = at(h, 1.0; branch=1)
        @test dist(p1, point(5.0 * cosh(1.0), 3.0 * sinh(1.0))) < 1e-10
    end

    @testset "Rotated hyperbola" begin
        angle = pi / 4
        h = Hyperbola(point(1.0, 2.0), 4.0, 2.0, angle)
        @test center(h) ≈ point(1.0, 2.0)
        @test a_coeff(h) ≈ 4.0
        @test b_coeff(h) ≈ 2.0
        @test rotation_angle(h) ≈ angle
        p0 = at(h, 0.0; branch=1)
        expected = point(1.0, 2.0) + 4.0 * point(cos(angle), sin(angle))
        @test p0 ≈ expected
    end

    @testset "Intersect line with hyperbola" begin
        h = Hyperbola(point(0.0, 0.0), 5.0, 3.0, 0.0)
        line = Line(point(0.0, 0.0), point(1.0, 0.0))
        pts = intersect_line_curve(line, h)
        @test length(pts) == 2
        xs = sort([p[1] for (p, _) in pts])
        @test xs[1] ≈ -5.0
        @test xs[2] ≈ 5.0
    end

    @testset "Intersect vertical line" begin
        h = Hyperbola(point(0.0, 0.0), 5.0, 3.0, 0.0)
        line = Line(point(0.0, 0.0), point(0.0, 1.0))
        pts = intersect_line_curve(line, h)
        @test length(pts) == 0
    end

    @testset "Curve2D union type" begin
        h = Hyperbola(point(0.0, 0.0), 5.0, 3.0, 0.0)
        @test h isa Curve2D{Float64}
        @test eltype(h) == Float64
    end

    @testset "geom_length" begin
        h = Hyperbola(point(0.0, 0.0), 5.0, 3.0, 0.0)
        @test isinf(geom_length(h))
    end

    @testset "Bisector function" begin
        line = Line(point(0.0, 1.0), point(1.0, 0.0))
        pt = point(0.0, 0.0)
        p = bisector(line, pt)
        @test p isa Parabola{Float64}
        @test p isa Curve2D{Float64}
    end

    @testset "Distance to hyperbola" begin
        h = Hyperbola(point(0.0, 0.0), 5.0, 3.0, 0.0)
        d = distance(point(5.0, 0.0), h)
        @test d < 0.1
        d2 = distance(point(0.0, 0.0), h)
        @test d2 > 4.0
    end
end

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
            d_line = x[1]
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
            @test abs(dist(pt, p) - pt[1]) < 1e-10
        end
    end

    @testset "geom_length and eltype" begin
        p = point(0.3, 0.4)
        line = Line(point(0.0, 0.0), point(0.0, 1.0))
        par = Parabola(p, line)
        @test isinf(geom_length(par))
        @test eltype(par) == Float64
    end
end