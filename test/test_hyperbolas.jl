using Test
using BasicCompGeometry
import BasicCompGeometry: a_coeff, at, b_coeff, c_coeff, center, rotation_angle
using BasicCompGeometry: distance
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

    @testset "Point on hyperbola" begin
        h = Hyperbola(point(0.0, 0.0), 5.0, 3.0, 0.0)
        p0 = at(h, 0.0; branch = 1)
        @test p0 ≈ point(5.0, 0.0)
        p0b = at(h, 0.0; branch = 2)
        @test p0b ≈ point(-5.0, 0.0)
        p1 = at(h, 1.0; branch = 1)
        @test dist(p1, point(5.0 * cosh(1.0), 3.0 * sinh(1.0))) < 1e-10
    end

    @testset "Rotated hyperbola" begin
        angle = pi / 4
        h = Hyperbola(point(1.0, 2.0), 4.0, 2.0, angle)
        @test center(h) ≈ point(1.0, 2.0)
        @test a_coeff(h) ≈ 4.0
        @test b_coeff(h) ≈ 2.0
        @test rotation_angle(h) ≈ angle
        p0 = at(h, 0.0; branch = 1)
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

    @testset "Distance to hyperbola" begin
        h = Hyperbola(point(0.0, 0.0), 5.0, 3.0, 0.0)
        d = distance(point(5.0, 0.0), h)
        @test d < 0.1
        d2 = distance(point(0.0, 0.0), h)
        @test d2 > 4.0
    end
end

@testset "Hyperbola convenience methods and intersections" begin
    h = Hyperbola([0.0, 0.0], 2.0, 1.0, 0.0)
    @test a_coeff(h) ≈ 2.0
    @test occursin("Hyperbola(f1=", sprint(show, h))
    line = Line(point(0.0, -1.0), point(1.0, 0.0))
    parabola = Hyperbola(line, point(0.0, 1.0))
    @test parabola isa Parabola
    @test at(parabola, 0.0) ≈ point(0.0, 0.0)
    asymptote = Line(point(0.0, 0.0), point(2.0, 1.0))
    @test isempty(intersect_line_curve(asymptote, h))

    t = range(-3.0, 3.0; length = 200)[130]
    other = Hyperbola(point(0.0, 2sinh(t)), 2.0, 1.0, 0.0)
    hits = intersect_hyperbolas(h, other)
    @test any(p -> p ≈ point(2cosh(t), sinh(t)), hits)
    @test all(p -> abs(abs(dist(p, other.f1) - dist(p, other.f2)) - other.d) < 0.01, hits)
    @test pairwise_intersections(h, other) == hits
    @test BasicCompGeometry.vertices([h, other]; unique = true) == hits
end
