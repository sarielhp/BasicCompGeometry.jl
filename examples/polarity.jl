#! /usr/bin/env julia

using Pkg
Pkg.activate(@__DIR__)

#
# examples/polarity.jl
#
# Demonstrate point-plane polarity with respect to the unit sphere.
# The polar of a point is a plane, and the polar of a plane is a point.
# This duality is involutive: polar(polar(x)) = x.
# In 2D, a plane is a line, so polar(Line) works the same as polar(Plane).

using BasicCompGeometry
using LinearAlgebra
using Printf

function verify_round_trip(p)
    pl = polar(p)
    q = polar(pl)
    err = dist(p, q)
    @printf("Point:      %s\n", p)
    if pl isa Line
        @printf("  Line p:   %s\n", pl.p)
        @printf("  Line u:   %s\n", pl.u)
    else
        @printf("  Plane p:  %s\n", pl.p)
        @printf("  Plane n:  %s\n", pl.n)
    end
    @printf("  Polar:    %s\n", q)
    @printf("  Error:    %.2e\n\n", err)
    return err
end

function same_plane(a::Plane{D,T}, b::Plane{D,T}) where {D,T}
    na, nb = a.n ./ norm(a.n), b.n ./ norm(b.n)
    da, db = dot(na, a.p), dot(nb, b.p)
    return isapprox(na, nb; atol=1e-12) && isapprox(da, db; atol=1e-12)
end

function main()
    println("="^60)
    println("Polarity Demo: Point-Plane Duality")
    println("="^60)

    println("\n--- 2D Polarity ---\n")
    p2 = point(2.0, 3.0)
    e2 = verify_round_trip(p2)

    p2b = point(-1.5, 4.0)
    e2b = verify_round_trip(p2b)

    println("--- 2D: polar(Line) works as polar(Plane) ---\n")
    l = Line(point(2.0, 1.0), point(1.0, 0.0))
    pl = Plane(l.p, l.u)
    pt_from_line = polar(l)
    pt_from_plane = polar(pl)
    @printf("Line:       p=%s, u=%s\n", l.p, l.u)
    @printf("  polar(Line):   %s\n", pt_from_line)
    @printf("  polar(Plane):  %s\n", pt_from_plane)
    @printf("  Match:         %s\n\n", pt_from_line ≈ pt_from_plane)

    println("--- 2D: polar of polar(Line) recovers the line ---\n")
    l = Line(point(2.0, 1.0), point(1.0, 0.0))
    ql = polar(l)
    pl_from_ql = polar(ql)
    @printf("Original Line:  p=%s, u=%s\n", l.p, l.u)
    @printf("Polar point:    %s\n", ql)
    @printf("Recovered line: p=%s, n=%s\n", pl_from_ql.p, pl_from_ql.n)
    @printf("Line eqn: dot([%s], x) = %.1f\n", l.u, dot(l.u, l.p))
    @printf("Recovered: dot([%s], x) = %.1f\n", pl_from_ql.n, dot(pl_from_ql.n, pl_from_ql.p))
    @printf("Match: %s\n\n",
        isapprox(l.u ./ norm(l.u), pl_from_ql.n ./ norm(pl_from_ql.n)) &&
        isapprox(dot(l.u ./ norm(l.u), l.p), dot(pl_from_ql.n ./ norm(pl_from_ql.n), pl_from_ql.p)))

    println("--- 3D Polarity ---\n")
    p3 = point(1.0, 2.0, 3.0)
    e3 = verify_round_trip(p3)

    p3b = point(-3.0, 0.5, 2.0)
    e3b = verify_round_trip(p3b)

    println("--- Higher Dimensions ---\n")
    p5 = Point{5,Float64}((1.0, 2.0, 3.0, 4.0, 5.0))
    e5 = verify_round_trip(p5)

    println("--- Plane->Point->Plane Round Trip ---\n")
    pl = Plane(point(1.0, 0.0, 0.0), point(1.0, 1.0, 0.0))
    q = polar(pl)
    pl2 = polar(q)
    @printf("Original plane: p=%s, n=%s\n", pl.p, pl.n)
    @printf("Polar point:    %s\n", q)
    @printf("Polar plane p:  %s\n", pl2.p)
    @printf("Polar plane n:  %s\n", pl2.n)
    @printf("Plane p match:  %s (equivalent representation)\n", same_plane(pl, pl2))

    max_err = max(e2, e2b, e3, e3b, e5)
    println("\n" * "="^60)
    @printf("Maximum round-trip error: %.2e\n", max_err)
    println("All round-trip tests passed!")

    return max_err
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end