using Test
using Random
using BasicCompGeometry

# Small-polygon fast path of `triangulate_ear_cut` (3 <= n <= 64, src/SmallEarCut.jl) and
# the bounding-box rejection in `_point_in_ccw_triangle`.

# Random simple CCW star-shaped polygon: sorted angles with every gap below pi.
function _small_star(rng, n)
    while true
        th = sort(rand(rng, n)) .* 2pi
        maximum(vcat(diff(th), th[1] + 2pi - th[end])) < pi || continue
        r = 0.5 .+ rand(rng, n)
        return [point(r[i] * cos(th[i]), r[i] * sin(th[i])) for i in 1:n]
    end
end

# Valid triangulation: n - 2 non-clockwise triangles covering the polygon's area and n - 3
# ordered diagonals.
function _small_valid(P, res; kernel = FloatingKernel())
    n = length(P)
    length(res.triangles) == n - 2 || return false
    length(res.diagonals) == n - 3 || return false
    all(d -> d[1] < d[2], res.diagonals) || return false
    tri(t) = [P[t[1]], P[t[2]], P[t[3]]]
    all(t -> orient_sign(P[t[1]], P[t[2]], P[t[3]]; kernel = kernel) >= 0, res.triangles) ||
        return false
    total = sum(twice_area(tri(t); kernel = kernel) for t in res.triangles)
    return kernel isa ExactKernel ? total == twice_area(P; kernel = kernel) :
           isapprox(total, twice_area(P); rtol = 1e-9)
end

@testset "Small-polygon ear cut (triangulate_ear_cut, n <= 64)" begin
    @testset "random simple polygons on both sides of the size cutoff" begin
        rng = Xoshiro(20261115)
        for n in (3, 4, 5, 6, 9, 17, 33, 63, 64, 65, 80), _ in 1:40
            P = _small_star(rng, n)
            @test _small_valid(P, triangulate_ear_cut(P))
            @test _small_valid(
                P,
                triangulate_ear_cut(P; kernel = ExactKernel());
                kernel = ExactKernel(),
            )
            @test _small_valid(P, triangulate_polygon(P; algorithm = :ear_cut))
        end
    end

    @testset "convex polygons are fanned from vertex 1" begin
        hexagon = [point(cos(2pi * i / 6), sin(2pi * i / 6)) for i in 0:5]
        res = triangulate_ear_cut(hexagon)
        @test res.triangles == [(1, i, i + 1) for i in 2:5]
        @test res.diagonals == [(1, i) for i in 3:5]
    end

    @testset "integer and rational coordinates" begin
        ints = [point(0, 0), point(4, 0), point(4, 4), point(2, 1), point(0, 4)]
        rats = [point(Rational(p.x), Rational(p.y)) for p in ints]
        @test _small_valid(ints, triangulate_ear_cut(ints))
        @test _small_valid(
            rats,
            triangulate_ear_cut(rats; kernel = ExactKernel());
            kernel = ExactKernel(),
        )
    end

    @testset "strictly simple quadrilateral from the 0.5.18 chain-swap fix" begin
        quad = [
            point(0.74, 0.21),
            point(-0.35, 0.56),
            point(-0.03, -0.56),
            point(-0.07, -0.42),
        ]
        @test _small_valid(quad, triangulate_ear_cut(quad))
    end

    @testset "validation" begin
        @test_throws ArgumentError triangulate_ear_cut([point(0.0, 0.0), point(1.0, 0.0)])
        cw = [point(0.0, 0.0), point(0.0, 1.0), point(1.0, 1.0), point(1.0, 0.0)]
        @test_throws ArgumentError triangulate_ear_cut(cw)
        consecutive = [point(0.0, 0.0), point(1.0, 0.0), point(1.0, 0.0), point(0.0, 1.0)]
        @test_throws ArgumentError triangulate_ear_cut(consecutive)
        closing = [point(0.0, 0.0), point(1.0, 0.0), point(0.0, 1.0), point(0.0, 0.0)]
        @test_throws ArgumentError triangulate_ear_cut(closing)
        @test_throws DomainError triangulate_ear_cut([
            point(NaN, 0.0),
            point(1.0, 0.0),
            point(0.0, 1.0),
        ])
    end

    @testset "non-simple and weakly simple input: valid result or ArgumentError" begin
        # A repeated vertex that is not consecutive (pinched bow tie), a zero-width spike, and
        # a vertex lying on another edge. The documented contract: either an ArgumentError or
        # non-clockwise triangles covering the polygon.
        pinch = [
            point(0.0, 0.0),
            point(2.0, 0.0),
            point(1.0, 1.0),
            point(2.0, 2.0),
            point(0.0, 2.0),
            point(1.0, 1.0),
        ]
        spike = [
            point(0.0, 0.0),
            point(2.0, 0.0),
            point(2.0, 1.0),
            point(1.0, 1.0),
            point(1.0, 3.0),
            point(1.0, 2.0),
            point(0.0, 1.0),
        ]
        on_edge = [
            point(0.0, 0.0),
            point(4.0, 0.0),
            point(4.0, 4.0),
            point(2.0, 0.0),
            point(0.0, 4.0),
        ]
        @test_throws ArgumentError triangulate_ear_cut(pinch)
        @test _small_valid(spike, triangulate_ear_cut(spike))
        @test _small_valid(on_edge, triangulate_ear_cut(on_edge))
        # A self-crossing quadrilateral (bow tie) has zero area and is rejected.
        bowtie = [point(0.0, 0.0), point(1.0, 1.0), point(1.0, 0.0), point(0.0, 1.0)]
        @test_throws ArgumentError triangulate_ear_cut(bowtie)
    end

    @testset "bounding-box rejection keeps closed-triangle semantics" begin
        inside = BasicCompGeometry._point_in_ccw_triangle
        a, b, c = point(0.0, 0.0), point(2.0, 0.0), point(0.0, 2.0)
        for k in (FloatingKernel(), FloatingKernel(Float32), ExactKernel())
            @test inside(k, point(0.5, 0.5), a, b, c)
            @test inside(k, point(1.0, 1.0), a, b, c)   # on the hypotenuse
            @test inside(k, point(2.0, 0.0), a, b, c)   # a corner
            @test !inside(k, point(1.5, 1.5), a, b, c)  # in the box, outside the triangle
            @test !inside(k, point(3.0, 0.0), a, b, c)  # outside the box
        end
        # Points just outside a corner. Float32 conversion rounds them onto the corner, so
        # the kernel's orientation tests see them as inside; the box test must agree.
        for q in (point(2.0 + 1e-12, 0.0), point(-1e-300, 0.0))
            @test inside(FloatingKernel(Float32), q, a, b, c)
            @test !inside(FloatingKernel(), q, a, b, c)
            @test !inside(ExactKernel(), q, a, b, c)
        end
    end
end
