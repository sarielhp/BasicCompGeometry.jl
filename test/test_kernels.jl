using Test
using Random
using BasicCompGeometry
import BasicCompGeometry: Points, compare_right_arc_to_point, compare_right_arcs

# Stateful external kernel proves that wrappers preserve the same object.
mutable struct CountingTestKernel <: AbstractKernel
    turns::Int
    comparisons::Int
    distances::Int
end
function BasicCompGeometry.turn_sign(k::CountingTestKernel, p::Point{2}, q::Point{2}, r::Point{2})
    k.turns += 1
    return turn_sign(ExactKernel(), p, q, r)
end
function BasicCompGeometry.compare_xy(k::CountingTestKernel, p::Point{2}, q::Point{2})
    k.comparisons += 1
    return compare_xy(ExactKernel(), p, q)
end
function BasicCompGeometry.compare_distance(k::CountingTestKernel, p::Point{D},
    q::Point{D}, radius::Real) where {D}
    k.distances += 1
    return compare_distance(ExactKernel(), p, q, radius)
end
struct UnsupportedTestKernel <: AbstractKernel end

# Independent reference: choose dyadic rational enclosing intervals for each
# radical using integer square roots. Use only strict, certified separation.
function certified_radical_order(d, a, b)
    for bits in (16, 64, 256)
        scale = big(1) << bits
        lower_a, lower_b = isqrt(big(a) * scale^2), isqrt(big(b) * scale^2)
        lower = big(d) * scale + lower_a - lower_b - 1
        upper = big(d) * scale + lower_a + 1 - lower_b
        lower > 0 && return 1
        upper < 0 && return -1
    end
    return nothing
end

@testset "Pluggable predicate kernels" begin
    exact = ExactKernel()
    p, q, r = point(0, 0), point(1, 0), point(0, 1)
    @test (@inferred turn_sign(p, q, r)) === 1.0
    @test (@inferred turn_sign(p, q, r; kernel=FloatingKernel(Float32))) === 1.0f0
    @test (@inferred turn_sign(p, q, r; kernel=exact)) == big(1)
    @test is_left_turn(p, q, r; kernel=exact)
    @test is_right_turn(p, r, q; kernel=exact)
    @test is_left_eq_turn(p, p, q; kernel=exact)
    @test is_right_eq_turn(p, p, q; kernel=exact)
    @test is_collinear(p, p, q; kernel=exact)
    @test compare_xy(p, q) == -1

    # Cancellation of two exactly representable but large products.
    n = 2^27
    a, b, c = point(0.0, 0.0), point(Float64(n), Float64(n - 1)),
              point(Float64(n + 1), Float64(n))
    @test turn_sign(a, b, c) == 0.0
    @test turn_sign(a, b, c; kernel=exact) == 1
    @test length(convex_hull([a, b, c]; kernel=exact)) == 3
    @test length(convex_hull([a, b, c])) == 2

    # Widen before subtraction, including extremes of the machine integer type.
    low, high = typemin(Int), typemax(Int)
    expected = (big(high) - low)^2
    @test turn_sign(point(low, low), point(high, low), point(low, high); kernel=exact) == expected
    @test compare_squared_distance(point(low, 0), point(high, 0), expected; kernel=exact) == 0
    @test compare_distance(point(low, 0), point(high, 0), big(high) - low; kernel=exact) == 0
    @test compare_distance(point(0, 0), point(3, 4), 5; kernel=exact) == 0
    @test compare_distance(point(0, 0), point(3, 4), 4; kernel=exact) == 1
    @test compare_distance(point(0, 0), point(3, 4), 6; kernel=exact) == -1
    @test is_inside(point(3, 4), Circle(point(0, 0), 5); kernel=exact)

    # Finite floats are interpreted as their exact stored binary values.
    @test compare_squared_distance(point(0.1, 0.0), point(0.0, 0.0), 1//100; kernel=exact) == 1
    @test compare_distance(point(1//10, 0//1), point(0, 0), 1//10; kernel=exact) == 0
    @test turn_sign(point(0//1, 0//1), point(1//3, 0//1), point(0//1, 1//7); kernel=exact) == 1//21
    @test compare_distance(point(-0.0, 0.0), point(0.0, 0.0), 0; kernel=exact) == 0

    # Radical equality, domain endpoints, zero radius, and huge coordinates.
    @test compare_right_arcs(point(0, 0), point(1, 3), 0, 5; kernel=exact) == 0
    @test compare_right_arc_to_point(point(0, 0), point(3, 4), 5; kernel=exact) == 0
    @test compare_right_arc_to_point(point(0, 0), point(-3, 4), 5; kernel=exact) == 1
    @test compare_right_arcs(point(0, 0), point(0, 0), 5, 5; kernel=exact) == 0
    @test compare_right_arcs(point(0, 0), point(1, 0), 0, 0; kernel=exact) == -1
    @test compare_right_arcs(point(low, 0), point(high, 0), 0, 1; kernel=exact) == -1
    @test compare_right_arc_to_point(point(1//3, 0//1), point(2//3, 0//1), 1//3; kernel=exact) == 0
    @test compare_right_arcs(point(0, 0), point(1, 0), 0, 1) == -1

    for kernel in (FloatingKernel(), exact)
        @test_throws DomainError compare_distance(p, q, -1; kernel=kernel)
        @test_throws DomainError compare_right_arcs(p, q, 2, 1; kernel=kernel)
        @test_throws DomainError turn_sign(point(Inf, 0.0), p, q; kernel=kernel)
        @test_throws DomainError compare_distance(point(NaN, 0.0), p, 1; kernel=kernel)
        @test_throws DomainError compare_distance(point(1//0, 0//1), p, 1; kernel=kernel)
    end
    @test_throws MethodError is_left_turn(p, q, r; kernel=UnsupportedTestKernel())

    # Algorithm entry and point-sequence wrapper both forward the kernel.
    tracked = CountingTestKernel(0, 0, 0)
    hull = convex_hull(PntSeq([a, b, c]); kernel=tracked)
    @test length(hull) == 3
    @test tracked.turns > 0 && tracked.comparisons > 0
    @test is_inside(point(3, 4), Circle(point(0, 0), 5); kernel=tracked)
    @test tracked.distances == 1
    @test eltype(Points(hull)) == Point{2,Float64}

    rng = MersenneTwister(714)
    checked = 0
    for _ in 1:2_000
        radius = rand(rng, 1:10_000)
        dy1, dy2 = rand(rng, -radius:radius), rand(rng, -radius:radius)
        dx = rand(rng, -2radius:2radius)
        a2, b2 = big(radius)^2 - big(dy1)^2, big(radius)^2 - big(dy2)^2
        reference = certified_radical_order(dx, a2, b2)
        isnothing(reference) && continue
        first, second = point(dx, dy1), point(0, dy2)
        @test compare_right_arcs(first, second, 0, radius; kernel=exact) == reference
        @test compare_right_arcs(second, first, 0, radius; kernel=exact) == -reference
        checked += 1
    end
    @test checked > 1_900
end
