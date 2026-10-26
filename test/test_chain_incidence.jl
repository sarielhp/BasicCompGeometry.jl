using Test
using Random
using BasicCompGeometry
import BasicCompGeometry as CG

@testset "Static RMQ" begin
    @test_throws BoundsError range_argmin(build_rmq(Int[]), 1, 1)
    for backend in (:blocked, :sparse), n in (1, 2, 5, 17, 70), trial = 1:3
        keys = rand(MersenneTwister(n + trial), 1:4, n)
        idx = build_rmq(keys; backend)
        for lo = 1:n, hi = lo:n
            m = minimum(keys[lo:hi])
            @test range_argmin(idx, lo, hi) == findlast(==(m), keys[1:hi]) |> (i -> max(i, lo))
        end
        @test_throws BoundsError range_argmin(idx, 0, 1)
        @test_throws BoundsError range_argmin(idx, 2, 1)
        @test_throws BoundsError range_argmin(idx, 1, n + 1)
    end
    # custom order, equivalent but unequal keys
    keys = [(3, :a), (1, :b), (2, :c), (1, :d), (5, :e)]
    idx = build_rmq(keys; less = (a, b) -> a[1] < b[1])
    @test range_argmin(idx, 1, 5) == 4
    @test range_argmin(idx, 1, 3) == 2
    @test range_argmin(idx, 3, 3) == 3
    # large n: all intervals of a long array spot-checked
    big = rand(MersenneTwister(1), 1:50, 3000)
    idx = build_rmq(big)
    rng = MersenneTwister(2)
    for _ = 1:2000
        lo = rand(rng, 1:3000); hi = rand(rng, lo:3000)
        m = minimum(@view big[lo:hi])
        @test range_argmin(idx, lo, hi) == lo - 1 + findlast(==(m), @view big[lo:hi])
    end
end

struct CountingVec{T} <: AbstractVector{T}
    v::Vector{T}
    reads::Base.RefValue{Int}
end
Base.size(c::CountingVec) = size(c.v)
Base.getindex(c::CountingVec, i::Int) = (c.reads[] += 1; c.v[i])

@testset "Chain index" begin
    occ(P, i) = P[mod(i, length(P))+1]
    function scan(P, lo, hi)
        n = length(P)
        vs = [(i, occ(P, i)) for i = lo:hi]
        lm = vs[1]
        for q in vs
            (q[2][1], q[2][2]) <= (lm[2][1], lm[2][2]) && (lm = q)
        end
        xs = [q[2][1] for q in vs]; ys = [q[2][2] for q in vs]
        lm, (xmin = minimum(xs), xmax = maximum(xs), ymin = minimum(ys), ymax = maximum(ys))
    end
    rng = MersenneTwister(7)
    for n = 1:9, trial = 1:4
        P = [(rand(rng, -2:2), rand(rng, -2:2) // rand(rng, 1:2)) for _ = 1:n]
        idx = build_chain_index(P)
        for lo = 0:n, hi = lo:n
            w = chain_window(idx, lo, hi)
            lm, bb = scan(P, lo, hi)
            @test nedges(w) == hi - lo
            @test CG.leftmost(w) == (occurrence = lm[1], point = lm[2])
            @test CG.bounds(w) == bb
        end
    end
    P = [(0, 0), (1, 0), (0, 0)]
    w = chain_window(build_chain_index(P), 0, 3)
    @test CG.leftmost(w).occurrence == 3
    @test CG.leftmost(chain_window(build_chain_index(P), 3, 3)).occurrence == 3
    @test CG.vertex(build_chain_index(P), 3) == (0, 0)

    Q = [(0, 0), (4, 2), (2, 6)]
    w = chain_window(build_chain_index(Q), 0, 3)
    @test CG.edge(w, 1) == ((0, 0), (4, 2))
    @test CG.edge(w, 3; direction = -1) == ((0, 0), (2, 6))
    @test CG.edge(w, 3) == ((2, 6), (0, 0))
    @test_throws BoundsError CG.edge(w, 0)
    @test_throws ArgumentError CG.edge(w, 1; direction = 0)
    @test CG.point_at(w, 1 // 2, 0) == (2 // 1, 1 // 1)
    @test CG.point_at(w, 1, 0) == CG.point_at(w, 1, 1) == (4 // 1, 2 // 1)
    @test CG.point_at(w, 3, 3) == (0 // 1, 0 // 1)
    @test CG.point_at(w, 5 // 2) == (1 // 1, 3 // 1)
    @test_throws ArgumentError CG.point_at(w, 3 // 2, 0)
    @test_throws BoundsError CG.point_at(w, 3, 4)
    @test_throws BoundsError CG.point_at(w, 3, -1)
    @test_throws BoundsError chain_window(build_chain_index(Q), 2, 4)
    @test_throws BoundsError chain_window(build_chain_index(Q), 2, 1)
    @test_throws ArgumentError build_chain_index(Tuple{Int,Int}[])
    @test_throws ArgumentError build_chain_index([(1.0, 2)])
    # overflow safety
    big = typemax(Int64)
    w2 = chain_window(build_chain_index([(big, big), (-big, -big)]), 0, 1)
    @test CG.point_at(w2, 1 // 2, 0) == (0 // 1, 0 // 1)
    # SVector points
    wp = chain_window(build_chain_index([Point(0, 0), Point(2, 2)]), 0, 1)
    @test CG.point_at(wp, 1 // 2, 0) == Point(1 // 1, 1 // 1)

    # no per-query scans
    N = 20000
    cv = CountingVec([(rand(rng, 1:100), rand(rng, 1:100)) for _ = 1:N], Ref(0))
    idx = build_chain_index(cv)
    cv.reads[] = 0
    w = chain_window(idx, 3, N - 5)
    CG.bounds(w); CG.leftmost(w)
    @test cv.reads[] < 200
end

@testset "Vertical rays" begin
    R = Rational{BigInt}
    h = vertical_ray_segment((0, 0), 1, (-1, 1), (1, 3))
    @test h isa PointHit && h.t == 2 && h.point == (0 // 1, 2 // 1)
    @test vertical_ray_segment((0, 0), -1, (-1, 1), (1, 3)) isa NoHit
    @test vertical_ray_segment((0, 0), 1, (1, 3), (-1, 1)).t == 2
    @test vertical_ray_segment((0, 2), 1, (-1, 1), (1, 3)).t == 0
    @test vertical_ray_segment((0, 0), 1, (0, 1), (5, 5)).t == 1       # endpoint
    @test vertical_ray_segment((6, 0), 1, (0, 1), (5, 5)) isa NoHit
    @test vertical_ray_segment((0, 0), 1, (0, 0), (0, 0)) == PointHit((0 // 1, 0 // 1), R(0))
    @test vertical_ray_segment((1, 0), 1, (0, 0), (0, 0)) isa NoHit
    # vertical overlaps
    @test vertical_ray_segment((0, 0), 1, (0, 2), (0, 5)) == OverlapHit(R(2), R(5))
    @test vertical_ray_segment((0, 0), 1, (0, 5), (0, 2)) == OverlapHit(R(2), R(5))
    @test vertical_ray_segment((0, 0), 1, (0, -3), (0, 4)) == OverlapHit(R(0), R(4))
    @test vertical_ray_segment((0, 0), -1, (0, -3), (0, 4)) == OverlapHit(R(0), R(3))
    @test vertical_ray_segment((0, 0), 1, (0, -3), (0, -1)) isa NoHit
    @test vertical_ray_segment((0, 0), 1, (0, -3), (0, 0)) == PointHit((0 // 1, 0 // 1), R(0))
    @test vertical_ray_segment((0, 0), 1, (1, 0), (1, 5)) isa NoHit
    @test_throws ArgumentError vertical_ray_segment((0, 0), 0, (0, 0), (1, 1))
    @test_throws ArgumentError vertical_ray_segment((0.0, 0), 1, (0, 0), (1, 1))
    # huge rationals
    b = big(10)^40
    h = vertical_ray_segment((1 // b, 0), 1, (0, 0), (1, 1))
    @test h.t == 1 // b
    # collection scan
    segs = [((0, 5), (2, 5), :far), ((0, 2), (2, 2), "near"), ((1, 1), (1, 1) .+ (1, 1), :third), ((0, 2), (4, 2), :tie), ((-1, 0), (1, 0), :origin)]
    r = first_positive_vertical_contact(segs, (1, 0), 1)
    @test r.carrier == :third && r.t == 1 && r.record_index == 3
    r = first_positive_vertical_contact(segs[[1, 2, 4]], (1, 0), 1)
    @test r.carrier == "near" && r.record_index == 2 && r.point == (1 // 1, 2 // 1)
    @test first_positive_vertical_contact(segs, (1, 0), -1) === nothing
    @test first_positive_vertical_contact([((0, 0), (2, 0), 1)], (1, 0), 1) === nothing
    @test_throws ArgumentError first_positive_vertical_contact([((0, 0), (0, 2), 1)], (0, -1), 1)
end

@testset "Angular sectors" begin
    ccw(v, u, d; kw...) = open_sector_contains(v, u, d; kw...)
    # acute sector x-axis to y-axis
    @test ccw((1, 0), (0, 1), (1, 1))
    @test !ccw((1, 0), (0, 1), (-1, -1))
    @test !ccw((1, 0), (0, 1), (1, 0))
    @test !ccw((1, 0), (0, 1), (0, 1))
    @test ccw((1, 0), (0, 1), (5, 7))
    # reflex: from y-axis ccw to x-axis passes through -x, -y
    @test ccw((0, 1), (1, 0), (-1, 0))
    @test ccw((0, 1), (1, 0), (-1, -1))
    @test !ccw((0, 1), (1, 0), (1, 1))
    # straight: upper half plane
    @test ccw((1, 0), (-1, 0), (0, 1))
    @test !ccw((1, 0), (-1, 0), (0, -1))
    @test !ccw((1, 0), (-1, 0), (3, 0))
    @test !ccw((1, 0), (-1, 0), (-3, 0))
    # slit tip
    @test ccw((1, 0), (2, 0), (-1, 0))
    @test ccw((1, 0), (2, 0), (0, 1))
    @test !ccw((1, 0), (2, 0), (3, 0))
    @test !ccw((1, 0), (2, 0), (1, 0); coincident = :empty)
    @test !ccw((1, 0), (2, 0), (-1, 0); coincident = :empty)
    @test_throws ArgumentError ccw((1, 0), (2, 0), (1, 1); coincident = :bogus)
    @test_throws ArgumentError ccw((0, 0), (1, 0), (1, 1))
    @test_throws ArgumentError ccw((1, 0), (1, 1), (0, 0))
    @test ccw(CG.Point(1, 0), CG.Point(0, 1), CG.Point(1, 1))
    # scaling invariance, exhaustive on small grid
    vs = [(x, y) for x = -2:2, y = -2:2 if (x, y) != (0, 0)]
    for v in vs, u in vs, d in vs
        @test ccw(v, u, d) == ccw(3 .* v, 5 .* u, 7 .* d)
    end

    # symbolic incoming vs exact rational perturbation
    function perturbed(v, u, s, eps)
        ccw(v, u, (-eps, -s // 1))
    end
    for v in vs, u in vs, s in (-1, 1), coin in (:full, :empty)
        expected = perturbed(v, u, s, 1 // 10^6) # grid thresholds are >> 1e-6
        # same-direction/empty handled by shared semantic
        @test incoming_sector_contains(v, u, s; coincident = coin) ==
              open_sector_contains(v, u, (-1 // 10^6, -s // 1); coincident = coin)
        coin == :full && @test incoming_sector_contains(v, u, s) == expected
    end
    @test_throws ArgumentError incoming_sector_contains((1, 0), (0, 1), 2)
    @test_throws ArgumentError incoming_sector_contains((0, 0), (0, 1), 1)
    # slit tip contains the incoming direction
    @test incoming_sector_contains((0, 1), (0, 3), 1)
end
