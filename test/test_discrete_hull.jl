using Test
using Random
using BasicCompGeometry
import BasicCompGeometry as CG

# Independent verification. H equals ConvexHull(C ∩ Z²) iff (a) every vertex of
# H is a lattice point of C, (b) H is strictly convex and counterclockwise (no
# collinear vertices), and (c) every lattice point of C lies in H. The brute
# force enumerates C ∩ Z² from exact membership tests only.
lcross(u, v) = Int128(u[1]) * v[2] - Int128(u[2]) * v[1]

function is_discrete_hull(H, pts, inC)
    isempty(pts) && return isempty(H)
    all(inC, H) || return false
    n = length(H)
    allunique(H) || return false
    if n >= 3
        for i = 1:n
            a, b, c = H[i], H[mod1(i + 1, n)], H[mod1(i + 2, n)]
            lcross(b - a, c - b) > 0 || return false
        end
        return all(p -> all(lcross(H[mod1(i + 1, n)] - H[i], p - H[i]) >= 0 for i = 1:n), pts)
    elseif n == 2
        a, b = H
        return all(p -> lcross(b - a, p - a) == 0 && 0 <= (p - a)' * (b - a) <= (b - a)' * (b - a),
                   pts)
    else
        return all(==(H[1]), pts)
    end
end

lowest_leftmost(H) = isempty(H) || H[1] == H[argmin([(v[2], v[1]) for v in H])]

function box_points(inC, x0, x1, y0, y1)
    return [Point2I(x, y) for x = x0:x1 for y = y0:y1 if inC(Point2I(x, y))]
end

function disk_points(D)
    c, r = CG.center(D), CG.radius(D)
    return box_points(p -> p ∈ D, floor(Int, c[1] - r) - 1, ceil(Int, c[1] + r) + 1,
                      floor(Int, c[2] - r) - 1, ceil(Int, c[2] + r) + 1)
end

# Column extremes of an integer-center, integer-radius disk by exact isqrt
# (independent of chord_range).
function disk_column_extremes(cx::Int, cy::Int, r::Int)
    pts = Point2I[]
    for x = cx-r:cx+r
        h = isqrt(r^2 - (x - cx)^2)
        push!(pts, Point2I(x, cy - h), Point2I(x, cy + h))
    end
    return pts
end

@testset "Discrete hull (Har-Peled 1998)" begin
    @testset "chord_range is exact" begin
        rng = MersenneTwister(7)
        bodies = [LatticeDisk((rand(rng) * 6 - 3, rand(rng) * 6 - 3), 1 + 9rand(rng)) for _ = 1:30]
        append!(bodies, [LatticeDisk((0, 0), 5), LatticeDisk((1 // 2, 0), 5 // 2)])
        for D in bodies, _ = 1:40
            o = Point2I(rand(rng, -12:12), rand(rng, -12:12))
            d = Point2I(rand(rng, -5:5), rand(rng, -5:5))
            d == Point2I(0, 0) && continue
            r = chord_range(D, o, d)
            ts = [t for t = -40:40 if o + t * d ∈ D]
            if r === nothing
                @test isempty(ts)
            else
                @test ts == [t for t = r[1]:r[2]]
            end
        end
        # Line tangent to x² + y² = 25 at the lattice point (3, 4): chord is {0}.
        @test chord_range(LatticeDisk((0, 0), 5), Point2I(3, 4), Point2I(4, -3)) == (0, 0)
        @test chord_range(LatticeDisk((0, 0), 5), Point2I(6, 0), Point2I(0, 1)) === nothing
    end

    @testset "find_edge_direction: slope convergents" begin
        # Half-plane y ≤ (8/5) x, clipped to a large box: from the origin, the
        # first vector clockwise of (0,1) in the quadrant cone is (5, 8) - the
        # last convergent of 8/5 (Figure 2.2 of the paper).
        P = ConvexPolygonBody([(0, 0), (1000, 0), (1000, 1600)])
        @test find_edge_direction(P, Point2I(0, 0), Point2I(0, 1), Point2I(-1, 0)) === nothing
        e = find_edge_direction(P, Point2I(0, 0), Point2I(1, 0), Point2I(0, 1))
        @test e == Point2I(1, 0)
        Q = ConvexPolygonBody([(0, 0), (1000, 1600), (0, 1000)])
        @test find_edge_direction(Q, Point2I(0, 0), Point2I(1, 0), Point2I(0, 1)) == Point2I(5, 8)
        # Slopes just below 355/113 and other near-rationals: compare with the
        # brute-force smallest y/x over the triangle's lattice points.
        for (num, den, eps) in ((355, 113, 1 // 7), (355, 113, 0), (13, 8, 1 // 1000), (1, 1, 1 // 3))
            k = 600 ÷ den
            R = ConvexPolygonBody([(0, 0), (den * k, num * k - eps), (0, num * k + 50)])
            best = nothing
            for x = 1:den*k, y = 0:num*k+50
                q = Point2I(x, y)
                q ∈ R || continue
                (best === nothing || lcross(q, best) > 0) && (best = q)
            end
            best = best .÷ gcd(best[1], best[2])
            @test find_edge_direction(R, Point2I(0, 0), Point2I(1, 0), Point2I(0, 1)) == best
        end
    end

    @testset "degenerate disks" begin
        @test discrete_hull(LatticeDisk((1 // 2, 1 // 2), 1 // 2)) == Point2I[]
        @test discrete_hull(LatticeDisk((3, -2), 0)) == [Point2I(3, -2)]
        @test discrete_hull(LatticeDisk((3, -2), 1 // 2)) == [Point2I(3, -2)]
        @test discrete_hull(LatticeDisk((1 // 2, 0), 1 // 2)) == [Point2I(0, 0), Point2I(1, 0)]
        @test discrete_hull(LatticeDisk((0, 0), 1)) ==
              [Point2I(0, -1), Point2I(1, 0), Point2I(0, 1), Point2I(-1, 0)]
        @test discrete_hull(LatticeDisk((0, 0), 5)) == [
            Point2I(0, -5), Point2I(3, -4), Point2I(4, -3), Point2I(5, 0),
            Point2I(4, 3), Point2I(3, 4), Point2I(0, 5), Point2I(-3, 4),
            Point2I(-4, 3), Point2I(-5, 0), Point2I(-4, -3), Point2I(-3, -4),
        ]
    end

    @testset "random disks vs brute force" begin
        rng = MersenneTwister(1998)
        for trial = 1:400
            r = trial <= 100 ? 3rand(rng) : 2 + 60rand(rng)
            c = (100randn(rng), 100randn(rng))
            D = trial % 4 == 0 ? LatticeDisk(round.(Int, c), round(Int, r)) : LatticeDisk(c, r)
            H = discrete_hull(D)
            @test is_discrete_hull(H, disk_points(D), p -> p ∈ D)
            @test lowest_leftmost(H)
        end
    end

    @testset "rational convex polygons (thin, non-connected) vs brute force" begin
        rng = MersenneTwister(42)
        ratl(x) = round(Int, 64x) // 64
        for trial = 1:300
            # Thin random triangles stress the continued-fraction search.
            a = (ratl(30randn(rng)), ratl(30randn(rng)))
            b = (ratl(30randn(rng)), ratl(30randn(rng)))
            w = trial % 3 == 0 ? 0.05 : 3.0
            m = ((a[1] + b[1]) / 2 + ratl(w * randn(rng)), (a[2] + b[2]) / 2 + ratl(w * randn(rng)))
            o = (b[1] - a[1]) * (m[2] - a[2]) - (b[2] - a[2]) * (m[1] - a[1])
            o == 0 && continue
            P = ConvexPolygonBody(o > 0 ? [a, b, m] : [a, m, b])
            xs = [v[1] for v in CG.vertices(P)]
            ys = [v[2] for v in CG.vertices(P)]
            pts = box_points(p -> p ∈ P, floor(Int, minimum(xs)), ceil(Int, maximum(xs)),
                             floor(Int, minimum(ys)), ceil(Int, maximum(ys)))
            H = discrete_hull(P)
            @test is_discrete_hull(H, pts, p -> p ∈ P)
            @test lowest_leftmost(H)
            isempty(H) || @test discrete_hull(P; start = H[end]) == H
        end
    end

    @testset "generic float oracle (ellipses)" begin
        rng = MersenneTwister(3)
        for _ = 1:100
            cx, cy = 10randn(rng), 10randn(rng)
            A, B, θ = 1 + 20rand(rng), 1 + 20rand(rng), π * rand(rng)
            cθ, sθ = cos(θ), sin(θ)
            # Ellipse ((u⋅e1)/A)² + ((u⋅e2)/B)² ≤ 1, u = x - c.
            function chord(o, d)
                u = (o[1] - cx, o[2] - cy)
                f(v) = ((v[1] * cθ + v[2] * sθ) / A, (-v[1] * sθ + v[2] * cθ) / B)
                U, Dd = f(u), f(d)
                a = Dd[1]^2 + Dd[2]^2
                b = U[1] * Dd[1] + U[2] * Dd[2]
                c = U[1]^2 + U[2]^2 - 1
                disc = b^2 - a * c
                disc < 0 && return nothing
                s = sqrt(disc)
                return ((-b - s) / a, (-b + s) / a)
            end
            E = FloatChordBody(chord)
            inE(p) = (r = chord(Point2F(p), Point2F(1, 0)); r !== nothing && r[1] <= 0 <= r[2])
            seed = Point2I(round(Int, cx), round(Int, cy))
            H = discrete_hull(E; seed = seed)
            pts = box_points(inE, floor(Int, cx - 22), ceil(Int, cx + 22),
                             floor(Int, cy - 22), ceil(Int, cy + 22))
            @test is_discrete_hull(H, pts, inE)
        end
        @test_throws ArgumentError discrete_hull(FloatChordBody((o, d) -> nothing))
    end

    @testset "large disk vs independent column hull" begin
        for (cx, cy, r) in ((0, 0, 100_000), (17, -5, 99_991), (3, 4, 31_623))
            D = LatticeDisk((cx, cy), r)
            H = discrete_hull(D)
            @test H == lattice_convex_hull(disk_column_extremes(cx, cy, r))
            @test all(p -> p ∈ D, H)
            # Balog–Bárány: |H| = Θ(r^(2/3)); the paper measures ≈ 3.45 r^(2/3).
            @test 3.0 < length(H) / r^(2 / 3) < 4.0
        end
        D = LatticeDisk((0.3141, 0.2718), 250_000.5)
        @test discrete_hull(D) == discrete_hull_columns(D)
    end

    @testset "lattice_convex_hull" begin
        rng = MersenneTwister(5)
        for _ = 1:200
            pts = [Point2I(rand(rng, -9:9), rand(rng, -9:9)) for _ = 1:rand(rng, 1:30)]
            H = lattice_convex_hull(pts)
            @test is_discrete_hull(H, pts, p -> p in pts)
        end
    end
end
