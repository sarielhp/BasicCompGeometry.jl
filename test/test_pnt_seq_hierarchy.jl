using BasicCompGeometry
import BasicCompGeometry: hausdorff_dist_subseg, dist_point_segment, cardin
using BasicCompGeometry.PntSeqHierarchy
using Random
using Test

const H = BasicCompGeometry.PntSeqHierarchy

function walk(n; seed = 1)
    rng = MersenneTwister(seed)
    cur = point(0.0, 0.0)
    pts = Point2F[]
    for _ = 1:n
        cur = cur + point(rand(rng) - 0.5, rand(rng) - 0.5) * (1 + 3 * rand(rng))
        push!(pts, cur)
    end
    return PntSeq(pts)
end

# Largest distance from a vertex of P to the polyline through P[idx].
function one_sided(P, idx)
    m = 0.0
    for t = 1:(length(idx)-1), k = idx[t]:idx[t+1]
        m = max(m, dist_point_segment(P[k], P[idx[t]], P[idx[t+1]]))
    end
    return m
end

@testset "PntSeqHierarchy" begin
    @testset "Palette: midpoint rule" begin
        P = walk(60)
        pal = palette(P, hausdorff_dist_subseg)
        @test pal isa AbstractVector{Float64}
        @test length(pal) == 60
        @test all(>=(0.0), pal)
        @test pal[1] == 0.0 && pal[60] == 0.0

        # Hand computed recursion for the midpoint rule.
        expect = zeros(60)
        function rec(i, j)
            m = (i + j) >> 1
            (m == i || m == j) && return 0.0
            l, r = rec(i, m), rec(m, j)
            expect[m] = max(l, r)
            return hausdorff_dist_subseg(P, i:j)
        end
        rec(1, 60)
        @test collect(pal) == expect

        # A plain vector is an implicit midpoint palette.
        a, b = palette_extract(P, collect(pal), 0.5), palette_extract(P, pal, 0.5)
        @test a[2] == b[2] && a[3] == b[3] && a[1].pnts == b[1].pnts
    end

    @testset "Extraction" begin
        for rule in (MidpointSplit(), FarthestSplit())
            P = walk(80; seed = 7)
            pal = palette(P, hausdorff_dist_subseg; rule = rule)
            # No vertex is below a zero threshold: everything is kept.
            Q, idx, _ = palette_extract(P, pal, 0.0)
            @test idx == 1:80
            @test collect(Q) == collect(P)
            # Everything is below an enormous threshold: only the root split remains.
            Q, idx, _ = palette_extract(P, pal, 1e9)
            @test length(idx) == 3 && idx[1] == 1 && idx[end] == 80
            # Intermediate thresholds give increasing indices, a subset of P.
            for w in (0.2, 0.5, 1.0, 2.0)
                Q, idx, _ = palette_extract(P, pal, w)
                @test issorted(idx; lt = <=)
                @test collect(Q) == [P[i] for i in idx]
            end
            # Finer thresholds never keep fewer vertices.
            sizes =
                [length(palette_extract(P, pal, w)[2]) for w in (4.0, 2.0, 1.0, 0.5, 0.1)]
            @test issorted(sizes)
        end
        P = walk(10)
        pal = palette(P, hausdorff_dist_subseg)
        @test palette_extract(walk(4), palette(walk(4), hausdorff_dist_subseg), 1.0)[2] ==
              1:4
        @test palette_level(P, pal, 0)[2] == [1, 10]
        @test palette_level(P, pal, 1)[2] == [1, 5, 10]
        @test palette_level(P, pal, 100)[2] == 1:10
    end

    @testset "Farthest split is Douglas-Peucker" begin
        P = PntSeq([
            point(0.0, 0.0),
            point(1.0, 0.1),
            point(2.0, 3.0),
            point(3.0, 0.2),
            point(4.0, 0.0),
        ])
        pal = palette(P, hausdorff_dist_subseg; rule = FarthestSplit())
        @test palette_level(P, pal, 1)[2] == [1, 3, 5]
        @test pal[3] ≥ pal[2] && pal[3] ≥ pal[4]
    end

    @testset "Degenerate curves" begin
        for n = 0:3
            P = n == 0 ? PntSeq2F() : walk(n)
            @test length(palette(P, hausdorff_dist_subseg)) == n
            @test all(iszero, palette(P, hausdorff_dist_subseg; rule = FarthestSplit()))
        end
    end

    @testset "Hierarchy" begin
        P = walk(300; seed = 3)
        h = hierarchy_init(P)
        @test h isa SimplificationHierarchy
        @test length(h.polys) == 1 && cardin(h.polys[1]) == 2
        @test h.len ≈ geom_length(P)
        hierarchy_build!(h)
        @test issorted(h.widths; rev = true)
        @test issorted(map(cardin, h.polys))
        for (Q, w) in zip(h.polys, h.widths)
            idx = [findfirst(==(q), collect(P)) for q in Q]
            @test one_sided(P, idx) ≤ w + 1e-12
        end

        for w in (3.0, 0.7, 0.1)
            Q, wq = hierarchy_approx(h, w)
            idx = [findfirst(==(q), collect(P)) for q in Q]
            @test one_sided(P, idx) ≤ wq + 1e-12
        end
        @test issorted(h.widths; rev = true)
        @test hierarchy_find([4.0, 1.0, 0.25], 1.5, 4.0) == 2
        @test hierarchy_find([4.0, 1.0, 0.25], 0.1, 4.0) == -1
        @test hierarchy_find([4.0, 1.0, 0.25], 5.0, 4.0) == 1
        @test hierarchy_find([4.0, 1.0, 0.25], 0.3, 1.1) == -1   # 0.25 is too narrow for 0.3
        @test hierarchy_find([4.0, 1.0, 0.25], 0.3, 1.3) == 3
        @test hierarchy_find([4.0, 1.0], 2.5, 4.0) == 2          # first level of width <= w

        redirect_stdout(devnull) do
            hierarchy_print(h)
        end
    end

    @testset "Push keeps levels sorted and drops dominated ones" begin
        P = walk(20)
        h = hierarchy_init(P)
        h.widths[1] = 10.0
        seq(n) = PntSeq(P.pnts[1:n])
        hierarchy_push!(h, seq(5), 6.0)
        @test h.widths == [10.0, 6.0]
        hierarchy_push!(h, seq(9), 8.0)          # not dominated: sorted into place
        @test h.widths == [10.0, 8.0, 6.0]
        @test map(cardin, h.polys) == [2, 9, 5]
        hierarchy_push!(h, seq(3), 5.0)          # last level is narrower-or-equal and larger: dropped
        @test h.widths == [10.0, 8.0, 5.0]
        @test map(cardin, h.polys) == [2, 9, 3]
        @test_throws AssertionError hierarchy_push!(h, PntSeq([P[1]]), 1.0)
    end

    @testset "Custom measure" begin
        P = walk(100; seed = 11)
        # Measure the number of skipped vertices instead of a distance.
        err(Q, rng) = float(length(rng) - 2)
        simplify(Q, w) = (
            idx = collect(1:max(1, floor(Int, w)+1):length(Q));
            idx[end] == length(Q) || push!(idx, length(Q));
            (PntSeq([Q[i] for i in idx]), idx)
        )
        certify(Q, R, idx) = maximum(diff(idx)) - 1.0
        m = HierarchyMeasure(err, simplify, certify, MidpointSplit())
        h = hierarchy_init(P, m)
        @test h.widths == [98.0]
        hierarchy_build!(h)
        @test issorted(h.widths; rev = true)
        @test hierarchy_approx(h, 20.0)[2] ≤ 100.0
    end
    @testset "Float32 curves and integer rejection" begin
        P32 = PntSeq([Point{2,Float32}(Float32(x), Float32(sin(x))) for x = 0:0.05:6])
        h = hierarchy_init(P32)
        hierarchy_build!(h)
        Q, w = hierarchy_approx(h, 0.3)
        @test Q isa PntSeq{2,Float32}
        @test w >= 0
        Pint = PntSeq([Point{2,Int}(i, i^2 % 7) for i = 1:20])
        @test_throws MethodError palette(Pint, hausdorff_dist_subseg)
        @test_throws MethodError hierarchy_init(Pint)
    end

    @testset "Build stops when no progress is made" begin
        line = PntSeq([point(Float64(i), 2.0 * i) for i = 1:12])   # collinear: all widths zero
        h = hierarchy_init(line)
        hierarchy_build!(h)
        @test length(h.polys) <= 3
        @test issorted(h.widths; rev = true)
    end

    @testset "Concurrent approx and push" begin
        P = walk(400; seed = 5)
        h = hierarchy_init(P)
        hierarchy_build!(h)
        tasks = [
            Threads.@spawn(hierarchy_approx(h, w)) for w in repeat([2.0, 0.9, 0.4, 0.15], 8)
        ]
        results = fetch.(tasks)
        @test all(r -> r[2] >= 0 && cardin(r[1]) >= 2, results)
        @test issorted(h.widths; rev = true)
        @test length(h.polys) == length(h.widths)
    end
end
