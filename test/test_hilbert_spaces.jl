using Test
using LinearAlgebra
using Random
using BasicCompGeometry
using BasicCompGeometry.MetricSpace: greedy_permutation_naive, PointsSpace
using BasicCompGeometry:
    VDiff,
    cache_stats,
    inner_diff,
    is_generator,
    is_vertex_curve,
    line_params,
    local_frame,
    split_edges,
    vdist2,
    vscale

# The coordinates of a virtual point of a space built from the coordinate vectors X.
realize(X, p::VPoint) = (1 - p.t) * X[p.a] + p.t * X[p.b]

# Euclidean line parameters of q against the line through a and b.
function line_params_xyz(q, a, b)
    d = b - a
    t0 = dot(q - a, d) / dot(d, d)
    return t0, sum(abs2, a + t0 * d - q)
end

# A random virtual point on a random pair of the given vertices.
function rand_vpoint(rng, S, verts)
    a, b = rand(rng, verts), rand(rng, verts)
    return VPoint(S, a, b, rand(rng, (0.0, 1.0, rand(rng))))
end

@testset "Hilbert spaces and virtual points" begin
    @testset "VPoint: canonical form, equality, validation" begin
        S = HilbertSpace([rand(3) for _ in 1:10])
        @test VPoint(S, 7, 3, 0.75) == VPoint(S, 3, 7, 0.25)
        @test hash(VPoint(S, 7, 3, 0.75)) == hash(VPoint(S, 3, 7, 0.25))
        @test VPoint(S, 3, 7, 0.0) == vpoint(S, 3)
        @test VPoint(S, 3, 7, 1.0) == vpoint(S, 7)
        @test VPoint(S, 5, 5, 0.4) == vpoint(S, 5)
        @test VPoint(S, 3, 7, 1.0).a == 7 && VPoint(S, 3, 7, 1.0).t == 0
        @test sizeof(VPoint(S, 1, 2, 0.5)) == 24
        @test isbitstype(VPoint{typeof(S)}) == false
        @test_throws ArgumentError VPoint(S, 1, 2, 1.5)
        @test_throws ArgumentError VPoint(S, 1, 2, NaN)
        @test_throws BoundsError VPoint(S, 0, 2, 0.5)
        @test_throws BoundsError vpoint(S, 11)
        T = HilbertSpace([rand(3) for _ in 1:10])
        @test VPoint(S, 1, 2, 0.5) != VPoint(T, 1, 2, 0.5)
        @test_throws ArgumentError dist(vpoint(S, 1), vpoint(T, 2))
        x = VPoint(S, 1, 2, 0.5)
        @test deepcopy(x).space === S
        @test deepcopy([x, x])[1].space === S
        @test !isnan(x) && isfinite(x)
        @test sprint(show, x) == "VPoint(1, 2, 0.5)"
        @test sprint(show, S) == "HilbertSpace with 10 vertices"
    end

    @testset "agreement with coordinates in R^$d" for d in (2, 5, 50)
        rng = Xoshiro(d)
        X = [randn(rng, d) for _ in 1:30]
        S = HilbertSpace(X)
        scale = maximum(norm(x - y) for x in X, y in X)
        tol = 1e-12 * scale^2
        for _ in 1:300
            verts = rand(rng, 1:30, 4)
            x, y = rand_vpoint(rng, S, verts[1:2]), rand_vpoint(rng, S, verts[3:4])
            z, w = rand_vpoint(rng, S, verts[1:2]), rand_vpoint(rng, S, verts[3:4])
            rx, ry, rz, rw = realize.(Ref(X), (x, y, z, w))
            @test dist_sq(x, y) ≈ sum(abs2, rx - ry) atol = tol
            @test dist(x, y) ≈ norm(rx - ry) atol = 1e-6 * scale
            @test inner_diff(x, y, z, w) ≈ dot(rx - ry, rz - rw) atol = tol
            @test dot(x - y, z - w) == inner_diff(x, y, z, w)
            @test (x - y) isa VDiff
            # The frame is isometric: all pairwise distances agree.
            F = local_frame(x, y, z, w)
            R = (rx, ry, rz, rw)
            for i in 1:4, j in (i+1):4
                @test sum(abs2, F[i] - F[j]) ≈ sum(abs2, R[i] - R[j]) atol = tol
            end
            vs = unique((x.a, x.b, y.a, y.b))
            @test vscale(x, y) ≈ maximum(norm(X[i] - X[j]) for i in vs, j in vs) rtol =
                1e-12
        end
        # Line parameters: three input vertices (Kahan) and the general case (frame).
        for _ in 1:200
            i, j, k = randperm(rng, 30)[1:3]
            t0, h2 = line_params(vpoint(S, i), vpoint(S, j), vpoint(S, k))
            e0, e2 = line_params_xyz(X[i], X[j], X[k])
            @test t0 ≈ e0 atol = 1e-12 * scale^2 / norm(X[j] - X[k])^2
            @test h2 ≈ e2 atol = tol
            q, a, b = VPoint(S, i, j, rand(rng)), vpoint(S, k), VPoint(S, i, j, rand(rng))
            dist(a, b) > 1e-3 * scale || continue
            t0, h2 = line_params(q, a, b)
            e0, e2 = line_params_xyz(realize(X, q), realize(X, a), realize(X, b))
            @test t0 ≈ e0 atol = 1e-10
            @test h2 ≈ e2 atol = tol
        end
    end

    @testset "polynomial kernel against its feature map" begin
        # (x·y + 1)^2 on R^3 has an explicit feature map into R^10.
        k(x, y) = (dot(x, y) + 1)^2
        r2 = sqrt(2)
        phi(x) = [
            1,
            r2 * x[1],
            r2 * x[2],
            r2 * x[3],
            x[1]^2,
            x[2]^2,
            x[3]^2,
            r2 * x[1] * x[2],
            r2 * x[1] * x[3],
            r2 * x[2] * x[3],
        ]
        rng = Xoshiro(7)
        X = [randn(rng, 3) for _ in 1:15]
        F = phi.(X)
        S = HilbertSpace(X; kernel = k)
        for _ in 1:200
            verts = rand(rng, 1:15, 4)
            x, y = rand_vpoint(rng, S, verts[1:2]), rand_vpoint(rng, S, verts[3:4])
            @test dist_sq(x, y) ≈ sum(abs2, realize(F, x) - realize(F, y)) rtol = 1e-9 atol =
                1e-9
        end
        @test dist(S, 1, 2) ≈ norm(F[1] - F[2])
    end

    @testset "precision near a segment (δ/D = $r)" for r in (1e-3, 1e-6, 0.0)
        rng = Xoshiro(11)
        d = 50
        for _ in 1:50
            Qm = qr(randn(rng, d, d)).Q   # a random orientation
            a, b = Qm * [1.0; zeros(d - 1)], Qm * [-1.0; zeros(d - 1)]
            q = Qm * [0.3; r; zeros(d - 2)]
            S = HilbertSpace([a, b, q])
            t0, h2 = line_params(vpoint(S, 3), vpoint(S, 1), vpoint(S, 2))
            # D = 2; the error is about ε D^2 / δ, and about √ε D for δ = 0.
            bound = r == 0 ? 1e-7 : 20 * eps() * 4 / r
            @test abs(sqrt(h2) - r) <= bound
            @test t0 ≈ 0.35 atol = 1e-12
        end
    end

    @testset "degenerate segments" begin
        X = [[0.0, 0.0], [0.0, 0.0], [1.0, 0.0], [3.0, 4.0]]
        S = HilbertSpace(X)
        @test dist(vpoint(S, 1), vpoint(S, 2)) == 0
        @test dist(VPoint(S, 1, 2, 0.5), vpoint(S, 4)) ≈ 5
        F = local_frame(VPoint(S, 1, 2, 0.5), vpoint(S, 3), vpoint(S, 4))
        @test norm(F[2] - F[3]) ≈ norm(X[4] - X[3])
        @test norm(F[1] - F[3]) ≈ 5
        @test_throws ArgumentError line_params(vpoint(S, 4), vpoint(S, 1), vpoint(S, 2))
        @test_throws ArgumentError local_frame(
            vpoint(S, 1),
            vpoint(S, 2),
            vpoint(S, 3),
            VPoint(S, 4, 1, 0.5),
            VPoint(HilbertSpace([zeros(2)]), 1, 1, 0.0),
        )
        T = HilbertSpace([randn(3) for _ in 1:6])
        @test_throws ArgumentError local_frame(
            VPoint(T, 1, 2, 0.5),
            VPoint(T, 3, 4, 0.5),
            VPoint(T, 5, 6, 0.5),
        )
    end

    @testset "convex combinations" begin
        X = [randn(4) for _ in 1:8]
        S = HilbertSpace(X)
        x = VPoint(S, 2, 6, 0.25)
        @test convex_comb(vpoint(S, 1), vpoint(S, 5), 0.3) == VPoint(S, 1, 5, 0.3)
        @test convex_comb(vpoint(S, 1), vpoint(S, 1), 0.3) == vpoint(S, 1)
        @test convex_comb(x, vpoint(S, 6), 0.5) == VPoint(S, 2, 6, 0.625)
        @test convex_comb(vpoint(S, 6), x, 0.5) == VPoint(S, 2, 6, 0.625)
        @test convex_comb(x, VPoint(S, 6, 2, 0.25), 0.5) == VPoint(S, 2, 6, 0.5)
        y = convex_comb(x, vpoint(S, 2), 0.2)
        @test realize(X, y) ≈ 0.8 * realize(X, x) + 0.2 * X[2]
        @test_throws ArgumentError convex_comb(x, VPoint(S, 3, 4, 0.5), 0.5)
        @test_throws ArgumentError convex_comb(x, vpoint(S, 3), 0.5)
    end

    @testset "curves" begin
        rng = Xoshiro(3)
        X = [randn(rng, 6) for _ in 1:12]
        S = HilbertSpace(X)
        P = vcurve(S, [1, 4, 4, 9, 2])
        @test P isa GPntSeq{AnyDim,Float64,VPoint{typeof(S)}}
        @test length(P) == 5 && is_vertex_curve(P)
        @test all(is_generator, P.pnts)
        @test !is_generator(VPoint(S, 1, 4, 0.5))
        @test is_generator(VPoint(S, 1, 4, 1.0)) && VPoint(S, 1, 4, 1.0) == vpoint(S, 4)
        L = sum(norm(X[j] - X[i]) for (i, j) in ((1, 4), (4, 4), (4, 9), (9, 2)))
        @test geom_length(P) ≈ L
        @test BasicCompGeometry.prefix_lengths(P)[end] ≈ L
        Q = split_edges(P)
        @test length(Q) == 7 && !is_vertex_curve(Q)
        @test Q[2] == VPoint(S, 1, 4, 0.5)
        @test geom_length(Q) ≈ L
        R = GPntSeq{AnyDim,Float64,VPoint{typeof(S)}}()
        push!(R, vpoint(S, 3))
        @test BasicCompGeometry.push_distinct!(R, vpoint(S, 3)) == false
        @test BasicCompGeometry.push_distinct!(R, VPoint(S, 3, 5, 0.5)) == true
        @test length(R) == 2
        @test_throws ArgumentError dim(P)
        @test_throws ArgumentError GPntSeq{2,Float64,VPoint{typeof(S)}}([vpoint(S, 1)])
        @test_throws ArgumentError GPntSeq{AnyDim,Float32,VPoint{typeof(S)}}([vpoint(S, 1)])
        # Re-basing: points of different pairs become input vertices of a new space.
        x, y = VPoint(S, 1, 2, 0.3), VPoint(S, 5, 7, 0.6)
        S2 = HilbertSpace([x, y])
        @test dist(S2, 1, 2) ≈ dist(x, y)
        @test dist(convex_comb(vpoint(S2, 1), vpoint(S2, 2), 0.5), vpoint(S2, 1)) ≈
              dist(x, y) / 2
    end

    @testset "cache" begin
        rng = Xoshiro(5)
        X = [randn(rng, 3) for _ in 1:600]
        calls = Ref(0)
        S = HilbertSpace(X, (x, y) -> (calls[] += 1; sum(abs2, x - y)))
        n = length(X)
        for _ in 1:2, i in 1:n, j in 1:n
            @test vdist2(S, i, j) == vdist2(S, j, i)
        end
        @test calls[] == n * (n - 1) ÷ 2
        st = cache_stats(S)
        @test st.oracle_calls == calls[] && st.uncached == 0
        @test st.blocks == 6 && st.bytes == 6 * 256^2 * 8   # 3 x 3 blocks, upper triangle
        @test allocs(vdist2, S, 17, 444) == 0
        x, y = VPoint(S, 1, 300, 0.5), VPoint(S, 2, 599, 0.25)
        @test allocs(dist, x, y) == 0
        @test allocs(local_frame, x, y) == 0
        empty!(S)
        @test cache_stats(S).blocks == 0
        @test vdist2(S, 1, 2) ≈ sum(abs2, X[1] - X[2])
        @test calls[] == n * (n - 1) ÷ 2 + 1
        # A budget of no blocks: values are computed every time and counted.
        B = HilbertSpace(X; max_bytes = 0)
        @test vdist2(B, 1, 2) == vdist2(B, 1, 2) ≈ sum(abs2, X[1] - X[2])
        @test cache_stats(B) == (oracle_calls = 2, uncached = 2, blocks = 0, bytes = 0)
        @test_throws ArgumentError HilbertSpace(X; max_bytes = -1)
        bad = HilbertSpace(1:3, (i, j) -> -1.0)
        @test_throws ArgumentError vdist2(bad, 1, 2)
        @test_throws ArgumentError HilbertSpace([:a, :b])
    end

    @testset "threads share a space" begin
        rng = Xoshiro(9)
        X = [randn(rng, 4) for _ in 1:700]
        calls = Threads.Atomic{Int}(0)
        S = HilbertSpace(X, (x, y) -> (Threads.atomic_add!(calls, 1); sum(abs2, x - y)))
        n = length(X)
        pairs = [(i, j) for i in 1:n for j in (i+1):n]
        out = zeros(length(pairs))
        tasks = map(1:8) do t
            Threads.@spawn for k in shuffle(Xoshiro(t), eachindex(pairs))
                out[k] = vdist2(S, pairs[k]...)
            end
        end
        foreach(wait, tasks)
        @test all(
            k -> out[k] == sum(abs2, X[pairs[k][1]] - X[pairs[k][2]]),
            eachindex(pairs),
        )
        @test calls[] >= length(pairs)
        @test calls[] == cache_stats(S).oracle_calls
        @test cache_stats(S).blocks == 6
    end

    @testset "finite metric space" begin
        rng = Xoshiro(13)
        X = [point(randn(rng), randn(rng)) for _ in 1:40]
        S = HilbertSpace(X)
        @test size(S) == 40
        @test dist(S, 3, 3) == 0
        @test dist(S, 3, 8) ≈ dist(X[3], X[8])
        @test greedy_permutation_naive(S, 40) ==
              greedy_permutation_naive(PointsSpace(X), 40)
    end
end
