# Property-based and reference tests on randomized (but seeded, deterministic)
# inputs. Fixed small hand-written cases miss bugs that only show up on
# degenerate or numerically delicate configurations; these tests check
# invariants and compare against brute-force or exact references instead.

using Test
using Random
using LinearAlgebra
using BasicCompGeometry
using BasicCompGeometry: Points, area, at, dist_segment_segment
using BasicCompGeometry.ConvexHull3D: convex_hull_3d, verify_convex_hull_3d
using BasicCompGeometry.VirtArray: orig_index
using BasicCompGeometry.LongestConvexSubset: compute_largest_convex_subset

@testset "Properties" begin

    # ------------------------------------------------------------------
    # Generators and small references
    # ------------------------------------------------------------------

    # Star-shaped (about the origin) simple polygon: jittered angles keep
    # every angular gap below pi, radii in [rmin, 1].
    function rand_star(rng, n; rmin = 0.25)
        off = 2π * rand(rng)
        angs = [off + 2π * (i - 1 + 0.9 * rand(rng)) / n for i in 1:n]
        return [
            (r = rmin + (1 - rmin) * rand(rng); point(r * cos(a), r * sin(a))) for a in angs
        ]
    end

    # y-monotone simple polygon (right chain x in [1, 2], left chain x in [-1, 0]).
    function rand_ymono(rng, n)
        ys = sort(rand(rng, n))
        left = [rand(rng, Bool) for _ in 1:n]
        right_chain = [point(1.0 + rand(rng), ys[i]) for i in 1:n if !left[i]]
        left_chain = [point(-rand(rng), ys[i]) for i in n:-1:1 if left[i]]
        return vcat(
            [point(0.0, ys[1] - 0.1)],
            right_chain,
            [point(0.0, ys[n] + 0.1)],
            left_chain,
        )
    end

    # x-monotone comb: flat bottom, alternating tall teeth and low valleys.
    function rand_comb(rng, k)
        top = [
            point(
                2k - j + 0.4 * (rand(rng) - 0.5),
                isodd(j) ? 0.1 + 0.5rand(rng) : 1 + rand(rng),
            ) for j in 0:2k
        ]
        return vcat([point(-0.5, 0.0), point(2k + 0.5, 0.0)], top)
    end

    # Random convex polygon (CCW) as a vector of points.
    rand_convex(rng, n) =
        collect(Points(convex_hull([point(randn(rng), randn(rng)) for _ in 1:n])))

    tri_area2(P, t) = twice_area([P[t[1]], P[t[2]], P[t[3]]])

    cross2(o, a, b) = (a[1] - o[1]) * (b[2] - o[2]) - (a[2] - o[2]) * (b[1] - o[1])

    # Point in a CCW convex polygon, with absolute tolerance on the turn test.
    function in_convex(p, C; tol = 1e-9)
        m = length(C)
        for i in 1:m
            cross2(C[i], C[mod1(i + 1, m)], p) < -tol && return false
        end
        return true
    end

    # Every point of A is within tol of some point of B.
    near_subset(A, B; tol = 1e-9) = all(any(norm(a - b) <= tol for b in B) for a in A)

    # Reference Minkowski sum of two convex vertex sets: hull of pairwise sums.
    ref_minkowski(A, B) = collect(Points(convex_hull([a + b for a in A for b in B])))

    # Strictly convex CCW polygon (no collinear or reflex turns).
    strictly_convex_ccw(C) =
        length(C) < 3 || all(
            cross2(C[i], C[mod1(i + 1, length(C))], C[mod1(i + 2, length(C))]) > 0 for
            i in 1:length(C)
        )

    tri_algs =
        (:sweep, :seidel, :ear_cut, :ear_cut_z_order, :ear_cut_kd_tree, :ear_cut_u_kd_tree)

    # Full combinatorial + geometric check of a triangulation of the CCW polygon P.
    function triangulation_ok(P, t; strict = true)
        n = length(P)
        length(t.triangles) == n - 2 || return false
        length(t.diagonals) == n - 3 || return false
        A = twice_area(P)
        s = zero(A)
        directed = Dict{Tuple{Int,Int},Int}()
        for tr in t.triangles
            all(1 <= i <= n for i in tr) || return false
            allunique(tr) || return false
            a = tri_area2(P, tr)
            (strict ? a > 0 : a >= 0) || return false
            s += a
            for (u, v) in ((tr[1], tr[2]), (tr[2], tr[3]), (tr[3], tr[1]))
                directed[(u, v)] = get(directed, (u, v), 0) + 1
            end
        end
        # Areas tile the polygon.
        if A isa AbstractFloat
            isapprox(s, A; rtol = 1e-9) || return false
        else
            s == A || return false
        end
        # Each polygon edge is used once (in CCW direction); every other edge
        # is a diagonal used once in each direction.
        for (e, c) in directed
            c == 1 || return false
            u, v = e
            if v == mod1(u + 1, n)
                continue
            end
            haskey(directed, (v, u)) || return false
            (min(u, v), max(u, v)) in t.diagonals || return false
        end
        for i in 1:n
            haskey(directed, (i, mod1(i + 1, n))) || return false
        end
        return true
    end

    # ------------------------------------------------------------------
    @testset "Polygon triangulation" begin
        for k in 1:60
            rng = Xoshiro(k)
            polys = (
                rand_star(rng, rand(rng, 3:40)),
                rand_ymono(rng, rand(rng, 2:25)),
                rand_comb(rng, rand(rng, 1:12)),
            )
            for P in polys
                for alg in tri_algs
                    t = triangulate_polygon(P; algorithm = alg, rng = Xoshiro(k))
                    @test triangulation_ok(P, t)
                end
                @test triangulation_ok(P, triangulate_ear_cut_u_kd_tree(P))
            end
            # x-monotone inputs for triangulate_monotone: the comb, and the
            # y-monotone polygon rotated by -90 degrees.
            C = polys[3]
            @test is_x_monotone(C)
            @test triangulation_ok(C, triangulate_monotone(C))
            Y = [point(p.y, -p.x) for p in polys[2]]
            @test is_x_monotone(Y)
            @test triangulation_ok(Y, triangulate_monotone(Y))
        end

        # Exact rational input with ExactKernel: integer star polygons.
        exact = ExactKernel()
        nexact = 0
        for k in 1:40
            rng = Xoshiro(1000 + k)
            P0 = rand_star(rng, rand(rng, 4:16))
            Pi = unique([point(round(Int, 50p.x), round(Int, 50p.y)) for p in P0])
            n = length(Pi)
            n >= 3 || continue
            # Star-shaped about the origin with strictly increasing angles.
            all(cross2(point(0, 0), Pi[i], Pi[mod1(i + 1, n)]) > 0 for i in 1:n) || continue
            P = [point(big(p.x) // 1, big(p.y) // 1) for p in Pi]
            for alg in tri_algs
                t = triangulate_polygon(
                    P;
                    kernel = exact,
                    algorithm = alg,
                    rng = Xoshiro(k),
                )
                @test triangulation_ok(P, t; strict = false)
            end
            nexact += 1
        end
        @test nexact >= 20
    end

    # ------------------------------------------------------------------
    @testset "Minkowski sums" begin
        @testset "convex vs hull of pairwise sums" begin
            for k in 1:150
                rng = Xoshiro(k)
                A = rand_convex(rng, rand(rng, 3:15))
                B = rand_convex(rng, rand(rng, 3:15))
                (length(A) >= 3 && length(B) >= 3) || continue
                M = minkowski_sum(A, B)
                R = ref_minkowski(A, B)
                @test twice_area(M) ≈ twice_area(R) rtol = 1e-10
                @test is_convex(M)
                @test twice_area(M) > 0
                @test near_subset(R, M)
                @test near_subset(M, R)
                # Orientation and starting vertex of the inputs do not matter.
                M2 = minkowski_sum(reverse(B), circshift(A, rand(rng, 0:length(A))))
                @test twice_area(M2) ≈ twice_area(R) rtol = 1e-10
            end
        end

        # Each piece must be convex and equal to (triangle of P) + Q; the
        # pieces must cover P + Q (sampled).
        function check_decomposition(P, Q, alg, seed; cover = true)
            Random.seed!(seed)
            pieces = minkowski_sum_decomposition(P, Q; algorithm = alg)
            Random.seed!(seed)
            t = triangulate_polygon(P; algorithm = alg)
            ok = length(pieces) == length(t.triangles)
            Qc = ref_minkowski(Q, [point(0.0, 0.0)])
            for (piece, tr) in zip(pieces, t.triangles)
                ok &= is_convex(piece)
                ok &= twice_area(piece) > 0
                R = ref_minkowski([P[tr[1]], P[tr[2]], P[tr[3]]], Qc)
                ok &= isapprox(twice_area(piece), twice_area(R); rtol = 1e-9)
                ok &= near_subset(R, piece) && near_subset(piece, R)
            end
            if cover && ok
                rng = Xoshiro(seed)
                for _ in 1:30
                    tr = t.triangles[rand(rng, 1:length(t.triangles))]
                    w = rand(rng, 3)
                    w ./= sum(w)
                    p = w[1] * P[tr[1]] + w[2] * P[tr[2]] + w[3] * P[tr[3]]
                    q = Q[rand(rng, 1:length(Q))]
                    ok &= any(in_convex(p + q, pc) for pc in pieces)
                end
            end
            return ok
        end

        @testset "decomposition of simple P with convex Q" begin
            for k in 1:25
                rng = Xoshiro(500 + k)
                P = rand_star(rng, rand(rng, 4:14))
                Q = rand_convex(rng, rand(rng, 3:8))
                is_convex(P) && continue
                for alg in tri_algs
                    @test check_decomposition(P, Q, alg, k)
                end
            end
        end

        @testset "decomposition: triangulators emitting clockwise triples" begin
            # Degenerate (spike) quadrilateral on which :ear_cut emits a
            # clockwise index triple; every piece must still be convex and
            # equal to the Minkowski sum of its triangle with Q.
            P = [
                point(0.74, 0.21),
                point(-0.35, 0.56),
                point(-0.03, -0.56),
                point(-0.07, -0.42),
            ]
            sq = [point(0.0, 0.0), point(1.0, 0.0), point(1.0, 1.0), point(0.0, 1.0)]
            for alg in tri_algs
                pieces = minkowski_sum_decomposition(P, sq; algorithm = alg)
                @test all(is_convex, pieces)
                @test all(pc -> twice_area(pc) > 0, pieces)
            end
            @test check_decomposition(P, sq, :ear_cut, 1; cover = false)
            # Both operands non-convex: |T(P)| * |T(Q)| convex pieces.
            for k in 1:10
                rng = Xoshiro(700 + k)
                P = rand_star(rng, rand(rng, 4:8))
                Q = rand_star(rng, rand(rng, 4:8))
                (is_convex(P) || is_convex(Q)) && continue
                pieces = minkowski_sum_decomposition(P, Q; algorithm = :ear_cut)
                @test length(pieces) == (length(P) - 2) * (length(Q) - 2)
                @test all(is_convex, pieces)
                @test all(pc -> twice_area(pc) > 0, pieces)
            end
        end
    end

    # ------------------------------------------------------------------
    @testset "Segment-segment distance vs exact reference" begin
        Q = Rational{BigInt}
        ex(p) = (Q(p[1]), Q(p[2]))
        orient_ex(a, b, c) =
            sign((b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]))
        on_box(a, b, c) =
            min(a[1], b[1]) <= c[1] <= max(a[1], b[1]) &&
            min(a[2], b[2]) <= c[2] <= max(a[2], b[2])
        function seg_intersect_ex(a, b, c, d)
            o1, o2, o3, o4 = orient_ex(a, b, c),
            orient_ex(a, b, d),
            orient_ex(c, d, a),
            orient_ex(c, d, b)
            (o1 * o2 < 0 && o3 * o4 < 0) && return true
            (o1 == 0 && on_box(a, b, c)) && return true
            (o2 == 0 && on_box(a, b, d)) && return true
            (o3 == 0 && on_box(c, d, a)) && return true
            (o4 == 0 && on_box(c, d, b)) && return true
            return false
        end
        function pt_seg_d2_ex(p, a, b)
            dx, dy = b[1] - a[1], b[2] - a[2]
            L2 = dx^2 + dy^2
            t =
                L2 == 0 ? zero(Q) :
                clamp(((p[1] - a[1]) * dx + (p[2] - a[2]) * dy) / L2, 0, 1)
            return (a[1] + t * dx - p[1])^2 + (a[2] + t * dy - p[2])^2
        end
        function ref_dist(ap, aq, bp, bq)
            a, b, c, d = ex(ap), ex(aq), ex(bp), ex(bq)
            seg_intersect_ex(a, b, c, d) && return 0.0
            d2 = min(
                pt_seg_d2_ex(a, c, d),
                pt_seg_d2_ex(b, c, d),
                pt_seg_d2_ex(c, a, b),
                pt_seg_d2_ex(d, a, b),
            )
            return Float64(sqrt(BigFloat(d2)))
        end
        # |computed - exact| <= S * (1e-10 + c * eps / sin(angle)), S the
        # coordinate scale. The distance is 1-Lipschitz in the endpoints, so a
        # backward-stable method needs only the first term; the second term
        # (c > 0) admits the O(eps / angle) error of solving the normal
        # equations for nearly parallel segments.
        # `same_dir = true` only tries the endpoint orders that keep the two
        # segments pointing the same way (no single segment reversed).
        function check(ap, aq, bp, bq; c = 0.0, same_dir = false)
            S = max(
                maximum(abs, ap),
                maximum(abs, aq),
                maximum(abs, bp),
                maximum(abs, bq),
                1e-300,
            )
            u, v = aq - ap, bq - bp
            sinθ = abs(u[1] * v[2] - u[2] * v[1]) / max(norm(u) * norm(v), 1e-300)
            r = ref_dist(ap, aq, bp, bq)
            tol = S * (1e-10 + (c == 0 ? 0.0 : c * eps() / max(sinθ, 1e-300)))
            ok = true
            orders = if same_dir
                ((ap, aq, bp, bq), (aq, ap, bq, bp), (bp, bq, ap, aq), (bq, bp, aq, ap))
            else
                ((ap, aq, bp, bq), (aq, ap, bp, bq), (bp, bq, ap, aq), (bq, bp, aq, ap))
            end
            for (x1, x2, y1, y2) in orders
                d = dist_segment_segment(x1, x2, y1, y2)
                ok &= abs(d - r) <= tol
            end
            return ok
        end
        rot(p, φ) = point(cos(φ) * p[1] - sin(φ) * p[2], sin(φ) * p[1] + cos(φ) * p[2])

        rng = Xoshiro(42)
        # Generic random segments at many scales.
        for _ in 1:300
            s = 10.0^rand(rng, -3:5)
            pts = [s * point(2rand(rng) - 1, 2rand(rng) - 1) for _ in 1:4]
            @test check(pts...)
        end
        crossing_small_ok = true
        for L in (1.0, 1e2, 1e4, 1e6), θ in (1e-3, 1e-5, 1e-6, 1e-7, 1e-9)
            for _ in 1:4
                φ = 2π * rand(rng)
                ctr = L * point(rand(rng), rand(rng))
                dir = rot(point(L, 0.0), φ)
                # Nearly parallel, crossing near the middle of both segments.
                a1, a2 = ctr - dir, ctr + dir
                b1, b2 = ctr - rot(dir, θ), ctr + rot(dir, θ)
                θ >= 1e-6 && @test check(a1, a2, b1, b2; c = 8.0, same_dir = true)
                crossing_small_ok &= check(a1, a2, b1, b2; c = 8.0)
                # Nearly parallel, disjoint: b is a shifted sideways by `off`
                # and tilted away, so the distance is attained at a1, b1.
                off = L * 10.0^rand(rng, -9:-2)
                nrm = rot(dir, π / 2) / L
                b1 = a1 + off * nrm
                @test check(a1, a2, b1, b1 + rot(2dir, θ))
            end
            # Reported case: symmetric about the origin, in this orientation.
            θ >= 1e-7 && @test dist_segment_segment(
                point(-L, 0.0),
                point(L, 0.0),
                point(-L * cos(θ), -L * sin(θ)),
                point(L * cos(θ), L * sin(θ)),
            ) <= L * (1e-10 + 8 * eps() / θ)
        end
        # Solving the normal equations (Cramer's rule) failed here: errors of
        # order eps / angle^2, dependent on endpoint order, and the endpoint
        # distance for angles below ~1e-8. The 2D path now tests for a crossing.
        @test crossing_small_ok
        offcenter_ok = true
        for θ in (1e-4, 1e-5), _ in 1:25
            φ = 2π * rand(rng)
            dir = rot(point(2.0, 0.0), φ)
            a1 = point(randn(rng), randn(rng))
            x = a1 + rand(rng, (1e-3, 1e-2, 0.3, 0.7)) * dir
            b1 = x - rand(rng) * rot(dir, θ)
            offcenter_ok &= check(a1, a1 + dir, b1, b1 + rot(dir, θ); c = 16.0)
        end
        @test offcenter_ok
        # Zero-length segments, shared endpoints, collinear (overlapping or not).
        for _ in 1:50
            p, q, r = (point(rand(rng, -5:5) / 4, rand(rng, -5:5) / 4) for _ in 1:3)
            @test check(p, p, q, r)
            @test check(p, p, q, q)
            @test check(p, q, q, r)
            u = point(rand(rng, -3:3) / 2, rand(rng, -3:3) / 2)
            t1, t2, t3, t4 = rand(rng, -4:4, 4) ./ 2
            @test check(p + t1 * u, p + t2 * u, p + t3 * u, p + t4 * u)
        end
    end

    # ------------------------------------------------------------------
    @testset "Vertical decomposition invariance" begin
        view = BBox(point(0.0, 0.0), point(1.0, 1.0))
        function rand_lines(rng, ngen, nhor)
            ls = Line{2,Float64}[]
            for _ in 1:ngen
                θ = π * (rand(rng) - 0.5) * 0.95
                push!(ls, Line(point(rand(rng), rand(rng)), point(cos(θ), sin(θ))))
            end
            for _ in 1:nhor
                push!(ls, Line(point(rand(rng), rand(rng)), point(1.0, 0.0)))
            end
            return ls
        end
        # A horizontal line re-parametrised: different anchor and direction,
        # identical y everywhere.
        redo_horizontal(rng, l) =
            Line(point(rand(rng), l.p.y), point(-2.0 * rand(rng, 1:4), 0.0))
        sum_area(d) = sum(area(t) for t in trapezoids(d))

        for k in 1:12
            rng = Xoshiro(900 + k)
            ls = rand_lines(rng, rand(rng, 15:30), rand(rng, 3:8))
            d0 = vertical_decomposition(ls, view)
            n0 = length(d0)
            @test sum_area(d0) ≈ 1.0 rtol = 1e-9
            for t in trapezoids(d0)
                @test t.x_left <= t.x_right
            end
            # Permutation invariance.
            for _ in 1:3
                dp = vertical_decomposition(shuffle(rng, ls), view)
                @test length(dp) == n0
                @test sum_area(dp) ≈ 1.0 rtol = 1e-9
            end
            # Duplicating (re-parametrised) horizontal lines changes nothing.
            hs = [l for l in ls if l.u.y == 0]
            dup = vcat(ls, [redo_horizontal(rng, l) for l in hs])
            for perm in (dup, shuffle(rng, dup))
                dd = vertical_decomposition(perm, view)
                @test length(dd) == n0
                @test sum_area(dd) ≈ 1.0 rtol = 1e-9
            end
            # Exact duplicates of arbitrary lines change nothing either.
            dd = vertical_decomposition(shuffle(rng, vcat(ls, ls[1:3])), view)
            @test length(dd) == n0
        end
    end

    # ------------------------------------------------------------------
    @testset "Visibility polygon and polygon kernel" begin
        for k in 1:40
            rng = Xoshiro(1300 + k)
            P = rand_star(rng, rand(rng, 3:30))
            A = twice_area(P)
            Ker = polygon_kernel(P)
            @test length(Ker) >= 3
            @test in_convex(point(0.0, 0.0), Ker)
            @test is_convex(Ker)
            @test twice_area(Ker) > 0
            @test twice_area(Ker) <= A * (1 + 1e-12)
            # Kernel = intersection of the inner edge halfplanes (sampled).
            n = length(P)
            for _ in 1:40
                x = point(2rand(rng) - 1, 2rand(rng) - 1)
                margins = [cross2(P[i], P[mod1(i + 1, n)], x) for i in 1:n]
                if minimum(margins) > 1e-9
                    @test in_convex(x, Ker; tol = 0.0)
                elseif minimum(margins) < -1e-9
                    @test !in_convex(x, Ker; tol = -1e-12)
                end
            end
            # From any kernel point the whole polygon is visible.
            for q in (point(0.0, 0.0), sum(Ker) / length(Ker))
                V = visibility_polygon(P, q)
                @test twice_area(V) ≈ A rtol = 1e-9
            end
            # From polygon vertices lying in the kernel (boundary vantage
            # points), and from points within rounding distance of them.
            for i in 1:n
                q = P[i]
                all(cross2(P[j], P[mod1(j + 1, n)], q) >= 0 for j in 1:n) || continue
                @test twice_area(visibility_polygon(P, q)) ≈ A rtol = 1e-9
                for s in (1e-16, 1e-13)
                    qs = q + s * point(randn(rng), randn(rng))
                    @test twice_area(visibility_polygon(P, qs)) ≈ A rtol = 1e-9
                end
            end
            # From the kernel vertices, which often sit within an ulp of a
            # polygon vertex or edge.
            for q in Ker
                @test twice_area(visibility_polygon(P, q)) ≈ A rtol = 1e-9
            end
        end
        # Non-star-shaped combs: visibility polygon is inside P, star-shaped
        # about q, and no larger than P.
        for k in 1:20
            rng = Xoshiro(1400 + k)
            C = rand_comb(rng, rand(rng, 2:8))
            q = point(2k % 7 + 0.3 * rand(rng), 0.05)
            q = point(clamp(q.x, 0.0, length(C) / 2 - 2.0), q.y)
            V = visibility_polygon(C, q)
            @test 0 < twice_area(V) <= twice_area(C) * (1 + 1e-12)
            # Inside C, up to rounding of vertices computed on its edges.
            near_bd(v) =
                minimum(
                    BasicCompGeometry.dist_point_segment(
                        v,
                        C[i],
                        C[mod1(i + 1, length(C))],
                    ) for i in 1:length(C)
                ) <= 1e-9
            @test all(is_inside(v, C) || near_bd(v) for v in V)
            @test in_convex(q, polygon_kernel(V); tol = 1e-9)
        end
    end

    # ------------------------------------------------------------------
    @testset "Convex hulls" begin
        @testset "2D" begin
            for k in 1:100
                rng = Xoshiro(1500 + k)
                n = rand(rng, 3:60)
                pts = if isodd(k)
                    [point(randn(rng), randn(rng)) for _ in 1:n]
                else
                    # Integer grid: many duplicates and collinear triples.
                    [point(rand(rng, 0:6), rand(rng, 0:6)) for _ in 1:n]
                end
                H = collect(Points(convex_hull(pts)))
                length(unique(pts)) >= 3 || continue
                all(cross2(pts[1], pts[2], p) == 0 for p in pts) && continue
                @test strictly_convex_ccw(H)
                @test all(h in pts for h in H)
                @test allunique(H)
                tol = eltype(eltype(pts)) <: Integer ? 0 : 1e-12
                @test all(in_convex(p, H; tol = tol) for p in pts)
                lo = minimum(p -> (p[1], p[2]), pts)
                @test H[1] == point(lo...)
            end
        end
        @testset "3D" begin
            for k in 1:20
                rng = Xoshiro(1600 + k)
                n = rand(rng, 4:120)
                on_sphere = isodd(k)
                pts = [
                    (
                        v = point(randn(rng), randn(rng), randn(rng));
                        on_sphere ? v / norm(v) : v
                    ) for _ in 1:n
                ]
                Random.seed!(k)
                verts, faces = convex_hull_3d(pts)
                @test all(1 <= v <= n for v in verts)
                @test allunique(verts)
                ok, msg = verify_convex_hull_3d(pts, verts, faces)
                @test ok
                # Every input point on the inner side of every face.
                H = pts[verts]
                for f in faces
                    nrm = cross(H[f[2]] - H[f[1]], H[f[3]] - H[f[1]])
                    @test norm(nrm) > 0
                    @test all(dot(nrm, p - H[f[1]]) <= 1e-9 * norm(nrm) for p in pts)
                end
                # Euler characteristic of a triangulated sphere.
                V, F = length(verts), length(faces)
                @test 3F == 2 * (3V - 6)
                # Extreme points in random directions are hull vertices.
                for _ in 1:20
                    u = point(randn(rng), randn(rng), randn(rng))
                    @test maximum(dot(u, p) for p in H) ≈ maximum(dot(u, p) for p in pts)
                end
                on_sphere && @test V == n
            end
        end
    end

    # ------------------------------------------------------------------
    @testset "Largest convex subset vs brute force" begin
        function in_convex_position(S)
            m = length(S)
            m <= 2 && return true
            for i in 1:m, j in 1:m, l in 1:m, r in 1:m
                # S[r] inside triangle (i, j, l) => not in convex position.
                length(unique((i, j, l, r))) == 4 || continue
                o = cross2(S[i], S[j], S[l])
                o == 0 && continue
                s1 = cross2(S[i], S[j], S[r]) * o
                s2 = cross2(S[j], S[l], S[r]) * o
                s3 = cross2(S[l], S[i], S[r]) * o
                (s1 >= 0 && s2 >= 0 && s3 >= 0) && return false
            end
            return true
        end
        ntested = 0
        for k in 1:60
            rng = Xoshiro(1700 + k)
            n = rand(rng, 3:9)
            pts = unique([point(rand(rng, 0:20), rand(rng, 0:20)) for _ in 1:n])
            n = length(pts)
            # General position: no three collinear.
            any(
                cross2(pts[i], pts[j], pts[l]) == 0 for i in 1:n for j in (i+1):n for
                l in (j+1):n
            ) && continue
            best = 0
            for mask in 1:(2^n-1)
                S = [pts[i] for i in 1:n if (mask >> (i - 1)) & 1 == 1]
                length(S) > best && in_convex_position(S) && (best = length(S))
            end
            res = collect(Points(compute_largest_convex_subset(PntSeq(pts))))
            @test length(res) == best
            @test all(p in pts for p in res)
            @test allunique(res)
            @test in_convex_position(res)
            ntested += 1
        end
        @test ntested >= 30
    end

    # ------------------------------------------------------------------
    @testset "Range minimum queries vs brute force" begin
        # Reference: the greatest index among equivalent minima.
        function ref_argmin(keys, less, lo, hi)
            best = lo
            for i in (lo+1):hi
                less(keys[best], keys[i]) || (best = i)
            end
            return best
        end
        for k in 1:40
            rng = Xoshiro(1800 + k)
            n = rand(rng, [1:20; 60:70; 120:300])
            keys = isodd(k) ? rand(rng, 1:5, n) : randn(rng, n)
            for less in (isless, (a, b) -> a > b)
                idx_b = build_rmq(keys; less = less)
                idx_s = build_rmq(keys; less = less, backend = :sparse)
                queries = if n <= 70
                    [(lo, hi) for lo in 1:n for hi in lo:n]
                else
                    [(sort(rand(rng, 1:n, 2))...,) for _ in 1:2000]
                end
                ok_b = ok_s = true
                for (lo, hi) in queries
                    r = ref_argmin(keys, less, lo, hi)
                    ok_b &= range_argmin(idx_b, lo, hi) == r
                    ok_s &= range_argmin(idx_s, lo, hi) == r
                end
                @test ok_b
                @test ok_s
            end
        end
    end

    # ------------------------------------------------------------------
    @testset "WSPD covers every pair exactly once" begin
        for k in 1:12
            rng = Xoshiro(1900 + k)
            n = rand(rng, 2:40)
            D = isodd(k) ? 2 : 3
            pts = [point((randn(rng) for _ in 1:D)...) for _ in 1:n]
            if k % 4 == 0  # clustered input
                pts = [p * (rand(rng) < 0.5 ? 1e-3 : 1.0) for p in pts]
            end
            for sep in (0.25, 1.0, 2.0)
                W = WSPD.init(PntSeq(pts), sep)
                finals = WSPD.expand!(W)
                cover = Dict{Tuple{Int,Int},Int}()
                sep_ok = true
                for pair in finals
                    L = [orig_index(W.tree.PS, i) for i in pair.left.r]
                    R = [orig_index(W.tree.PS, i) for i in pair.right.r]
                    sep_ok &= isempty(intersect(L, R))
                    # Geometric separation, recomputed from the points.
                    dmin = minimum(dist(pts[i], pts[j]) for i in L for j in R)
                    dL = maximum(dist(pts[i], pts[j]) for i in L for j in L)
                    dR = maximum(dist(pts[i], pts[j]) for i in R for j in R)
                    sep_ok &= max(dL, dR) <= sep * dmin * (1 + 1e-9)
                    for i in L, j in R
                        key = (min(i, j), max(i, j))
                        cover[key] = get(cover, key, 0) + 1
                    end
                end
                @test sep_ok
                @test length(cover) == n * (n - 1) ÷ 2
                @test all(==(1), values(cover))
            end
        end
    end

    # ------------------------------------------------------------------
    @testset "Approximate diameter within factor" begin
        for k in 1:30
            rng = Xoshiro(2000 + k)
            n = rand(rng, 2:200)
            D = rand(rng, 2:3)
            pts = [point((randn(rng) for _ in 1:D)...) for _ in 1:n]
            isodd(k) && (pts = [p * (1 + 0.01 * rand(rng)) / norm(p) for p in pts])
            dref = maximum(dist(p, q) for p in pts for q in pts)
            P = PntSeq(pts)
            @test exact_diameter(P) ≈ dref
            for ε in (1.0, 0.25, 0.05, 0.01)
                d = approx_diameter(P, ε)
                @test d <= dref * (1 + 1e-12)
                @test d >= dref / (1 + ε) * (1 - 1e-12)
            end
        end
    end

    # ------------------------------------------------------------------
    @testset "Point-sequence sampling lies on the curve" begin
        for k in 1:30
            rng = Xoshiro(2100 + k)
            m = rand(rng, 2:15)
            # x-monotone polyline (x strictly increasing): arc-length position
            # of a point on it is determined by its x coordinate.
            xs = cumsum(0.05 .+ rand(rng, m))
            V = [point(xs[i], 3randn(rng)) for i in 1:m]
            isodd(k) && m >= 3 && (V[2] = V[1] + point(1e-12, 0.0))  # near-zero edge
            xs = [v.x for v in V]
            P = PntSeq(V)
            plen = vcat(0.0, cumsum([dist(V[i], V[i+1]) for i in 1:(m-1)]))
            total = plen[end]
            dist_to_edge(p, a, b) = BasicCompGeometry.dist_point_segment(p, a, b)
            function arclen_of(p)
                i = clamp(searchsortedlast(xs, p.x), 1, m - 1)
                return plen[i] + dist(V[i], p), dist_to_edge(p, V[i], V[i+1])
            end
            for ns in (2, 3, rand(rng, 4:60))
                S = collect(Points(sample_uniformly(P, ns)))
                @test length(S) == ns
                @test S[1] == V[1]
                @test S[end] == V[end]
                ok = true
                for (j, p) in enumerate(S)
                    s, off = arclen_of(p)
                    ok &= off <= 1e-9 * (1 + total)
                    ok &= abs(s - (j - 1) / (ns - 1) * total) <= 1e-9 * (1 + total)
                end
                @test ok
            end
            @test at(P, 0.0) == V[1]
            @test at(P, 1.0) == V[end]
            ok = true
            for t in rand(rng, 20)
                s, off = arclen_of(at(P, t))
                ok &= off <= 1e-9 * (1 + total)
                ok &= abs(s - t * total) <= 1e-9 * (1 + total)
            end
            @test ok
        end
    end
end
