using Test
using BasicCompGeometry

@testset "Monotone Polygon Triangulation" begin
    exact = ExactKernel()

    @testset "Size preconditions and k = 3" begin
        # k < 3 throws ArgumentError
        @test_throws ArgumentError triangulate_monotone([point(0, 0), point(1, 0)])
        @test_throws ArgumentError triangulate_monotone(Point2F[])

        # k = 3 returns 0 diagonals and 1 triangle
        tri_pts = [
            point(big(0) // 1, big(0) // 1),
            point(big(2) // 1, -big(1) // 1),
            point(big(1) // 1, big(1) // 1),
        ]
        res3 = triangulate_monotone(tri_pts; kernel = exact)
        @test isempty(res3.diagonals)
        @test length(res3.triangles) == 1
        @test res3.triangles[1] == (1, 2, 3)

        # Destructuring support
        d, t = res3
        @test isempty(d)
        @test length(t) == 1
    end

    @testset "Monotone chains and merge" begin
        p = [
            point(big(0) // 1, big(0) // 1),
            point(big(1) // 1, -big(1) // 1),
            point(big(3) // 1, -big(1) // 1),
            point(big(4) // 1, big(0) // 1),
            point(big(2) // 1, big(2) // 1),
        ]
        upper_idx, lower_idx = monotone_chain_indices(p; kernel = exact)
        @test upper_idx == [1, 5, 4]
        @test lower_idx == [1, 2, 3, 4]

        upper, lower = monotone_chains(p; kernel = exact)
        @test upper == [p[1], p[5], p[4]]
        @test lower == [p[1], p[2], p[3], p[4]]

        # Both chains are strictly increasing in x
        @test all(upper[i].x < upper[i + 1].x for i in 1:length(upper)-1)
        @test all(lower[i].x < lower[i + 1].x for i in 1:length(lower)-1)

        # Merged chain has length k
        merged = merge_chains(upper, lower)
        @test length(merged) == length(p)
        @test all(merged[i].x < merged[i + 1].x for i in 1:length(merged)-1)
    end

    @testset "Triangulation correctness across shapes" begin
        test_polygons = [
            # 1. Convex quad
            [
                point(big(0) // 1, big(0) // 1),
                point(big(1) // 1, -big(1) // 1),
                point(big(2) // 1, big(0) // 1),
                point(big(1) // 1, big(1) // 1),
            ],
            # 2. Reflex vertex on upper chain
            [
                point(big(0) // 1, big(0) // 1),
                point(big(3) // 1, -big(1) // 1),
                point(big(5) // 1, big(0) // 1),
                point(big(4) // 1, big(2) // 1),
                point(big(2) // 1, big(1) // 1),
                point(big(1) // 1, big(3) // 1),
            ],
            # 3. Reflex vertex on lower chain
            [
                point(big(0) // 1, big(0) // 1),
                point(big(1) // 1, -big(3) // 1),
                point(big(2) // 1, -big(1) // 1),
                point(big(4) // 1, -big(2) // 1),
                point(big(5) // 1, big(0) // 1),
                point(big(3) // 1, big(1) // 1),
            ],
            # 4. Comb / alternating chain shape
            [
                point(big(0) // 1, big(0) // 1),
                point(big(1) // 1, -big(1) // 1),
                point(big(2) // 1, big(0) // 1),
                point(big(3) // 1, -big(1) // 1),
                point(big(4) // 1, big(0) // 1),
                point(big(2) // 1, big(2) // 1),
            ],
            # 5. Long monotone polygon with 8 vertices
            [
                point(big(0) // 1, big(0) // 1),
                point(big(1) // 1, -big(2) // 1),
                point(big(3) // 1, -big(3) // 1),
                point(big(5) // 1, -big(1) // 1),
                point(big(7) // 1, big(0) // 1),
                point(big(6) // 1, big(2) // 1),
                point(big(4) // 1, big(1) // 1),
                point(big(2) // 1, big(3) // 1),
            ],
        ]

        for (idx, poly) in enumerate(test_polygons)
            k = length(poly)
            res = triangulate_monotone(poly; kernel = exact)

            # 1. Counts: k - 3 diagonals, k - 2 triangles
            @test length(res.diagonals) == k - 3
            @test length(res.triangles) == k - 2

            # 2. Distinct diagonals and sorted endpoint order (i < j)
            @test length(unique(res.diagonals)) == k - 3
            for (u, v) in res.diagonals
                @test u < v
                @test 1 <= u <= k
                @test 1 <= v <= k
                # None is a boundary edge
                d = abs(u - v)
                @test d > 1 && d < k - 1
            end

            # 3. stack_diagonals helper matches
            @test stack_diagonals(poly; kernel = exact) == res.diagonals

            # 4. No proper crossings among diagonals
            diag_segs = [Segment(poly[u], poly[v]) for (u, v) in res.diagonals]
            for i in 1:length(diag_segs)-1
                for j in i+1:length(diag_segs)
                    @test !proper_crossing(diag_segs[i], diag_segs[j]; kernel = exact)
                end
            end

            # 5. Triangles are all counterclockwise
            for (i1, i2, i3) in res.triangles
                @test orient_sign(exact, poly[i1], poly[i2], poly[i3]) > 0
            end

            # 6. Exact area coverage: sum of triangle areas equals polygon area
            poly_area = twice_area(poly; kernel = exact)
            tri_area_sum = sum(
                twice_area([poly[t[1]], poly[t[2]], poly[t[3]]]; kernel = exact) for
                t in res.triangles
            )
            @test poly_area == tri_area_sum

            # 7. Triangle edges: internal triangle edges must match the returned diagonals
            edge_counts = Dict{Tuple{Int,Int},Int}()
            for (i1, i2, i3) in res.triangles
                for (a, b) in ((i1, i2), (i2, i3), (i3, i1))
                    edge = a < b ? (a, b) : (b, a)
                    edge_counts[edge] = get(edge_counts, edge, 0) + 1
                end
            end
            # Boundary edges appear in exactly 1 triangle; internal diagonals appear in exactly 2 triangles
            for i in 1:k
                b_edge = i < mod1(i + 1, k) ? (i, mod1(i + 1, k)) : (mod1(i + 1, k), i)
                @test get(edge_counts, b_edge, 0) == 1
            end
            for d in res.diagonals
                @test get(edge_counts, d, 0) == 2
            end
        end
    end
end
