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

    @testset "Monotonicity predicates (is_x_monotone, is_y_monotone)" begin
        p_x = [
            point(big(0) // 1, big(0) // 1),
            point(big(3) // 1, -big(1) // 1),
            point(big(5) // 1, big(0) // 1),
            point(big(4) // 1, big(2) // 1),
            point(big(2) // 1, big(1) // 1),
            point(big(1) // 1, big(3) // 1),
        ]
        @test is_x_monotone(p_x)
        @test !is_y_monotone(p_x)

        p_y = [point(-p.y, p.x) for p in p_x]
        @test !is_x_monotone(p_y)
        @test is_y_monotone(p_y)

        # Convex quad is both x-monotone and y-monotone
        p_convex = [
            point(big(0) // 1, big(0) // 1),
            point(big(2) // 1, -big(1) // 1),
            point(big(4) // 1, big(0) // 1),
            point(big(2) // 1, big(1) // 1),
        ]
        @test is_x_monotone(p_convex)
        @test is_y_monotone(p_convex)

        # U-shape is neither
        p_u = [
            point(big(0) // 1, big(0) // 1),
            point(big(4) // 1, big(0) // 1),
            point(big(4) // 1, big(4) // 1),
            point(big(3) // 1, big(4) // 1),
            point(big(3) // 1, big(1) // 1),
            point(big(1) // 1, big(1) // 1),
            point(big(1) // 1, big(4) // 1),
            point(big(0) // 1, big(4) // 1),
        ]
        @test !is_x_monotone(p_u)
        @test !is_y_monotone(p_u)

        @test !is_x_monotone(Point2F[])
        @test !is_y_monotone(Point2F[])
        @test !is_x_monotone([point(0, 0), point(1, 0)])
    end

    @testset "General Polygon Triangulation (triangulate_polygon)" begin
        @testset "Preconditions and error handling" begin
            @test_throws ArgumentError triangulate_polygon([point(0, 0), point(1, 0)])
            @test_throws ArgumentError triangulate_polygon(Point2F[])

            # Clockwise polygon throws ArgumentError
            cw_pts = [
                point(big(0) // 1, big(0) // 1),
                point(big(1) // 1, big(1) // 1),
                point(big(2) // 1, -big(1) // 1),
            ]
            @test_throws ArgumentError triangulate_polygon(cw_pts; kernel = exact)

            # k = 3
            ccw_pts = [
                point(big(0) // 1, big(0) // 1),
                point(big(2) // 1, -big(1) // 1),
                point(big(1) // 1, big(1) // 1),
            ]
            res3 = triangulate_polygon(ccw_pts; kernel = exact)
            @test isempty(res3.diagonals)
            @test length(res3.triangles) == 1
            @test res3.triangles[1] == (1, 2, 3)
            @test res3 isa PolygonTriangulation
            @test res3 isa MonotoneTriangulation
        end

        test_general_polygons = [
            # 1. U-shaped polygon (8 vertices, neither x- nor y-monotone)
            (
                "U-shape",
                [
                    point(big(0) // 1, big(0) // 1),
                    point(big(4) // 1, big(0) // 1),
                    point(big(4) // 1, big(4) // 1),
                    point(big(3) // 1, big(4) // 1),
                    point(big(3) // 1, big(1) // 1),
                    point(big(1) // 1, big(1) // 1),
                    point(big(1) // 1, big(4) // 1),
                    point(big(0) // 1, big(4) // 1),
                ],
            ),
            # 2. L-shaped polygon (6 vertices)
            (
                "L-shape",
                [
                    point(big(0) // 1, big(0) // 1),
                    point(big(3) // 1, big(0) // 1),
                    point(big(3) // 1, big(1) // 1),
                    point(big(1) // 1, big(1) // 1),
                    point(big(1) // 1, big(3) // 1),
                    point(big(0) // 1, big(3) // 1),
                ],
            ),
            # 3. Zigzag polygon (9 vertices)
            (
                "Zigzag",
                [
                    point(big(0) // 1, big(0) // 1),
                    point(big(6) // 1, big(0) // 1),
                    point(big(6) // 1, big(2) // 1),
                    point(big(2) // 1, big(2) // 1),
                    point(big(2) // 1, big(4) // 1),
                    point(big(6) // 1, big(4) // 1),
                    point(big(6) // 1, big(6) // 1),
                    point(big(0) // 1, big(6) // 1),
                    point(big(0) // 1, big(3) // 1),
                ],
            ),
            # 4. Comb polygon (10 vertices)
            (
                "Comb",
                [
                    point(big(0) // 1, big(0) // 1),
                    point(big(4) // 1, big(0) // 1),
                    point(big(4) // 1, big(3) // 1),
                    point(big(3) // 1, big(3) // 1),
                    point(big(3) // 1, big(1) // 1),
                    point(big(2) // 1, big(1) // 1),
                    point(big(2) // 1, big(3) // 1),
                    point(big(1) // 1, big(3) // 1),
                    point(big(1) // 1, big(1) // 1),
                    point(big(0) // 1, big(1) // 1),
                ],
            ),
            # 5. Star polygon (10 vertices)
            (
                "Star10",
                reverse([
                    point(big(0) // 1, big(10) // 1),
                    point(big(3) // 1, big(3) // 1),
                    point(big(10) // 1, big(3) // 1),
                    point(big(5) // 1, -big(2) // 1),
                    point(big(7) // 1, -big(9) // 1),
                    point(big(0) // 1, -big(5) // 1),
                    point(-big(7) // 1, -big(9) // 1),
                    point(-big(5) // 1, -big(2) // 1),
                    point(-big(10) // 1, big(3) // 1),
                    point(-big(3) // 1, big(3) // 1),
                ]),
            ),
            # 6. Y-monotone polygon (delegation check)
            (
                "Y-monotone",
                [point(-p.y, p.x) for p in [
                    point(big(0) // 1, big(0) // 1),
                    point(big(3) // 1, -big(1) // 1),
                    point(big(5) // 1, big(0) // 1),
                    point(big(4) // 1, big(2) // 1),
                    point(big(2) // 1, big(1) // 1),
                    point(big(1) // 1, big(3) // 1),
                ]],
            ),
            # 7. Polygon with collinear vertex
            (
                "CollinearEdge",
                [
                    point(big(0) // 1, big(0) // 1),
                    point(big(2) // 1, big(0) // 1),
                    point(big(4) // 1, big(0) // 1),
                    point(big(4) // 1, big(4) // 1),
                    point(big(0) // 1, big(4) // 1),
                ],
            ),
        ]

        for (name, poly) in test_general_polygons
            k = length(poly)
            res = triangulate_polygon(poly; kernel = exact)

            # 1. Counts: k - 3 diagonals, k - 2 triangles
            @test length(res.diagonals) == k - 3
            @test length(res.triangles) == k - 2

            # 2. Distinct diagonals and sorted endpoint order (i < j)
            @test length(unique(res.diagonals)) == k - 3
            for (u, v) in res.diagonals
                @test u < v
                @test 1 <= u <= k
                @test 1 <= v <= k
                d = abs(u - v)
                @test d > 1 && d < k - 1
            end

            # 3. polygon_diagonals helper matches
            @test polygon_diagonals(poly; kernel = exact) == res.diagonals

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
            for i in 1:k
                b_edge = i < mod1(i + 1, k) ? (i, mod1(i + 1, k)) : (mod1(i + 1, k), i)
                @test get(edge_counts, b_edge, 0) == 1
            end
            for d in res.diagonals
                @test get(edge_counts, d, 0) == 2
            end
        end

        # Test floating point kernel on larger star polygon (100 vertices)
        n = 100
        star100 = [
            point(cos(2 * pi * i / n) * (i % 2 == 0 ? 10.0 : 4.0),
                  sin(2 * pi * i / n) * (i % 2 == 0 ? 10.0 : 4.0))
            for i in 0:n-1
        ]
        res100 = triangulate_polygon(star100)
        @test length(res100.diagonals) == n - 3
        @test length(res100.triangles) == n - 2
        poly_area100 = twice_area(star100)
        tri_sum100 = sum(twice_area([star100[t[1]], star100[t[2]], star100[t[3]]]) for t in res100.triangles)
        @test isapprox(poly_area100, tri_sum100; rtol = 1e-10)

        # Invalid algorithm argument throws ArgumentError
        @test_throws ArgumentError triangulate_polygon(star100; algorithm = :bogus)
    end

    @testset "Seidel's Randomized Polygon Triangulation (triangulate_seidel)" begin
        @testset "Preconditions and error handling" begin
            @test_throws ArgumentError triangulate_seidel([point(0, 0), point(1, 0)])
            @test_throws ArgumentError triangulate_seidel(Point2F[])

            cw_pts = [
                point(big(0) // 1, big(0) // 1),
                point(big(1) // 1, big(1) // 1),
                point(big(2) // 1, -big(1) // 1),
            ]
            @test_throws ArgumentError triangulate_seidel(cw_pts; kernel = exact)

            ccw_pts = [
                point(big(0) // 1, big(0) // 1),
                point(big(2) // 1, -big(1) // 1),
                point(big(1) // 1, big(1) // 1),
            ]
            res3 = triangulate_seidel(ccw_pts; kernel = exact)
            @test isempty(res3.diagonals)
            @test length(res3.triangles) == 1
            @test res3.triangles[1] == (1, 2, 3)
            @test res3 isa PolygonTriangulation
        end

        test_seidel_shapes = [
            (
                "U-shape",
                [
                    point(big(0) // 1, big(0) // 1),
                    point(big(4) // 1, big(0) // 1),
                    point(big(4) // 1, big(4) // 1),
                    point(big(3) // 1, big(4) // 1),
                    point(big(3) // 1, big(1) // 1),
                    point(big(1) // 1, big(1) // 1),
                    point(big(1) // 1, big(4) // 1),
                    point(big(0) // 1, big(4) // 1),
                ],
            ),
            (
                "L-shape",
                [
                    point(big(0) // 1, big(0) // 1),
                    point(big(3) // 1, big(0) // 1),
                    point(big(3) // 1, big(1) // 1),
                    point(big(1) // 1, big(1) // 1),
                    point(big(1) // 1, big(3) // 1),
                    point(big(0) // 1, big(3) // 1),
                ],
            ),
            (
                "Zigzag",
                [
                    point(big(0) // 1, big(0) // 1),
                    point(big(6) // 1, big(0) // 1),
                    point(big(6) // 1, big(2) // 1),
                    point(big(2) // 1, big(2) // 1),
                    point(big(2) // 1, big(4) // 1),
                    point(big(6) // 1, big(4) // 1),
                    point(big(6) // 1, big(6) // 1),
                    point(big(0) // 1, big(6) // 1),
                    point(big(0) // 1, big(3) // 1),
                ],
            ),
            (
                "Comb",
                [
                    point(big(0) // 1, big(0) // 1),
                    point(big(4) // 1, big(0) // 1),
                    point(big(4) // 1, big(3) // 1),
                    point(big(3) // 1, big(3) // 1),
                    point(big(3) // 1, big(1) // 1),
                    point(big(2) // 1, big(1) // 1),
                    point(big(2) // 1, big(3) // 1),
                    point(big(1) // 1, big(3) // 1),
                    point(big(1) // 1, big(1) // 1),
                    point(big(0) // 1, big(1) // 1),
                ],
            ),
            (
                "Star10",
                reverse([
                    point(big(0) // 1, big(10) // 1),
                    point(big(3) // 1, big(3) // 1),
                    point(big(10) // 1, big(3) // 1),
                    point(big(5) // 1, -big(2) // 1),
                    point(big(7) // 1, -big(9) // 1),
                    point(big(0) // 1, -big(5) // 1),
                    point(-big(7) // 1, -big(9) // 1),
                    point(-big(5) // 1, -big(2) // 1),
                    point(-big(10) // 1, big(3) // 1),
                    point(-big(3) // 1, big(3) // 1),
                ]),
            ),
            (
                "CollinearEdge",
                [
                    point(big(0) // 1, big(0) // 1),
                    point(big(2) // 1, big(0) // 1),
                    point(big(4) // 1, big(0) // 1),
                    point(big(4) // 1, big(4) // 1),
                    point(big(0) // 1, big(4) // 1),
                ],
            ),
        ]

        for (name, poly) in test_seidel_shapes
            k = length(poly)
            # Direct triangulate_seidel with seeded RNG
            rng = Xoshiro(0x12345678)
            res = triangulate_seidel(poly; kernel = exact, rng = rng)

            @test length(res.diagonals) == k - 3
            @test length(res.triangles) == k - 2

            # All diagonals are internal and non-crossing
            @test length(unique(res.diagonals)) == k - 3
            for (u, v) in res.diagonals
                @test u < v
                @test 1 <= u <= k
                @test 1 <= v <= k
                d = abs(u - v)
                @test d > 1 && d < k - 1
            end

            diag_segs = [Segment(poly[u], poly[v]) for (u, v) in res.diagonals]
            for i in 1:length(diag_segs)-1
                for j in i+1:length(diag_segs)
                    @test !proper_crossing(diag_segs[i], diag_segs[j]; kernel = exact)
                end
            end

            # All triangles are CCW
            for (i1, i2, i3) in res.triangles
                @test orient_sign(exact, poly[i1], poly[i2], poly[i3]) > 0
            end

            # Area conservation
            poly_area = twice_area(poly; kernel = exact)
            tri_area_sum = sum(
                twice_area([poly[t[1]], poly[t[2]], poly[t[3]]]; kernel = exact) for
                t in res.triangles
            )
            @test poly_area == tri_area_sum

            # Via triangulate_polygon(algorithm=:seidel)
            res_opt = triangulate_polygon(poly; algorithm = :seidel, kernel = exact)
            @test length(res_opt.diagonals) == k - 3
            @test length(res_opt.triangles) == k - 2
        end

        # Test floating point kernel on 100-gon star
        n = 100
        star100 = [
            point(cos(2 * pi * i / n) * (i % 2 == 0 ? 10.0 : 4.0),
                  sin(2 * pi * i / n) * (i % 2 == 0 ? 10.0 : 4.0))
            for i in 0:n-1
        ]
        res_seidel100 = triangulate_seidel(star100)
        @test length(res_seidel100.diagonals) == n - 3
        @test length(res_seidel100.triangles) == n - 2
        poly_area100 = twice_area(star100)
        tri_sum100 = sum(twice_area([star100[t[1]], star100[t[2]], star100[t[3]]]) for t in res_seidel100.triangles)
        @test isapprox(poly_area100, tri_sum100; rtol = 1e-10)
    end
end
