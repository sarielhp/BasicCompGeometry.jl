using Test
using Random
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
        @test all(upper[i].x < upper[i+1].x for i in 1:(length(upper)-1))
        @test all(lower[i].x < lower[i+1].x for i in 1:(length(lower)-1))

        # Merged chain has length k
        merged = merge_chains(upper, lower)
        @test length(merged) == length(p)
        @test all(merged[i].x < merged[i+1].x for i in 1:(length(merged)-1))
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
            for i in 1:(length(diag_segs)-1)
                for j in (i+1):length(diag_segs)
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
                [
                    point(-p.y, p.x) for p in [
                        point(big(0) // 1, big(0) // 1),
                        point(big(3) // 1, -big(1) // 1),
                        point(big(5) // 1, big(0) // 1),
                        point(big(4) // 1, big(2) // 1),
                        point(big(2) // 1, big(1) // 1),
                        point(big(1) // 1, big(3) // 1),
                    ]
                ],
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
            for i in 1:(length(diag_segs)-1)
                for j in (i+1):length(diag_segs)
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
            point(
                cos(2 * pi * i / n) * (i % 2 == 0 ? 10.0 : 4.0),
                sin(2 * pi * i / n) * (i % 2 == 0 ? 10.0 : 4.0),
            ) for i in 0:(n-1)
        ]
        res100 = triangulate_polygon(star100)
        @test res100.diagonals == triangulate_sweep(star100).diagonals
        @test length(res100.diagonals) == n - 3
        @test length(res100.triangles) == n - 2
        poly_area100 = twice_area(star100)
        tri_sum100 = sum(
            twice_area([star100[t[1]], star100[t[2]], star100[t[3]]]) for
            t in res100.triangles
        )
        @test isapprox(poly_area100, tri_sum100; rtol = 1e-10)

        res100_ear = triangulate_polygon(star100; algorithm = :ear_cut)
        @test triangulate_ear_cut(star100).diagonals == res100_ear.diagonals
        @test_throws ArgumentError triangulate_ear_cut([point(0.0, 0.0), point(1.0, 0.0)])
        @test_throws ArgumentError triangulate_ear_cut([
            point(0.0, 0.0),
            point(0.0, 1.0),
            point(1.0, 0.0),
        ])
        @test triangulate_ear_cut([point(0.0, 0.0), point(1.0, 0.0), point(0.0, 1.0)]).triangles ==
              [(1, 2, 3)]
        @test length(res100_ear.diagonals) == n - 3
        @test length(res100_ear.triangles) == n - 2
        tri_sum100_ear = sum(
            twice_area([star100[t[1]], star100[t[2]], star100[t[3]]]) for
            t in res100_ear.triangles
        )
        @test isapprox(poly_area100, tri_sum100_ear; rtol = 1e-10)

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
            for i in 1:(length(diag_segs)-1)
                for j in (i+1):length(diag_segs)
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
            point(
                cos(2 * pi * i / n) * (i % 2 == 0 ? 10.0 : 4.0),
                sin(2 * pi * i / n) * (i % 2 == 0 ? 10.0 : 4.0),
            ) for i in 0:(n-1)
        ]
        res_seidel100 = triangulate_seidel(star100)
        @test length(res_seidel100.diagonals) == n - 3
        @test length(res_seidel100.triangles) == n - 2
        poly_area100 = twice_area(star100)
        tri_sum100 = sum(
            twice_area([star100[t[1]], star100[t[2]], star100[t[3]]]) for
            t in res_seidel100.triangles
        )
        @test isapprox(poly_area100, tri_sum100; rtol = 1e-10)

        @testset "Trapezoidal map invariants, spirals, collinearity, and O(n log* n) scaling" begin
            function _verify_seidel_result(poly, res; kernel = FloatingKernel())
                m = length(poly)
                @test length(res.diagonals) == m - 3
                @test length(res.triangles) == m - 2
                @test length(unique(res.diagonals)) == m - 3
                edges = [Segment(poly[i], poly[mod1(i + 1, m)]) for i in 1:m]
                diag_segs = [Segment(poly[u], poly[v]) for (u, v) in res.diagonals]
                for i in 1:length(diag_segs)
                    for j in (i+1):length(diag_segs)
                        @test !proper_crossing(diag_segs[i], diag_segs[j]; kernel = kernel)
                    end
                    for e in edges
                        @test !proper_crossing(diag_segs[i], e; kernel = kernel)
                    end
                end
                for (a, b, c) in res.triangles
                    @test orient_sign(kernel, poly[a], poly[b], poly[c]) > 0
                end
                poly_a = twice_area(poly; kernel = kernel)
                tri_a = sum(
                    twice_area([poly[a], poly[b], poly[c]]; kernel = kernel) for
                    (a, b, c) in res.triangles
                )
                if kernel isa ExactKernel
                    @test poly_a == tri_a
                else
                    @test isapprox(poly_a, tri_a; rtol = 1e-9)
                end
            end

            function _verify_trap_invariants(poly, rng; kernel = FloatingKernel())
                m = BasicCompGeometry._seidel_build_trapezoidal_map(kernel, poly, rng)
                traps = m.traps
                dag = m.dag
                alive_cnt = 0
                for (t_id, T) in enumerate(traps)
                    T.alive || continue
                    alive_cnt += 1
                    @test dag[T.sink].kind == BasicCompGeometry._SEIDEL_DAG_LEAF
                    @test dag[T.sink].key == t_id
                    @test BasicCompGeometry._seidel_cmp_pt(
                        kernel,
                        m.pts,
                        T.leftp,
                        T.rightp,
                    ) < 0
                    @test (T.l0 == 0 && T.l1 == 0) || (T.l0 != T.l1)
                    @test (T.r0 == 0 && T.r1 == 0) || (T.r0 != T.r1)
                    for nb in (T.l0, T.l1)
                        if nb != 0
                            @test traps[nb].alive
                            @test traps[nb].r0 == t_id || traps[nb].r1 == t_id
                            @test traps[nb].rightp == T.leftp
                        end
                    end
                    for nb in (T.r0, T.r1)
                        if nb != 0
                            @test traps[nb].alive
                            @test traps[nb].l0 == t_id || traps[nb].l1 == t_id
                            @test traps[nb].leftp == T.rightp
                        end
                    end
                end
                @test alive_cnt <= 3 * m.n + 1

                # Verify extracted pieces are CCW and strictly compare_xy-monotone
                diags = BasicCompGeometry._seidel_extract_diagonals(m)
                pieces = BasicCompGeometry._seidel_extract_pieces(length(poly), diags)
                @test sum(length(pc) - 2 for pc in pieces) == length(poly) - 2
                for pc in pieces
                    @test twice_area([poly[idx] for idx in pc]; kernel = kernel) > 0
                end
            end

            # Vertical-edge quad and multi-side collinear triangle
            vquad =
                [point(4//1, 5//1), point(0//1, 4//1), point(0//1, 0//1), point(4//1, 0//1)]
            _verify_trap_invariants(vquad, Xoshiro(7); kernel = exact)
            _verify_seidel_result(
                vquad,
                triangulate_seidel(vquad; kernel = exact, rng = Xoshiro(7));
                kernel = exact,
            )

            collinear_tri = [
                point(0//1, 0//1),
                point(1//1, 0//1),
                point(2//1, 0//1),
                point(3//1, 0//1),
                point(5//2, 1//1),
                point(2//1, 2//1),
                point(3//2, 3//1),
                point(1//1, 2//1),
                point(1//2, 1//1),
            ]
            for seed in 1:5
                _verify_trap_invariants(collinear_tri, Xoshiro(seed); kernel = exact)
                _verify_seidel_result(
                    collinear_tri,
                    triangulate_seidel(collinear_tri; kernel = exact, rng = Xoshiro(seed));
                    kernel = exact,
                )
            end

            # Spiral polygon
            outer = Point2F[]
            inner = Point2F[]
            for i in 0:32
                θ = 2 * pi * i / 8
                r_out = 2.0 + 0.5 * i
                r_in = 1.0 + 0.5 * i
                push!(outer, point(r_out * cos(θ), r_out * sin(θ)))
                push!(inner, point(r_in * cos(θ), r_in * sin(θ)))
            end
            spiral = vcat(outer, reverse(inner))
            for seed in 1:5
                _verify_trap_invariants(spiral, Xoshiro(seed))
                _verify_seidel_result(
                    spiral,
                    triangulate_seidel(spiral; rng = Xoshiro(seed)),
                )
            end

            # Orthogonal snake with many shared x and y coordinates
            snake = [point(0.0, 0.0), point(10.0, 0.0)]
            y_cur = 0.0
            for _ in 1:12
                push!(snake, point(10.0, y_cur + 1.0))
                push!(snake, point(2.0, y_cur + 1.0))
                push!(snake, point(2.0, y_cur + 2.0))
                push!(snake, point(10.0, y_cur + 2.0))
                y_cur += 2.0
            end
            push!(snake, point(10.0, y_cur + 1.0), point(0.0, y_cur + 1.0))
            _verify_trap_invariants(snake, Xoshiro(11))
            _verify_seidel_result(snake, triangulate_seidel(snake; rng = Xoshiro(11)))

            # Phase schedule sanity: n = 65536 yields 4 phases (<= log* n)
            bounds = BasicCompGeometry._seidel_phase_bounds(65536)
            @test bounds == [4096, 16384, 32768, 65536]

            # Double-comb scaling & custom kernel forwarding
            mutable struct _SeidelCountKernel <: AbstractKernel
                turns::Int
                cmps::Int
            end
            BasicCompGeometry.turn_sign(
                k::_SeidelCountKernel,
                p::Point{2},
                q::Point{2},
                r::Point{2},
            ) = (k.turns += 1; turn_sign(FloatingKernel(), p, q, r))
            BasicCompGeometry.compare_xy(k::_SeidelCountKernel, p::Point{2}, q::Point{2}) =
                (k.cmps += 1; compare_xy(FloatingKernel(), p, q))

            function _double_comb(m_teeth::Int)
                pts = Point2F[point(0.0, 0.0), point(20.0, 0.0)]
                for i in 1:m_teeth
                    push!(pts, point(20.0, 2.0 * i - 1.0), point(12.0, 2.0 * i))
                end
                push!(
                    pts,
                    point(20.0, 2.0 * m_teeth + 1.0),
                    point(0.0, 2.0 * m_teeth + 1.0),
                )
                for i in m_teeth:-1:1
                    push!(pts, point(0.0, 2.0 * i + 0.5), point(8.0, 2.0 * i - 0.5))
                end
                return pts
            end

            comb_small = _double_comb(250)  # n = 1004
            comb_large = _double_comb(2000) # n = 8004
            kc_small = _SeidelCountKernel(0, 0)
            kc_large = _SeidelCountKernel(0, 0)
            res_small = triangulate_seidel(comb_small; kernel = kc_small, rng = Xoshiro(1))
            res_large = triangulate_seidel(comb_large; kernel = kc_large, rng = Xoshiro(1))
            @test length(res_small.triangles) == length(comb_small) - 2
            @test length(res_large.triangles) == length(comb_large) - 2
            @test isapprox(
                twice_area(comb_large),
                sum(
                    twice_area([comb_large[a], comb_large[b], comb_large[c]]) for
                    (a, b, c) in res_large.triangles
                );
                rtol = 1e-8,
            )
            # Predicate work per vertex remains bounded by a small O(log* n) constant
            @test kc_small.turns / length(comb_small) < 20.0
            @test kc_large.turns / length(comb_large) < 20.0
            @test kc_large.cmps / length(comb_large) < 25.0

            # Ear clipping with reflex-vertex list stays O(n * (r + 1)) rather than O(n^2)
            function _non_monotone_few_reflex(n::Int)
                pts = Point2F[]
                for i in 0:(n-3)
                    th = 0.2 + 1.8 * pi * i / (n - 2)
                    push!(pts, point(cos(th), sin(th)))
                end
                # A zigzag notch in the gap between the arc's ends (breaks both monotonicities).
                push!(pts, point(0.5, 0.1), point(0.7, -0.1))
                return pts
            end
            for n_ear in (500, 4000)
                poly_ear = _non_monotone_few_reflex(n_ear)
                @test !is_x_monotone(poly_ear) && !is_y_monotone(poly_ear)
                ke = _SeidelCountKernel(0, 0)
                res_ear = triangulate_polygon(poly_ear; kernel = ke, algorithm = :ear_cut)
                @test length(res_ear.triangles) == n_ear - 2
                @test ke.turns / n_ear < 15.0
                @test all(
                    t -> orient_sign(poly_ear[t[1]], poly_ear[t[2]], poly_ear[t[3]]) >= 0,
                    res_ear.triangles,
                )
                @test isapprox(
                    sum(
                        twice_area([poly_ear[a], poly_ear[b], poly_ear[c]]) for
                        (a, b, c) in res_ear.triangles
                    ),
                    twice_area(poly_ear);
                    rtol = 1e-9,
                )
                @test_nowarn triangulate_sweep(poly_ear)
            end
            # r = 0 convex polygon with vertical edge also stays linear O(n)
            poly_vert_conv = Point2F[point(0.0, 1.0), point(0.0, -1.0)]
            for i in 1:1998
                th = -pi / 2 + pi * i / 1999
                push!(poly_vert_conv, point(cos(th), sin(th)))
            end
            kvc = _SeidelCountKernel(0, 0)
            res_vc = triangulate_polygon(poly_vert_conv; kernel = kvc, algorithm = :ear_cut)
            @test length(res_vc.triangles) == length(poly_vert_conv) - 2
            @test kvc.turns / length(poly_vert_conv) < 5.0
        end
    end

    @testset "Deterministic O(n log n) Plane-Sweep Triangulation (triangulate_sweep)" begin
        @testset "Preconditions and invalid inputs" begin
            @test_throws ArgumentError triangulate_sweep([point(0, 0), point(1, 0)])
            @test_throws ArgumentError triangulate_sweep(Point2F[])

            cw_pts = [
                point(big(0) // 1, big(0) // 1),
                point(big(1) // 1, big(1) // 1),
                point(big(2) // 1, -big(1) // 1),
            ]
            @test_throws ArgumentError triangulate_sweep(cw_pts; kernel = exact)

            dup_pts = [point(0.0, 0.0), point(2.0, 0.0), point(2.0, 0.0), point(0.0, 2.0)]
            @test_throws ArgumentError triangulate_sweep(dup_pts)

            spike_pts = [
                point(0.0, 0.0),
                point(4.0, 0.0),
                point(2.0, 0.0),
                point(4.0, 2.0),
                point(0.0, 2.0),
            ]
            @test_throws ArgumentError triangulate_sweep(spike_pts)

            # Non-simple self-crossing polygon with positive signed area
            bowtie = [point(0.0, 0.0), point(6.0, 0.0), point(2.0, 8.0), point(5.0, -2.0)]
            @test_throws ArgumentError triangulate_sweep(bowtie)
        end

        @testset "Sweep AVL tree structural invariants" begin
            function _verify_avl(t::BasicCompGeometry._SweepAVLTree)
                seen = Int[]
                function check(u::Int, p::Int)::Int
                    u == 0 && return 0
                    @test t.parent[u] == p
                    @test t.slot[t.key[u]] == u
                    hl = check(t.left[u], u)
                    push!(seen, t.key[u])
                    hr = check(t.right[u], u)
                    @test abs(hl - hr) <= 1
                    h = 1 + max(hl, hr)
                    @test t.height[u] == h
                    return h
                end
                check(t.root, 0)
                return seen
            end

            # Horizontal stack of 64 segments: edge i is y = i for x in [0, 100]
            m = 64
            pts = Point2F[]
            next_v = collect(1:2m)
            for i in 1:m
                push!(pts, point(Float64(i) * 1e-3, Float64(i)))
                push!(pts, point(100.0 + Float64(i) * 1e-3, Float64(i)))
                next_v[2i-1] = 2i
                next_v[2i] = 2i - 1
            end
            tree = BasicCompGeometry._SweepAVLTree(2m)
            rng = Xoshiro(99)
            perm = randperm(rng, m)
            active = Int[]
            for i in perm
                e = 2i - 1
                BasicCompGeometry._avl_insert!(tree, FloatingKernel(), pts, next_v, e, e)
                push!(active, i)
                sort!(active)
                seen = _verify_avl(tree)
                @test seen == [2k - 1 for k in active]
                @test tree.height[tree.root] <= ceil(Int, 1.44 * log2(length(active) + 2))
            end

            # Query points at various heights between edges
            for yq in (0.5, 1.5, 32.5, 64.5)
                q_idx = length(pts) + 1
                qpts = vcat(pts, [point(50.0, yq)])
                got = BasicCompGeometry._avl_below(
                    tree,
                    FloatingKernel(),
                    qpts,
                    next_v,
                    q_idx,
                )
                expected = yq < 1.0 ? 0 : 2 * floor(Int, yq) - 1
                @test got == expected
            end

            for i in randperm(rng, m)
                e = 2i - 1
                BasicCompGeometry._avl_delete!(tree, e)
                filter!(!=(i), active)
                seen = _verify_avl(tree)
                @test seen == [2k - 1 for k in active]
            end
            @test tree.root == 0
        end

        function _verify_sweep_result(poly, res; kernel = FloatingKernel())
            m = length(poly)
            @test length(res.diagonals) == m - 3
            @test length(res.triangles) == m - 2
            @test length(unique(res.diagonals)) == m - 3
            for (u, v) in res.diagonals
                @test 1 <= u < v <= m
                d = abs(u - v)
                @test d > 1 && d < m - 1
            end

            edges = [Segment(poly[i], poly[mod1(i + 1, m)]) for i in 1:m]
            diag_segs = [Segment(poly[u], poly[v]) for (u, v) in res.diagonals]
            for i in 1:length(diag_segs)
                for j in (i+1):length(diag_segs)
                    @test !proper_crossing(diag_segs[i], diag_segs[j]; kernel = kernel)
                end
                for e in edges
                    @test !proper_crossing(diag_segs[i], e; kernel = kernel)
                end
            end

            edge_counts = Dict{Tuple{Int,Int},Int}()
            for (a, b, c) in res.triangles
                @test orient_sign(kernel, poly[a], poly[b], poly[c]) > 0
                for (u, v) in ((a, b), (b, c), (c, a))
                    edge = u < v ? (u, v) : (v, u)
                    edge_counts[edge] = get(edge_counts, edge, 0) + 1
                end
            end
            for i in 1:m
                b_edge = i < mod1(i + 1, m) ? (i, mod1(i + 1, m)) : (mod1(i + 1, m), i)
                @test get(edge_counts, b_edge, 0) == 1
            end
            for d in res.diagonals
                @test get(edge_counts, d, 0) == 2
            end

            poly_a = twice_area(poly; kernel = kernel)
            tri_a = sum(
                twice_area([poly[a], poly[b], poly[c]]; kernel = kernel) for
                (a, b, c) in res.triangles
            )
            if kernel isa ExactKernel
                @test poly_a == tri_a
            else
                @test isapprox(poly_a, tri_a; rtol = 1e-9)
            end
        end

        # Vertical edges, shared-x jagged reflex/merge chains, collinear chains, spirals, snakes
        vquad = [
            point(4 // 1, 5 // 1),
            point(0 // 1, 4 // 1),
            point(0 // 1, 0 // 1),
            point(4 // 1, 0 // 1),
        ]
        _verify_sweep_result(
            vquad,
            triangulate_sweep(vquad; kernel = exact);
            kernel = exact,
        )

        vcollinear = [
            point(0 // 1, 3 // 1),
            point(0 // 1, 2 // 1),
            point(0 // 1, 1 // 1),
            point(0 // 1, 0 // 1),
            point(2 // 1, 0 // 1),
            point(2 // 1, 3 // 1),
        ]
        _verify_sweep_result(
            vcollinear,
            triangulate_sweep(vcollinear; kernel = exact);
            kernel = exact,
        )

        collinear_tri = [
            point(0 // 1, 0 // 1),
            point(1 // 1, 0 // 1),
            point(2 // 1, 0 // 1),
            point(3 // 1, 0 // 1),
            point(5 // 2, 1 // 1),
            point(2 // 1, 2 // 1),
            point(3 // 2, 3 // 1),
            point(1 // 1, 2 // 1),
            point(1 // 2, 1 // 1),
        ]
        _verify_sweep_result(
            collinear_tri,
            triangulate_sweep(collinear_tri; kernel = exact);
            kernel = exact,
        )

        shared_x_jagged = [
            point(0 // 1, 0 // 1),
            point(2 // 1, 0 // 1),
            point(1 // 1, 1 // 1),
            point(2 // 1, 2 // 1),
            point(1 // 1, 3 // 1),
            point(2 // 1, 4 // 1),
            point(0 // 1, 4 // 1),
            point(1 // 1, 2 // 1),
        ]
        _verify_sweep_result(
            shared_x_jagged,
            triangulate_sweep(shared_x_jagged; kernel = exact);
            kernel = exact,
        )

        u_poly = [
            point(0 // 1, 0 // 1),
            point(4 // 1, 0 // 1),
            point(4 // 1, 4 // 1),
            point(3 // 1, 4 // 1),
            point(3 // 1, 1 // 1),
            point(1 // 1, 1 // 1),
            point(1 // 1, 4 // 1),
            point(0 // 1, 4 // 1),
        ]
        res_u_opt = triangulate_polygon(u_poly; kernel = exact, algorithm = :sweep)
        _verify_sweep_result(u_poly, res_u_opt; kernel = exact)
        @test res_u_opt.diagonals == triangulate_sweep(u_poly; kernel = exact).diagonals

        outer = Point2F[]
        inner = Point2F[]
        for i in 0:32
            θ = 2 * pi * i / 8
            r_out = 2.0 + 0.5 * i
            r_in = 1.0 + 0.5 * i
            push!(outer, point(r_out * cos(θ), r_out * sin(θ)))
            push!(inner, point(r_in * cos(θ), r_in * sin(θ)))
        end
        spiral = vcat(outer, reverse(inner))
        _verify_sweep_result(spiral, triangulate_sweep(spiral))

        snake = [point(0.0, 0.0), point(10.0, 0.0)]
        y_cur = 0.0
        for _ in 1:12
            push!(snake, point(10.0, y_cur + 1.0))
            push!(snake, point(2.0, y_cur + 1.0))
            push!(snake, point(2.0, y_cur + 2.0))
            push!(snake, point(10.0, y_cur + 2.0))
            y_cur += 2.0
        end
        push!(snake, point(10.0, y_cur + 1.0), point(0.0, y_cur + 1.0))
        _verify_sweep_result(snake, triangulate_sweep(snake))

        # Horizontal comb with Theta(n) simultaneously active edges in the sweep AVL tree:
        # verifies O(n log n) predicate scaling rather than O(n^2)
        mutable struct _SweepCountKernel <: AbstractKernel
            turns::Int
            cmps::Int
        end
        BasicCompGeometry.turn_sign(
            k::_SweepCountKernel,
            p::Point{2},
            q::Point{2},
            r::Point{2},
        ) = (k.turns += 1; turn_sign(FloatingKernel(), p, q, r))
        BasicCompGeometry.compare_xy(k::_SweepCountKernel, p::Point{2}, q::Point{2}) =
            (k.cmps += 1; compare_xy(FloatingKernel(), p, q))

        function _horizontal_comb(m_prongs::Int)
            pts = Point2F[point(0.0, 0.0)]
            for i in 1:m_prongs
                y_bot = 2.0 * (i - 1)
                x_notch = 1.0 + (i - 1) / (m_prongs + 1)
                push!(pts, point(10.0, y_bot))
                push!(pts, point(10.0, y_bot + 1.0))
                push!(pts, point(x_notch, y_bot + 1.2))
                push!(pts, point(x_notch, y_bot + 1.8))
            end
            push!(pts, point(10.0, 2.0 * m_prongs))
            push!(pts, point(10.0, 2.0 * m_prongs + 1.0))
            push!(pts, point(0.0, 2.0 * m_prongs + 1.0))
            return pts
        end

        comb_s = _horizontal_comb(250)   # n = 1004
        comb_l = _horizontal_comb(2000)  # n = 8004
        ks = _SweepCountKernel(0, 0)
        kl = _SweepCountKernel(0, 0)
        rs = triangulate_sweep(comb_s; kernel = ks)
        rl = triangulate_sweep(comb_l; kernel = kl)
        @test length(rs.triangles) == length(comb_s) - 2
        @test length(rl.triangles) == length(comb_l) - 2
        @test isapprox(
            twice_area(comb_l),
            sum(
                twice_area([comb_l[a], comb_l[b], comb_l[c]]) for (a, b, c) in rl.triangles
            );
            rtol = 1e-8,
        )
        ns, nl = length(comb_s), length(comb_l)
        @test ks.turns / (ns * log2(ns)) < 6.0
        @test kl.turns / (nl * log2(nl)) < 6.0
        @test kl.cmps / (nl * log2(nl)) < 6.0
    end

    # Full structural check shared by the ear-cut testsets (counts, crossings, edge use, CCW, area).
    function _verify_triangulation_result(poly, res; kernel = FloatingKernel())
        m = length(poly)
        @test length(res.diagonals) == m - 3
        @test length(res.triangles) == m - 2
        @test length(unique(res.diagonals)) == m - 3
        @test all(d -> 1 <= d[1] < d[2] <= m && 1 < abs(d[1] - d[2]) < m - 1, res.diagonals)

        edges = [Segment(poly[i], poly[mod1(i + 1, m)]) for i in 1:m]
        diag_segs = [Segment(poly[u], poly[v]) for (u, v) in res.diagonals]
        no_diag_crossing = true
        no_edge_crossing = true
        for i in 1:length(diag_segs)
            di = diag_segs[i]
            for j in (i+1):length(diag_segs)
                if proper_crossing(di, diag_segs[j]; kernel = kernel)
                    no_diag_crossing = false
                    break
                end
            end
            for e in edges
                if proper_crossing(di, e; kernel = kernel)
                    no_edge_crossing = false
                    break
                end
            end
        end
        @test no_diag_crossing
        @test no_edge_crossing

        edge_counts = Dict{Tuple{Int,Int},Int}()
        all_ccw = true
        tri_sum = zero(twice_area([poly[1], poly[2], poly[3]]; kernel = kernel))
        for (a, b, c) in res.triangles
            pa, pb, pc = poly[a], poly[b], poly[c]
            orient_sign(kernel, pa, pb, pc) > 0 || (all_ccw = false)
            tri_sum += twice_area([pa, pb, pc]; kernel = kernel)
            for (u, v) in ((a, b), (b, c), (c, a))
                edge = u < v ? (u, v) : (v, u)
                edge_counts[edge] = get(edge_counts, edge, 0) + 1
            end
        end
        @test all_ccw
        @test all(
            i -> get(edge_counts, (min(i, mod1(i + 1, m)), max(i, mod1(i + 1, m))), 0) == 1,
            1:m,
        )
        @test all(d -> get(edge_counts, d, 0) == 2, res.diagonals)
        if kernel isa ExactKernel
            @test tri_sum == twice_area(poly; kernel = kernel)
        else
            @test isapprox(tri_sum, twice_area(poly; kernel = kernel); rtol = 1e-8)
        end
    end

    mutable struct _EarCutKdTreeCountKernel{K<:AbstractKernel} <: AbstractKernel
        inner::K
        orients::Int
    end
    BasicCompGeometry.orient_sign(
        k::_EarCutKdTreeCountKernel,
        a::Point{2},
        b::Point{2},
        c::Point{2},
    ) = (k.orients += 1; orient_sign(k.inner, a, b, c))
    BasicCompGeometry.turn_sign(
        k::_EarCutKdTreeCountKernel,
        a::Point{2},
        b::Point{2},
        c::Point{2},
    ) = (k.orients += 1; turn_sign(k.inner, a, b, c))

    for (ts_name, alg_sym, alg_fn) in (
        (
            "Ear-Cut Polygon Triangulation with the Kd-tree index (triangulate_ear_cut)",
            :ear_cut,
            triangulate_ear_cut,
        ),
        (
            "Unshrunk Kd-Tree Ear-Cut Polygon Triangulation (triangulate_ear_cut_u_kd_tree)",
            :ear_cut_u_kd_tree,
            triangulate_ear_cut_u_kd_tree,
        ),
    )
        @testset "$ts_name" begin
            @testset "Validation, dispatch, and edge cases" begin
                @test_throws ArgumentError alg_fn(Point2F[point(0.0, 0.0), point(1.0, 0.0)])
                cw_sq = Point2F[
                    point(0.0, 0.0),
                    point(0.0, 1.0),
                    point(1.0, 1.0),
                    point(1.0, 0.0),
                ]
                @test_throws ArgumentError alg_fn(cw_sq)
                @test_throws ArgumentError triangulate_polygon(cw_sq; algorithm = alg_sym)

                tri3 = Point2F[point(0.0, 0.0), point(2.0, 0.0), point(0.5, 1.5)]
                res3 = alg_fn(tri3)
                @test isempty(res3.diagonals)
                @test res3.triangles == [(1, 2, 3)]

                u_poly = Point2F[
                    point(0.0, 0.0),
                    point(4.0, 0.0),
                    point(4.0, 4.0),
                    point(3.0, 4.0),
                    point(3.0, 1.0),
                    point(1.0, 1.0),
                    point(1.0, 4.0),
                    point(0.0, 4.0),
                ]
                res_kd = alg_fn(u_poly)
                @test triangulate_polygon(u_poly; algorithm = alg_sym).diagonals ==
                      res_kd.diagonals
                @test triangulate_polygon(PntSeq(u_poly); algorithm = alg_sym).diagonals ==
                      res_kd.diagonals
                @test polygon_diagonals(u_poly; algorithm = alg_sym) == res_kd.diagonals
                _verify_triangulation_result(u_poly, res_kd)
            end

            @testset "ExactKernel, extreme scales, and collinear chains" begin
                exact = ExactKernel()
                collinear_poly = [
                    point(big(0) // 1, big(0) // 1),
                    point(big(1) // 1, big(0) // 1),
                    point(big(2) // 1, big(0) // 1),
                    point(big(3) // 1, big(0) // 1),
                    point(big(5) // 2, big(1) // 1),
                    point(big(2) // 1, big(2) // 1),
                    point(big(3) // 2, big(3) // 1),
                    point(big(1) // 1, big(2) // 1),
                    point(big(1) // 2, big(1) // 1),
                ]
                _verify_triangulation_result(
                    collinear_poly,
                    alg_fn(collinear_poly; kernel = exact);
                    kernel = exact,
                )

                exact_comb = Point{2,Rational{BigInt}}[]
                for i in 1:25
                    y0 = Rational{BigInt}(2 * (i - 1))
                    push!(exact_comb, point(big(10) // 1, y0))
                    push!(exact_comb, point(big(10) // 1, y0 + 1))
                    push!(exact_comb, point(big(1) // 1 + (i % 3) // 2, y0 + 6 // 5))
                    push!(exact_comb, point(big(1) // 1 + (i % 3) // 2, y0 + 9 // 5))
                end
                pushfirst!(exact_comb, point(big(0) // 1, big(0) // 1))
                push!(
                    exact_comb,
                    point(big(10) // 1, big(50) // 1),
                    point(big(0) // 1, big(50) // 1),
                )
                _verify_triangulation_result(
                    exact_comb,
                    alg_fn(exact_comb; kernel = exact);
                    kernel = exact,
                )

                huge = big(10)^320
                huge_comb =
                    [Point{2,Rational{BigInt}}(p.x * huge, p.y * huge) for p in exact_comb]
                _verify_triangulation_result(
                    huge_comb,
                    alg_fn(huge_comb; kernel = exact);
                    kernel = exact,
                )

                if alg_sym == :ear_cut_u_kd_tree
                    int_comb = Point{2,Int}[point(0, 0)]
                    for i in 1:30
                        y0 = 10 * (i - 1)
                        xn = 2 + (i % 4)
                        push!(
                            int_comb,
                            point(20, y0),
                            point(20, y0 + 4),
                            point(xn, y0 + 5),
                            point(xn, y0 + 9),
                        )
                    end
                    push!(int_comb, point(20, 300), point(0, 300))
                    _verify_triangulation_result(
                        int_comb,
                        alg_fn(int_comb; kernel = exact);
                        kernel = exact,
                    )
                end

                for s in (1e-12, 1e12)
                    comb_s = [point(Float64(p.x) * s, Float64(p.y) * s) for p in exact_comb]
                    _verify_triangulation_result(comb_s, alg_fn(comb_s))
                end
            end

            @testset "Subquadratic practical performance vs ear_cut" begin
                poly_conv =
                    Point2F[point(cos(2pi * i / 2000), sin(2pi * i / 2000)) for i in 0:1999]
                kc = _EarCutKdTreeCountKernel(FloatingKernel(), 0)
                res_conv = triangulate_polygon(poly_conv; kernel = kc, algorithm = alg_sym)
                @test length(res_conv.triangles) == length(poly_conv) - 2
                @test kc.orients / length(poly_conv) < 5.0

                poly_lowr = Point2F[
                    point(
                        (i == 1000 || i == 3000 ? 0.995 : 1.0) * cos(2pi * i / 4000),
                        (i == 1000 || i == 3000 ? 0.995 : 1.0) * sin(2pi * i / 4000),
                    ) for i in 0:3999
                ]
                klr = _EarCutKdTreeCountKernel(FloatingKernel(), 0)
                res_lr = triangulate_polygon(poly_lowr; kernel = klr, algorithm = alg_sym)
                @test length(res_lr.triangles) == 3998
                @test klr.orients / length(poly_lowr) < 20.0

                function _hcomb_kd(m_prongs::Int)
                    pts = Point2F[point(0.0, 0.0)]
                    for i in 1:m_prongs
                        yb = 2.0 * (i - 1)
                        xn = 1.2 + 1.4 * (i % 2) + 0.5 * sin(pi * i / (m_prongs + 1))
                        push!(
                            pts,
                            point(10.0, yb),
                            point(10.0, yb + 1.0),
                            point(xn, yb + 1.2),
                            point(xn, yb + 1.8),
                        )
                    end
                    push!(pts, point(10.0, 2.0 * m_prongs), point(0.0, 2.0 * m_prongs))
                    return pts
                end
                comb_1k = _hcomb_kd(300)
                comb_5k = _hcomb_kd(1200)
                kec1 = _EarCutKdTreeCountKernel(FloatingKernel(), 0)
                kec5 = _EarCutKdTreeCountKernel(FloatingKernel(), 0)
                kcl1 = _EarCutKdTreeCountKernel(FloatingKernel(), 0)

                r_ec1 = triangulate_polygon(comb_1k; kernel = kec1, algorithm = alg_sym)
                r_ec5 = triangulate_polygon(comb_5k; kernel = kec5, algorithm = alg_sym)
                BasicCompGeometry._ear_clipping_from_valid(
                    comb_1k,
                    length(comb_1k);
                    kernel = kcl1,
                )

                _verify_triangulation_result(comb_1k, r_ec1)
                @test length(r_ec5.triangles) == length(comb_5k) - 2
                @test isapprox(
                    twice_area(comb_5k),
                    sum(
                        twice_area([comb_5k[a], comb_5k[b], comb_5k[c]]) for
                        (a, b, c) in r_ec5.triangles
                    );
                    rtol = 1e-8,
                )
                @test kec1.orients / length(comb_1k) < 35.0
                @test kec5.orients / length(comb_5k) < 35.0
                @test (kec5.orients / length(comb_5k)) - (kec1.orients / length(comb_1k)) <
                      8.0
                @test kec1.orients * 20 < kcl1.orients

                rng_s = Xoshiro(77)
                star = [
                    let r = iseven(i) ? 10.0 : 1.5 + 1.5 * rand(rng_s)
                        point(r * cos(2pi * i / 1000), r * sin(2pi * i / 1000))
                    end for i in 0:999
                ]
                spiral = Point2F[]
                for k in 0:512
                    r = 1.0 + 0.5 * k / 32
                    push!(spiral, point(r * cos(2pi * k / 32), r * sin(2pi * k / 32)))
                end
                for k in 512:-1:0
                    r = 0.7 + 0.5 * k / 32
                    push!(spiral, point(r * cos(2pi * k / 32), r * sin(2pi * k / 32)))
                end
                coast = [
                    let th = 2pi * i / 2000
                        r =
                            10.0 +
                            2.0 * sin(13 * th) +
                            1.1 * cos(31 * th) +
                            0.6 * sin(73 * th)
                        point(r * cos(th), r * sin(th))
                    end for i in 0:1999
                ]
                _verify_triangulation_result(
                    star,
                    triangulate_polygon(star; algorithm = alg_sym),
                )
                _verify_triangulation_result(
                    spiral,
                    triangulate_polygon(spiral; algorithm = alg_sym),
                )
                _verify_triangulation_result(
                    coast,
                    triangulate_polygon(coast; algorithm = alg_sym),
                )
            end
        end
    end

    @testset "Deprecated ear-cut variants forward to the survivors" begin
        rng = Xoshiro(5)
        P = [
            let r = iseven(i) ? 2.0 : 0.5 + rand(rng)
                point(r * cos(2pi * i / 400), r * sin(2pi * i / 400))
            end for i in 0:399
        ]
        z = @test_deprecated triangulate_ear_cut_z_order(P)
        @test z.triangles == triangulate_ear_cut(P).triangles
        kd = @test_deprecated triangulate_ear_cut_kd_tree(P)
        @test kd.triangles == triangulate_ear_cut_u_kd_tree(P).triangles
        z = @test_deprecated triangulate_polygon(P; algorithm = :ear_cut_z_order)
        @test z.diagonals == triangulate_ear_cut(P).diagonals
        kd = @test_deprecated polygon_diagonals(P; algorithm = :ear_cut_kd_tree)
        @test kd == triangulate_ear_cut_u_kd_tree(P).diagonals
    end

    @testset "Degenerate input handling and cross-algorithm agreement" begin
        all_algs = (:sweep, :seidel, :ear_cut, :ear_cut_u_kd_tree)
        sq = [point(0.0, 0.0), point(1.0, 0.0), point(1.0, 1.0), point(0.0, 1.0)]

        # `:ear_clipping` (pre-0.5.15 name) stays accepted as an alias of `:ear_cut`.
        @test triangulate_polygon(sq; algorithm = :ear_clipping).triangles ==
              triangulate_polygon(sq; algorithm = :ear_cut).triangles

        @test occursin("Seidel", string(@doc(triangulate_seidel)))

        # Repeated vertices (consecutive, closing repeat, pinch) are rejected, never truncated.
        dup_polys = (
            [
                point(0.0, 0.0),
                point(1.0, 0.0),
                point(1.0, 0.0),
                point(1.0, 1.0),
                point(0.0, 1.0),
            ],
            [
                point(0.0, 0.0),
                point(2.0, 0.0),
                point(2.0, 2.0),
                point(0.0, 2.0),
                point(0.0, 0.0),
            ],
            [
                point(0.0, 0.0),
                point(2.0, 2.0),
                point(4.0, 0.0),
                point(4.0, 4.0),
                point(2.0, 2.0),
                point(0.0, 4.0),
            ],
        )
        for P in dup_polys, alg in all_algs
            @test_throws ArgumentError triangulate_polygon(P; algorithm = alg)
        end

        # A simple polygon that is numerically degenerate for Float64 predicates: ear clipping
        # must throw instead of returning a partial triangulation; ExactKernel succeeds.
        q14 = [
            (0.1, 0.3),
            (0.011707525606919664, 0.3469514532847201),
            (-0.42975484635848205, 0.5817087197083206),
            (-0.3884138252501218, 0.4464647920305202),
            (-0.5706092062126423, 0.3168312965290796),
            (-0.34707280414176156, 0.31122086435271973),
            (-0.3940242574264817, 0.22292838995963937),
            (-0.48792716399592195, 0.04634344117347866),
            (-0.3526832363181214, 0.08768446228183896),
            (-0.534878617280642, -0.04194903321960164),
            (-0.18170871970832062, -0.22975484635848203),
            (-0.08780581313888042, -0.05316989757232138),
            (-0.2643907619250411, 0.04073300899711885),
            (-0.25878032974868126, 0.2642694110679996),
        ]
        P14 = [point(x, y) for (x, y) in q14]
        for alg in (:ear_cut, :ear_cut_u_kd_tree)
            r = try
                triangulate_polygon(P14; algorithm = alg)
            catch e
                e
            end
            @test r isa ArgumentError || length(r.triangles) == length(P14) - 2
            @test length(
                triangulate_polygon(P14; algorithm = alg, kernel = ExactKernel()).triangles,
            ) == 12
        end

        # Seidel walks fail with a clear ArgumentError (not BoundsError) on inconsistent
        # floating-point predicates; the same points triangulate exactly with ExactKernel.
        k19 = [
            (1, 2),
            (0, 2),
            (1, 0),
            (2, 0),
            (6, 0),
            (6, 1),
            (5, 2),
            (5, 1),
            (3, 5),
            (5, 5),
            (6, 6),
            (3, 6),
            (1, 5),
            (0, 5),
            (1, 4),
            (3, 2),
            (4, 1),
            (2, 1),
            (0, 3),
        ]
        P19 = [point(a * (1 / 3), b * (1 / 3)) for (a, b) in k19]
        r19 = try
            triangulate_seidel(P19; rng = Xoshiro(1))
        catch e
            e
        end
        @test r19 isa ArgumentError || length(r19.triangles) == 17
        _verify_triangulation_result(
            P19,
            triangulate_seidel(P19; kernel = ExactKernel(), rng = Xoshiro(1));
            kernel = ExactKernel(),
        )

        # Randomized differential test on simple polygons over small integer grids (many
        # collinear vertices, shared coordinates, horizontal/vertical edges), exact predicates.
        _o(a, b, c) = sign((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x))
        _onseg(a, b, p) =
            _o(a, b, p) == 0 &&
            min(a.x, b.x) <= p.x <= max(a.x, b.x) &&
            min(a.y, b.y) <= p.y <= max(a.y, b.y)
        function _meet(a, b, c, d)
            o1, o2, o3, o4 = _o(a, b, c), _o(a, b, d), _o(c, d, a), _o(c, d, b)
            (o1 * o2 < 0 && o3 * o4 < 0) && return true
            return _onseg(a, b, c) || _onseg(a, b, d) || _onseg(c, d, a) || _onseg(c, d, b)
        end
        function _simple(P)
            m = length(P)
            for i in 1:m, j in (i+1):m
                a, b, c, d = P[i], P[mod1(i + 1, m)], P[j], P[mod1(j + 1, m)]
                if j == i + 1
                    (_onseg(c, d, a) || _onseg(a, b, d)) && return false
                elseif i == 1 && j == m
                    (_onseg(c, d, b) || _onseg(a, b, c)) && return false
                else
                    _meet(a, b, c, d) && return false
                end
            end
            return true
        end
        function _rand_simple_polygon(rng, n, G)
            P = unique([
                point(Float64(rand(rng, 0:G)), Float64(rand(rng, 0:G))) for _ in 1:n
            ])
            length(P) < 3 && return nothing
            shuffle!(rng, P)
            m = length(P)
            for _ in 1:500  # 2-opt untangling
                changed = false
                for i in 1:m, j in (i+2):m
                    (i == 1 && j == m) && continue
                    if _meet(P[i], P[i+1], P[j], P[mod1(j + 1, m)])
                        reverse!(P, i + 1, j)
                        changed = true
                    end
                end
                changed || break
            end
            _simple(P) || return nothing
            a2 = twice_area(P; kernel = ExactKernel())
            a2 == 0 && return nothing
            a2 < 0 && reverse!(P)
            return P
        end

        rng = Xoshiro(20261113)
        exact = ExactKernel()
        tested = 0
        while tested < 120
            P = _rand_simple_polygon(rng, rand(rng, 4:64), rand(rng, (4, 8, 16, 30)))
            P === nothing && continue
            tested += 1
            ref = BasicCompGeometry._ear_clipping_from_valid(P, length(P); kernel = exact)
            _verify_triangulation_result(P, ref; kernel = exact)
            _verify_triangulation_result(
                P,
                triangulate_ear_cut(P; kernel = exact);
                kernel = exact,
            )
            for alg in (:sweep, :seidel)
                _verify_triangulation_result(
                    P,
                    triangulate_polygon(
                        P;
                        algorithm = alg,
                        kernel = exact,
                        rng = Xoshiro(tested),
                    );
                    kernel = exact,
                )
            end
            # The spatial indexes only accelerate the same blocking test, so the indexed ear
            # cutters must reproduce the plain reflex-list ear clipper exactly (polygons with
            # more than 16 reflex vertices exercise the Z-order and Kd-tree indexes).
            for alg in (:ear_cut_u_kd_tree,)
                @test triangulate_polygon(P; algorithm = alg, kernel = exact).triangles ==
                      ref.triangles
            end
        end
    end

    @testset "Performance review #001 regressions" begin
        # Finding 1 & 2: Stack buffer reuse and direct CCW emission in triangulate_monotone
        z = vcat(
            [point(0.0, 0.0)],
            [point(2.0 * i - 1.0, 1.0) for i in 1:1000],
            [point(2001.0, 0.0)],
            [point(2.0 * i, -1.0) for i in 1000:-1:1],
        )
        res_float = triangulate_monotone(z; kernel = FloatingKernel())
        @test length(res_float.triangles) == length(z) - 2
        @test all(
            t -> orient_sign(FloatingKernel(), z[t[1]], z[t[2]], z[t[3]]) > 0,
            res_float.triangles,
        )

        # ExactKernel: all emitted triangles must be strictly CCW
        res_exact = triangulate_monotone(z; kernel = ExactKernel())
        @test length(res_exact.triangles) == length(z) - 2
        @test all(
            t -> orient_sign(ExactKernel(), z[t[1]], z[t[2]], z[t[3]]) > 0,
            res_exact.triangles,
        )

        # Stack buffer reuse prevents O(k) heap array reallocations on alternating chains
        triangulate_monotone(z; kernel = FloatingKernel())
        allocs_float = @allocated triangulate_monotone(z; kernel = FloatingKernel())
        @test allocs_float < 200_000

        # Finding 3: _check_distinct_vertices allocation and PntSeq support
        pts_20k =
            [point(cos(th), sin(th)) for th in range(0, 2pi, length = 20000)[1:(end-1)]]
        @test BasicCompGeometry._check_distinct_vertices(FloatingKernel(), pts_20k) ===
              nothing
        BasicCompGeometry._check_distinct_vertices(FloatingKernel(), pts_20k)
        allocs_20k =
            @allocated BasicCompGeometry._check_distinct_vertices(FloatingKernel(), pts_20k)
        @test allocs_20k < 300_000

        pseq_20k = PntSeq(pts_20k)
        @test BasicCompGeometry._check_distinct_vertices(FloatingKernel(), pseq_20k) ===
              nothing

        dup_pts = [point(0.0, 0.0), point(1.0, 0.0), point(1.0, 1.0), point(0.0, 0.0)]
        @test_throws ArgumentError BasicCompGeometry._check_distinct_vertices(
            FloatingKernel(),
            dup_pts,
        )
        @test_throws ArgumentError BasicCompGeometry._check_distinct_vertices(
            FloatingKernel(),
            PntSeq(dup_pts),
        )

        # Finding 4: Sweep partition without prev_v array allocation and with sizehinted diags
        poly_sweep =
            [point(cos(th), sin(th)) for th in range(0, 2pi, length = 200)[1:(end-1)]]
        res_sweep = triangulate_sweep(poly_sweep)
        @test length(res_sweep.triangles) == length(poly_sweep) - 2
        @test length(res_sweep.diagonals) == length(poly_sweep) - 3
        diags = BasicCompGeometry._sweep_partition_diagonals(FloatingKernel(), poly_sweep)
        @test diags isa Vector{Tuple{Int,Int}}
    end

    @testset "Weakly simple and near-degenerate polygons: valid triangulation or ArgumentError" begin
        all_algs = (:sweep, :seidel, :ear_cut, :ear_cut_u_kd_tree)
        exact = ExactKernel()
        tri_orient(P, t, kernel) = orient_sign(kernel, P[t[1]], P[t[2]], P[t[3]])

        # Each polygon edge is used once in CCW direction, every diagonal once each way.
        function edges_ok(n, triangles)
            directed = Dict{Tuple{Int,Int},Int}()
            for t in triangles, e in ((t[1], t[2]), (t[2], t[3]), (t[3], t[1]))
                directed[e] = get(directed, e, 0) + 1
            end
            all(c == 1 for c in values(directed)) || return false
            return all(
                e[2] == mod1(e[1] + 1, n) || haskey(directed, (e[2], e[1])) for
                e in keys(directed)
            )
        end

        # The contract: no triangle is clockwise under `kernel`, the triangles form a
        # triangulation of the boundary cycle and, with exact predicates, tile the polygon.
        function contract_ok(P, res, kernel)
            n = length(P)
            length(res.triangles) == n - 2 || return false
            all(t -> tri_orient(P, t, kernel) >= 0, res.triangles) || return false
            edges_ok(n, res.triangles) || return false
            kernel isa ExactKernel || return true
            s = sum(
                t -> abs(twice_area([P[t[1]], P[t[2]], P[t[3]]]; kernel = exact)),
                res.triangles,
            )
            return s == twice_area(P; kernel = exact)
        end

        # :ok, :throw (ArgumentError), or :violation.
        function outcome(P, alg, kernel)
            res = try
                triangulate_polygon(P; algorithm = alg, kernel = kernel, rng = Xoshiro(7))
            catch e
                e isa ArgumentError && return :throw
                rethrow()
            end
            return contract_ok(P, res, kernel) ? :ok : :violation
        end

        # Reported repro: vertex 4 is within rounding of edge 2 -> 3 (exactly on it in
        # decimal; exact orient(2, 3, 4) = +1, so the polygon is strictly simple). The y-monotone
        # fast path of :ear_cut rotated it and the local turn test that split the monotone
        # chains rounded to 0, swapping them: triangles [(3, 1, 4), (2, 1, 3)], one clockwise.
        P = [
            point(0.74, 0.21),
            point(-0.35, 0.56),
            point(-0.03, -0.56),
            point(-0.07, -0.42),
        ]
        @test orient_sign(exact, P[2], P[3], P[4]) > 0
        for alg in all_algs, kernel in (FloatingKernel(), exact)
            @test outcome(P, alg, kernel) === :ok
        end
        R = [point(p.y, -p.x) for p in P]
        @test is_x_monotone(R)
        @test contract_ok(R, triangulate_monotone(R), FloatingKernel())
        @test monotone_chain_indices(R) == ([3, 2], [3, 4, 1, 2])

        # Exact zero-width spike at the leftmost vertex: both chains leave it along the
        # same ray, so the old local turn test was exactly zero even with ExactKernel.
        S = [
            point(0.0, 0.0),
            point(2.0, 1.0),
            point(4.0, 0.0),
            point(6.0, 2.0),
            point(4.0, 5.0),
            point(1.0, 0.5),
        ]
        @test monotone_chain_indices(S; kernel = exact) == ([1, 6, 5, 4], [1, 2, 3, 4])
        for kernel in (FloatingKernel(), exact)
            @test contract_ok(S, triangulate_monotone(S; kernel = kernel), kernel)
            @test contract_ok(S, triangulate_ear_cut(S; kernel = kernel), kernel)
        end
        # Clockwise x-monotone input still triangulates into CCW triangles.
        Sr = reverse(S)
        res_r = triangulate_monotone(Sr)
        @test length(res_r.triangles) == length(Sr) - 2
        @test all(t -> tri_orient(Sr, t, exact) >= 0, res_r.triangles)

        # Vertex 3 touches edge 6 -> 7: Seidel's piece triangulation swapped the chains of a
        # piece with a zero-width spike and returned overlapping triangles.
        W = [
            point(-54.0, -42.0),
            point(-36.0, -44.0),
            point(-12.0, 0.0),
            point(-18.0, -46.0),
            point(0.0, -48.0),
            point(24.0, 0.0),
            point(-48.0, 0.0),
        ]
        for alg in all_algs, kernel in (FloatingKernel(), exact), seed in 1:5
            r = try
                triangulate_polygon(W; algorithm = alg, kernel = kernel, rng = Xoshiro(seed))
            catch e
                e
            end
            @test r isa ArgumentError || contract_ok(W, r, kernel)
        end

        # Strictly simple, but FloatingKernel predicates are inconsistent: ear clipping used
        # to emit the clockwise remaining triangle (4, 5, 7) unchecked.
        F = [
            point(-0.69, 0.59),
            point(0.89, -0.89),
            point(0.59, 0.69),
            point(0.46, 0.69),
            point(0.1, -0.15),
            point(0.33, 0.69),
            point(0.2, 0.69),
        ]
        for alg in all_algs
            @test outcome(F, alg, FloatingKernel()) !== :violation
            @test outcome(F, alg, exact) === :ok
        end

        # x-monotone vertex order whose chains cross (not simple): the monotone fast path
        # of :ear_cut returned clockwise triangles; now every algorithm rejects it.
        X = [
            point(0.0, 0.0),
            point(2.0, -1.0),
            point(4.0, 2.0),
            point(6.0, 0.0),
            point(3.0, 1.0),
        ]
        @test is_x_monotone(X) && twice_area(X; kernel = exact) > 0
        for alg in all_algs, kernel in (FloatingKernel(), exact)
            @test outcome(X, alg, kernel) !== :violation
        end
        @test_throws ArgumentError triangulate_monotone(X)

        # Seeded weakly simple polygons from random convex lattice polygons (coordinates
        # are multiples of 6, so midpoints and thirds stay integral and Float64 predicates
        # are exact): (a) a zero-width spike, a vertex inserted on its predecessor edge;
        # (b) a notch whose apex touches another edge. Ear cutting handles both.
        lerp(p, q, f) = point(p.x + f * (q.x - p.x), p.y + f * (q.y - p.y))
        function rand_convex_lattice(rng)
            pts = unique([
                point(6.0 * rand(rng, -10:10), 6.0 * rand(rng, -10:10)) for _ in 1:10
            ])
            H = [point(p[1], p[2]) for p in convex_hull(pts)]
            twice_area(H; kernel = exact) < 0 && reverse!(H)
            m = length(H)
            return [
                H[i] for i in 1:m if
                orient_sign(exact, H[mod1(i - 1, m)], H[i], H[mod1(i + 1, m)]) > 0
            ]
        end
        function weakly_simple(rng, H, spike::Bool)
            n = length(H)
            i = rand(rng, 1:n)
            spike && return vcat(H[1:i], [lerp(H[i], H[mod1(i - 1, n)], 0.5)], H[(i+1):end])
            j = mod1(i + rand(rng, 2:(n-1)), n)
            t = lerp(H[i], H[mod1(i + 1, n)], 0.5)
            x, y =
                lerp(H[j], H[mod1(j + 1, n)], 1 / 3), lerp(H[j], H[mod1(j + 1, n)], 2 / 3)
            return vcat(H[1:j], [x, t, y], H[(j+1):end])
        end
        counts = Dict{Symbol,Int}(:ok => 0, :throw => 0, :violation => 0)
        tested = 0
        for seed in 1:120
            rng = Xoshiro(seed)
            H = rand_convex_lattice(rng)
            length(H) >= 4 || continue
            Q = weakly_simple(rng, H, isodd(seed))
            allunique(Q) || continue
            tested += 1
            for alg in all_algs, kernel in (FloatingKernel(), exact)
                o = outcome(Q, alg, kernel)
                counts[o] += 1
                alg in (:sweep, :seidel) || @test o === :ok
            end
        end
        @test tested >= 100
        @test counts[:violation] == 0
    end
end
