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
        @test res100.diagonals == triangulate_sweep(star100).diagonals
        @test length(res100.diagonals) == n - 3
        @test length(res100.triangles) == n - 2
        poly_area100 = twice_area(star100)
        tri_sum100 = sum(twice_area([star100[t[1]], star100[t[2]], star100[t[3]]]) for t in res100.triangles)
        @test isapprox(poly_area100, tri_sum100; rtol = 1e-10)

        res100_ear = triangulate_polygon(star100; algorithm = :ear_cut)
        @test triangulate_ear_cut(star100).diagonals == res100_ear.diagonals
        @test_throws ArgumentError triangulate_ear_cut([point(0.0, 0.0), point(1.0, 0.0)])
        @test_throws ArgumentError triangulate_ear_cut([point(0.0, 0.0), point(0.0, 1.0), point(1.0, 0.0)])
        @test triangulate_ear_cut([point(0.0, 0.0), point(1.0, 0.0), point(0.0, 1.0)]).triangles == [(1, 2, 3)]
        @test length(res100_ear.diagonals) == n - 3
        @test length(res100_ear.triangles) == n - 2
        tri_sum100_ear = sum(twice_area([star100[t[1]], star100[t[2]], star100[t[3]]]) for t in res100_ear.triangles)
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

        @testset "Trapezoidal map invariants, spirals, collinearity, and O(n log* n) scaling" begin
            function _verify_seidel_result(poly, res; kernel = FloatingKernel())
                m = length(poly)
                @test length(res.diagonals) == m - 3
                @test length(res.triangles) == m - 2
                @test length(unique(res.diagonals)) == m - 3
                edges = [Segment(poly[i], poly[mod1(i + 1, m)]) for i in 1:m]
                diag_segs = [Segment(poly[u], poly[v]) for (u, v) in res.diagonals]
                for i in 1:length(diag_segs)
                    for j in i+1:length(diag_segs)
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
                tri_a = sum(twice_area([poly[a], poly[b], poly[c]]; kernel = kernel) for (a, b, c) in res.triangles)
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
                    @test BasicCompGeometry._seidel_cmp_pt(kernel, m.pts, T.leftp, T.rightp) < 0
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
            vquad = [point(4//1, 5//1), point(0//1, 4//1), point(0//1, 0//1), point(4//1, 0//1)]
            _verify_trap_invariants(vquad, Xoshiro(7); kernel = exact)
            _verify_seidel_result(vquad, triangulate_seidel(vquad; kernel = exact, rng = Xoshiro(7)); kernel = exact)

            collinear_tri = [
                point(0//1, 0//1), point(1//1, 0//1), point(2//1, 0//1), point(3//1, 0//1),
                point(5//2, 1//1), point(2//1, 2//1), point(3//2, 3//1), point(1//1, 2//1),
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
                _verify_seidel_result(spiral, triangulate_seidel(spiral; rng = Xoshiro(seed)))
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
            BasicCompGeometry.turn_sign(k::_SeidelCountKernel, p::Point{2}, q::Point{2}, r::Point{2}) =
                (k.turns += 1; turn_sign(FloatingKernel(), p, q, r))
            BasicCompGeometry.compare_xy(k::_SeidelCountKernel, p::Point{2}, q::Point{2}) =
                (k.cmps += 1; compare_xy(FloatingKernel(), p, q))

            function _double_comb(m_teeth::Int)
                pts = Point2F[point(0.0, 0.0), point(20.0, 0.0)]
                for i in 1:m_teeth
                    push!(pts, point(20.0, 2.0 * i - 1.0), point(12.0, 2.0 * i))
                end
                push!(pts, point(20.0, 2.0 * m_teeth + 1.0), point(0.0, 2.0 * m_teeth + 1.0))
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
                sum(twice_area([comb_large[a], comb_large[b], comb_large[c]]) for (a, b, c) in res_large.triangles);
                rtol = 1e-8,
            )
            # Predicate work per vertex remains bounded by a small O(log* n) constant
            @test kc_small.turns / length(comb_small) < 20.0
            @test kc_large.turns / length(comb_large) < 20.0
            @test kc_large.cmps / length(comb_large) < 25.0

            # Ear clipping with reflex-vertex list stays O(n * (r + 1)) rather than O(n^2)
            function _non_monotone_few_reflex(n::Int)
                pts = Point2F[]
                for i in 0:(n - 3)
                    th = 0.2 + 1.8 * pi * i / (n - 2)
                    push!(pts, point(cos(th), sin(th)))
                end
                push!(pts, point(0.6, 0.2), point(0.7, -0.3))
                return pts
            end
            for n_ear in (500, 4000)
                poly_ear = _non_monotone_few_reflex(n_ear)
                @test !is_x_monotone(poly_ear) && !is_y_monotone(poly_ear)
                ke = _SeidelCountKernel(0, 0)
                res_ear = triangulate_polygon(poly_ear; kernel = ke, algorithm = :ear_cut)
                @test length(res_ear.triangles) == n_ear - 2
                @test ke.turns / n_ear < 15.0
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

            dup_pts = [
                point(0.0, 0.0),
                point(2.0, 0.0),
                point(2.0, 0.0),
                point(0.0, 2.0),
            ]
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
            bowtie = [
                point(0.0, 0.0),
                point(6.0, 0.0),
                point(2.0, 8.0),
                point(5.0, -2.0),
            ]
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
                next_v[2i - 1] = 2i
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
                got = BasicCompGeometry._avl_below(tree, FloatingKernel(), qpts, next_v, q_idx)
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
                for j in i+1:length(diag_segs)
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
                twice_area([poly[a], poly[b], poly[c]]; kernel = kernel) for (a, b, c) in res.triangles
            )
            if kernel isa ExactKernel
                @test poly_a == tri_a
            else
                @test isapprox(poly_a, tri_a; rtol = 1e-9)
            end
        end

        # Vertical edges, shared-x jagged reflex/merge chains, collinear chains, spirals, snakes
        vquad = [point(4 // 1, 5 // 1), point(0 // 1, 4 // 1), point(0 // 1, 0 // 1), point(4 // 1, 0 // 1)]
        _verify_sweep_result(vquad, triangulate_sweep(vquad; kernel = exact); kernel = exact)

        vcollinear = [
            point(0 // 1, 3 // 1),
            point(0 // 1, 2 // 1),
            point(0 // 1, 1 // 1),
            point(0 // 1, 0 // 1),
            point(2 // 1, 0 // 1),
            point(2 // 1, 3 // 1),
        ]
        _verify_sweep_result(vcollinear, triangulate_sweep(vcollinear; kernel = exact); kernel = exact)

        collinear_tri = [
            point(0 // 1, 0 // 1), point(1 // 1, 0 // 1), point(2 // 1, 0 // 1), point(3 // 1, 0 // 1),
            point(5 // 2, 1 // 1), point(2 // 1, 2 // 1), point(3 // 2, 3 // 1), point(1 // 1, 2 // 1),
            point(1 // 2, 1 // 1),
        ]
        _verify_sweep_result(collinear_tri, triangulate_sweep(collinear_tri; kernel = exact); kernel = exact)

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
            point(0 // 1, 0 // 1), point(4 // 1, 0 // 1), point(4 // 1, 4 // 1),
            point(3 // 1, 4 // 1), point(3 // 1, 1 // 1), point(1 // 1, 1 // 1),
            point(1 // 1, 4 // 1), point(0 // 1, 4 // 1),
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
        BasicCompGeometry.turn_sign(k::_SweepCountKernel, p::Point{2}, q::Point{2}, r::Point{2}) =
            (k.turns += 1; turn_sign(FloatingKernel(), p, q, r))
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
            sum(twice_area([comb_l[a], comb_l[b], comb_l[c]]) for (a, b, c) in rl.triangles);
            rtol = 1e-8,
        )
        ns, nl = length(comb_s), length(comb_l)
        @test ks.turns / (ns * log2(ns)) < 6.0
        @test kl.turns / (nl * log2(nl)) < 6.0
        @test kl.cmps / (nl * log2(nl)) < 6.0
    end

    @testset "Z-Order Ear-Cut Polygon Triangulation (triangulate_ear_cut_z_order)" begin
        mutable struct _EarCutZOrderCountKernel{K<:AbstractKernel} <: AbstractKernel
            inner::K
            orients::Int
        end
        BasicCompGeometry.orient_sign(k::_EarCutZOrderCountKernel, a::Point{2}, b::Point{2}, c::Point{2}) =
            (k.orients += 1; orient_sign(k.inner, a, b, c))
        BasicCompGeometry.turn_sign(k::_EarCutZOrderCountKernel, a::Point{2}, b::Point{2}, c::Point{2}) =
            (k.orients += 1; turn_sign(k.inner, a, b, c))

        function _verify_ear_cut_z_order_result(poly, res; kernel = FloatingKernel())
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
                for j in (i + 1):length(diag_segs)
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
            for (a, b, c) in res.triangles
                orient_sign(kernel, poly[a], poly[b], poly[c]) > 0 || (all_ccw = false)
                for (u, v) in ((a, b), (b, c), (c, a))
                    edge = u < v ? (u, v) : (v, u)
                    edge_counts[edge] = get(edge_counts, edge, 0) + 1
                end
            end
            @test all_ccw
            @test all(i -> get(edge_counts, (i < mod1(i + 1, m) ? (i, mod1(i + 1, m)) : (mod1(i + 1, m), i)), 0) == 1, 1:m)
            @test all(d -> get(edge_counts, d, 0) == 2, res.diagonals)

            poly_a = twice_area(poly; kernel = kernel)
            tri_a = sum(twice_area([poly[a], poly[b], poly[c]]; kernel = kernel) for (a, b, c) in res.triangles)
            if kernel isa ExactKernel
                @test poly_a == tri_a
            else
                @test isapprox(poly_a, tri_a; rtol = 1e-8)
            end
        end

        @testset "Preconditions and trivial triangle" begin
            @test_throws ArgumentError triangulate_ear_cut_z_order([point(0.0, 0.0), point(1.0, 0.0)])
            @test_throws ArgumentError triangulate_ear_cut_z_order(Point2F[])
            cw = [point(0.0, 0.0), point(0.0, 1.0), point(1.0, 0.0)]
            @test_throws ArgumentError triangulate_ear_cut_z_order(cw)
            @test_throws ArgumentError triangulate_polygon(cw; algorithm = :ear_cut_z_order)

            tri = [point(0.0, 0.0), point(2.0, 0.0), point(0.0, 2.0)]
            res_tri = triangulate_ear_cut_z_order(tri)
            @test isempty(res_tri.diagonals)
            @test res_tri.triangles == [(1, 2, 3)]
        end

        @testset "Morton BIGMIN/LITMAX and _ZOrderReflexIndex invariants" begin
            rng = Xoshiro(9182)
            for _ in 1:400
                x0 = rand(rng, UInt32(0):UInt32(65535))
                x1 = rand(rng, UInt32(0):UInt32(65535))
                y0 = rand(rng, UInt32(0):UInt32(65535))
                y1 = rand(rng, UInt32(0):UInt32(65535))
                x0 > x1 && ((x0, x1) = (x1, x0))
                y0 > y1 && ((y0, y1) = (y1, y0))
                lo = BasicCompGeometry._morton_code(x0, y0)
                hi = BasicCompGeometry._morton_code(x1, y1)
                z = rand(rng, UInt32(0):typemax(UInt32))
                min_zx = lo & BasicCompGeometry._MORTON_X_MASK
                max_zx = hi & BasicCompGeometry._MORTON_X_MASK
                min_zy = lo & BasicCompGeometry._MORTON_Y_MASK
                max_zy = hi & BasicCompGeometry._MORTON_Y_MASK

                ok_b, bm = BasicCompGeometry._morton_bigmin(z, lo, hi)
                if z <= hi
                    @test ok_b
                    @test bm >= z && BasicCompGeometry._in_morton_box(bm, min_zx, max_zx, min_zy, max_zy)
                    if !BasicCompGeometry._in_morton_box(z, min_zx, max_zx, min_zy, max_zy)
                        @test bm > z
                    end
                else
                    @test !ok_b
                end

                ok_l, lm = BasicCompGeometry._morton_litmax(z, lo, hi)
                if z >= lo
                    @test ok_l
                    @test lm <= z && BasicCompGeometry._in_morton_box(lm, min_zx, max_zx, min_zy, max_zy)
                    if !BasicCompGeometry._in_morton_box(z, min_zx, max_zx, min_zy, max_zy)
                        @test lm < z
                    end
                else
                    @test !ok_l
                end
            end

            # Verify _ZOrderReflexIndex bucket, bitmask, and doubly-linked list invariants
            n_idx = 200
            pts_idx = [point(cos(2pi * i / n_idx), sin(2pi * i / n_idx)) for i in 0:(n_idx - 1)]
            refl = falses(n_idx)
            for i in 1:3:n_idx
                refl[i] = true
            end
            ok_b, min_x, min_y, inv_x, inv_y = BasicCompGeometry._compute_z_bounds(pts_idx)
            @test ok_b
            zidx = BasicCompGeometry._build_zorder_reflex_index(
                pts_idx,
                refl,
                count(refl),
                min_x,
                min_y,
                inv_x,
                inv_y,
            )

            function _verify_zidx(zidx)
                seq = Int[]
                cur = zidx.rz_head
                prev_c = 0
                while cur != 0
                    @test zidx.rz_prev[cur] == prev_c
                    @test zidx.in_rz[cur]
                    @test zidx.is_reflex[cur]
                    push!(seq, cur)
                    prev_c = cur
                    cur = zidx.rz_next[cur]
                end
                @test prev_c == zidx.rz_tail
                @test length(seq) == count(zidx.in_rz)
                @test issorted([zidx.z_code[v] for v in seq])

                seen_buckets = Set{Int}()
                for (pos, v) in enumerate(seq)
                    b = Int(zidx.z_code[v] >> zidx.shift) + 1
                    is_first = (pos == 1) || ((Int(zidx.z_code[seq[pos - 1]] >> zidx.shift) + 1) != b)
                    if is_first
                        push!(seen_buckets, b)
                        @test zidx.b_first[b] == v
                    end
                end
                for b in 1:length(zidx.b_first)
                    w = ((b - 1) >>> 6) + 1
                    bit = (zidx.l0[w] >>> ((b - 1) & 63)) & UInt64(1)
                    if b in seen_buckets
                        @test zidx.b_first[b] != 0
                        @test bit == UInt64(1)
                        @test ((zidx.l1 >>> (w - 1)) & UInt64(1)) == UInt64(1)
                    else
                        @test zidx.b_first[b] == 0
                        @test bit == UInt64(0)
                    end
                end
            end

            _verify_zidx(zidx)
            for v in [1, 7, 10, 64, 100, 199]
                BasicCompGeometry._remove_reflex!(zidx, v)
            end
            _verify_zidx(zidx)
            for v in [2, 10, 64, 150, 199]
                BasicCompGeometry._add_reflex!(zidx, v)
            end
            _verify_zidx(zidx)
        end

        @testset "ExactKernel, PntSeq, extreme scales, and fallbacks" begin
            u_poly = [
                point(0 // 1, 0 // 1), point(4 // 1, 0 // 1), point(4 // 1, 4 // 1),
                point(3 // 1, 4 // 1), point(3 // 1, 1 // 1), point(1 // 1, 1 // 1),
                point(1 // 1, 4 // 1), point(0 // 1, 4 // 1),
            ]
            res_u = triangulate_ear_cut_z_order(u_poly; kernel = exact)
            _verify_ear_cut_z_order_result(u_poly, res_u; kernel = exact)
            @test triangulate_polygon(u_poly; kernel = exact, algorithm = :ear_cut_z_order).diagonals == res_u.diagonals
            @test polygon_diagonals(u_poly; kernel = exact, algorithm = :ear_cut_z_order) == res_u.diagonals
            @test triangulate_polygon(PntSeq(u_poly); kernel = exact, algorithm = :ear_cut_z_order).diagonals == res_u.diagonals

            # Collinear chains + shared-x jagged reflex vertices (Rational{BigInt})
            collinear_poly = [
                point(big(0) // 1, big(0) // 1), point(big(1) // 1, big(0) // 1),
                point(big(2) // 1, big(0) // 1), point(big(3) // 1, big(0) // 1),
                point(big(5) // 2, big(1) // 1), point(big(2) // 1, big(2) // 1),
                point(big(3) // 2, big(3) // 1), point(big(1) // 1, big(2) // 1),
                point(big(1) // 2, big(1) // 1),
            ]
            _verify_ear_cut_z_order_result(collinear_poly, triangulate_ear_cut_z_order(collinear_poly; kernel = exact); kernel = exact)

            # Large exact comb (> 16 reflex vertices -> Z-order index active)
            exact_comb = Point{2,Rational{BigInt}}[]
            for i in 1:25
                y0 = Rational{BigInt}(2 * (i - 1))
                push!(exact_comb, point(big(10) // 1, y0))
                push!(exact_comb, point(big(10) // 1, y0 + 1))
                push!(exact_comb, point(big(1) // 1 + (i % 3) // 2, y0 + 6 // 5))
                push!(exact_comb, point(big(1) // 1 + (i % 3) // 2, y0 + 9 // 5))
            end
            pushfirst!(exact_comb, point(big(0) // 1, big(0) // 1))
            push!(exact_comb, point(big(10) // 1, big(50) // 1), point(big(0) // 1, big(50) // 1))
            _verify_ear_cut_z_order_result(exact_comb, triangulate_ear_cut_z_order(exact_comb; kernel = exact); kernel = exact)

            # BigInt > 10^308 non-finite Float64 span fallback
            huge = big(10)^320
            huge_comb = [Point{2,Rational{BigInt}}(p.x * huge, p.y * huge) for p in exact_comb]
            _verify_ear_cut_z_order_result(huge_comb, triangulate_ear_cut_z_order(huge_comb; kernel = exact); kernel = exact)

            # Extreme Float64 coordinate scales (tiny 1e-12 and huge 1e12)
            for s in (1e-12, 1e12)
                comb_s = [point(Float64(p.x) * s, Float64(p.y) * s) for p in exact_comb]
                _verify_ear_cut_z_order_result(comb_s, triangulate_ear_cut_z_order(comb_s))
            end
        end

        @testset "Near-linear practical performance vs ear_cut" begin
            # 1. r = 0 convex polygon stays O(n)
            poly_conv = Point2F[point(cos(2pi * i / 2000), sin(2pi * i / 2000)) for i in 0:1999]
            kc = _EarCutZOrderCountKernel(FloatingKernel(), 0)
            res_conv = triangulate_polygon(poly_conv; kernel = kc, algorithm = :ear_cut_z_order)
            @test length(res_conv.triangles) == length(poly_conv) - 2
            @test kc.orients / length(poly_conv) < 5.0

            # 2. Low-reflex r = 2 polygon (n = 4000)
            poly_lowr = Point2F[
                point((i == 1000 || i == 3000 ? 0.995 : 1.0) * cos(2pi * i / 4000),
                      (i == 1000 || i == 3000 ? 0.995 : 1.0) * sin(2pi * i / 4000))
                for i in 0:3999
            ]
            klr = _EarCutZOrderCountKernel(FloatingKernel(), 0)
            res_lr = triangulate_polygon(poly_lowr; kernel = klr, algorithm = :ear_cut_z_order)
            @test length(res_lr.triangles) == 3998
            @test klr.orients / length(poly_lowr) < 20.0

            # 3. High-reflex non-monotone horizontal comb (n = 1204 and n = 4804)
            function _hcomb(m_prongs::Int)
                pts = Point2F[point(0.0, 0.0)]
                for i in 1:m_prongs
                    yb = 2.0 * (i - 1)
                    xn = 1.2 + 1.4 * (i % 2) + 0.5 * sin(pi * i / (m_prongs + 1))
                    push!(pts, point(10.0, yb), point(10.0, yb + 1.0), point(xn, yb + 1.2), point(xn, yb + 1.8))
                end
                push!(pts, point(10.0, 2.0 * m_prongs), point(10.0, 2.0 * m_prongs + 1.0), point(0.0, 2.0 * m_prongs + 1.0))
                return pts
            end

            comb_1k = _hcomb(300)
            comb_5k = _hcomb(1200)
            kec1 = _EarCutZOrderCountKernel(FloatingKernel(), 0)
            kec5 = _EarCutZOrderCountKernel(FloatingKernel(), 0)
            kcl1 = _EarCutZOrderCountKernel(FloatingKernel(), 0)

            r_ec1 = triangulate_polygon(comb_1k; kernel = kec1, algorithm = :ear_cut_z_order)
            r_ec5 = triangulate_polygon(comb_5k; kernel = kec5, algorithm = :ear_cut_z_order)
            triangulate_polygon(comb_1k; kernel = kcl1, algorithm = :ear_cut)

            _verify_ear_cut_z_order_result(comb_1k, r_ec1)
            @test length(r_ec5.triangles) == length(comb_5k) - 2
            @test isapprox(
                twice_area(comb_5k),
                sum(twice_area([comb_5k[a], comb_5k[b], comb_5k[c]]) for (a, b, c) in r_ec5.triangles);
                rtol = 1e-8,
            )
            @test kec1.orients / length(comb_1k) < 10.0
            @test kec5.orients / length(comb_5k) < 10.0
            @test kec1.orients * 20 < kcl1.orients

            # 4. Random star (n = 1000) and tight spiral (n = 1026)
            rng_s = Xoshiro(77)
            star = [
                let r = iseven(i) ? 10.0 : 1.5 + 1.5 * rand(rng_s)
                    point(r * cos(2pi * i / 1000), r * sin(2pi * i / 1000))
                end
                for i in 0:999
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
            _verify_ear_cut_z_order_result(star, triangulate_polygon(star; algorithm = :ear_cut_z_order))
            _verify_ear_cut_z_order_result(spiral, triangulate_polygon(spiral; algorithm = :ear_cut_z_order))
        end
    end

    mutable struct _EarCutKdTreeCountKernel{K<:AbstractKernel} <: AbstractKernel
        inner::K
        orients::Int
    end
    BasicCompGeometry.orient_sign(k::_EarCutKdTreeCountKernel, a::Point{2}, b::Point{2}, c::Point{2}) =
        (k.orients += 1; orient_sign(k.inner, a, b, c))
    BasicCompGeometry.turn_sign(k::_EarCutKdTreeCountKernel, a::Point{2}, b::Point{2}, c::Point{2}) =
        (k.orients += 1; turn_sign(k.inner, a, b, c))

    function _verify_ear_cut_kd_tree_result(poly, res; kernel = FloatingKernel())
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
                for j in (i + 1):length(diag_segs)
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
            @test all(i -> get(edge_counts, (min(i, mod1(i + 1, m)), max(i, mod1(i + 1, m))), 0) == 1, 1:m)
            @test all(d -> get(edge_counts, d, 0) == 2, res.diagonals)
            if kernel isa ExactKernel
                @test tri_sum == twice_area(poly; kernel = kernel)
            else
                @test isapprox(tri_sum, twice_area(poly; kernel = kernel); rtol = 1e-8)
            end
    end

    for (ts_name, alg_sym, alg_fn) in (
        ("Kd-Tree Ear-Cut Polygon Triangulation (triangulate_ear_cut_kd_tree)", :ear_cut_kd_tree, triangulate_ear_cut_kd_tree),
        ("Unshrunk Kd-Tree Ear-Cut Polygon Triangulation (triangulate_ear_cut_u_kd_tree)", :ear_cut_u_kd_tree, triangulate_ear_cut_u_kd_tree),
    )
        @testset "$ts_name" begin
            @testset "Validation, dispatch, and edge cases" begin
                @test_throws ArgumentError alg_fn(Point2F[point(0.0, 0.0), point(1.0, 0.0)])
                cw_sq = Point2F[point(0.0, 0.0), point(0.0, 1.0), point(1.0, 1.0), point(1.0, 0.0)]
                @test_throws ArgumentError alg_fn(cw_sq)
                @test_throws ArgumentError triangulate_polygon(cw_sq; algorithm = alg_sym)

                tri3 = Point2F[point(0.0, 0.0), point(2.0, 0.0), point(0.5, 1.5)]
                res3 = alg_fn(tri3)
                @test isempty(res3.diagonals)
                @test res3.triangles == [(1, 2, 3)]

                u_poly = Point2F[
                    point(0.0, 0.0), point(4.0, 0.0), point(4.0, 4.0), point(3.0, 4.0),
                    point(3.0, 1.0), point(1.0, 1.0), point(1.0, 4.0), point(0.0, 4.0),
                ]
                res_kd = alg_fn(u_poly)
                @test triangulate_polygon(u_poly; algorithm = alg_sym).diagonals == res_kd.diagonals
                @test triangulate_polygon(PntSeq(u_poly); algorithm = alg_sym).diagonals == res_kd.diagonals
                @test polygon_diagonals(u_poly; algorithm = alg_sym) == res_kd.diagonals
                _verify_ear_cut_kd_tree_result(u_poly, res_kd)
            end

            @testset "ExactKernel, extreme scales, and collinear chains" begin
                exact = ExactKernel()
                collinear_poly = [
                    point(big(0) // 1, big(0) // 1), point(big(1) // 1, big(0) // 1),
                    point(big(2) // 1, big(0) // 1), point(big(3) // 1, big(0) // 1),
                    point(big(5) // 2, big(1) // 1), point(big(2) // 1, big(2) // 1),
                    point(big(3) // 2, big(3) // 1), point(big(1) // 1, big(2) // 1),
                    point(big(1) // 2, big(1) // 1),
                ]
                _verify_ear_cut_kd_tree_result(collinear_poly, alg_fn(collinear_poly; kernel = exact); kernel = exact)

                exact_comb = Point{2,Rational{BigInt}}[]
                for i in 1:25
                    y0 = Rational{BigInt}(2 * (i - 1))
                    push!(exact_comb, point(big(10) // 1, y0))
                    push!(exact_comb, point(big(10) // 1, y0 + 1))
                    push!(exact_comb, point(big(1) // 1 + (i % 3) // 2, y0 + 6 // 5))
                    push!(exact_comb, point(big(1) // 1 + (i % 3) // 2, y0 + 9 // 5))
                end
                pushfirst!(exact_comb, point(big(0) // 1, big(0) // 1))
                push!(exact_comb, point(big(10) // 1, big(50) // 1), point(big(0) // 1, big(50) // 1))
                _verify_ear_cut_kd_tree_result(exact_comb, alg_fn(exact_comb; kernel = exact); kernel = exact)

                huge = big(10)^320
                huge_comb = [Point{2,Rational{BigInt}}(p.x * huge, p.y * huge) for p in exact_comb]
                _verify_ear_cut_kd_tree_result(huge_comb, alg_fn(huge_comb; kernel = exact); kernel = exact)

                if alg_sym == :ear_cut_u_kd_tree
                    int_comb = Point{2,Int}[point(0, 0)]
                    for i in 1:30
                        y0 = 10 * (i - 1)
                        xn = 2 + (i % 4)
                        push!(int_comb, point(20, y0), point(20, y0 + 4), point(xn, y0 + 5), point(xn, y0 + 9))
                    end
                    push!(int_comb, point(20, 300), point(0, 300))
                    _verify_ear_cut_kd_tree_result(int_comb, alg_fn(int_comb; kernel = exact); kernel = exact)
                end

                for s in (1e-12, 1e12)
                    comb_s = [point(Float64(p.x) * s, Float64(p.y) * s) for p in exact_comb]
                    _verify_ear_cut_kd_tree_result(comb_s, alg_fn(comb_s))
                end
            end

            @testset "Subquadratic practical performance vs ear_cut" begin
                poly_conv = Point2F[point(cos(2pi * i / 2000), sin(2pi * i / 2000)) for i in 0:1999]
                kc = _EarCutKdTreeCountKernel(FloatingKernel(), 0)
                res_conv = triangulate_polygon(poly_conv; kernel = kc, algorithm = alg_sym)
                @test length(res_conv.triangles) == length(poly_conv) - 2
                @test kc.orients / length(poly_conv) < 5.0

                poly_lowr = Point2F[
                    point((i == 1000 || i == 3000 ? 0.995 : 1.0) * cos(2pi * i / 4000),
                          (i == 1000 || i == 3000 ? 0.995 : 1.0) * sin(2pi * i / 4000))
                    for i in 0:3999
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
                        push!(pts, point(10.0, yb), point(10.0, yb + 1.0), point(xn, yb + 1.2), point(xn, yb + 1.8))
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
                triangulate_polygon(comb_1k; kernel = kcl1, algorithm = :ear_cut)

                _verify_ear_cut_kd_tree_result(comb_1k, r_ec1)
                @test length(r_ec5.triangles) == length(comb_5k) - 2
                @test isapprox(
                    twice_area(comb_5k),
                    sum(twice_area([comb_5k[a], comb_5k[b], comb_5k[c]]) for (a, b, c) in r_ec5.triangles);
                    rtol = 1e-8,
                )
                @test kec1.orients / length(comb_1k) < 35.0
                @test kec5.orients / length(comb_5k) < 35.0
                @test (kec5.orients / length(comb_5k)) - (kec1.orients / length(comb_1k)) < 8.0
                @test kec1.orients * 20 < kcl1.orients

                rng_s = Xoshiro(77)
                star = [
                    let r = iseven(i) ? 10.0 : 1.5 + 1.5 * rand(rng_s)
                        point(r * cos(2pi * i / 1000), r * sin(2pi * i / 1000))
                    end
                    for i in 0:999
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
                        r = 10.0 + 2.0 * sin(13 * th) + 1.1 * cos(31 * th) + 0.6 * sin(73 * th)
                        point(r * cos(th), r * sin(th))
                    end
                    for i in 0:1999
                ]
                _verify_ear_cut_kd_tree_result(star, triangulate_polygon(star; algorithm = alg_sym))
                _verify_ear_cut_kd_tree_result(spiral, triangulate_polygon(spiral; algorithm = alg_sym))
                _verify_ear_cut_kd_tree_result(coast, triangulate_polygon(coast; algorithm = alg_sym))
            end
        end
    end
end
