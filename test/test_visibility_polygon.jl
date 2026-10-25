using Test
using BasicCompGeometry

@testset "Visibility Polygon" begin
    # 1. Basic Triangles
    @testset "Convex Triangles" begin
        tri = [point(0.0, 0.0), point(4.0, 0.0), point(1.0, 3.0)]
        # Interior query
        v_int = visibility_polygon(tri, point(1.0, 1.0))
        @test length(v_int) == 3
        @test twice_area(v_int) ≈ twice_area(tri)

        # Boundary edge query
        v_edge = visibility_polygon(tri, point(2.0, 0.0))
        @test length(v_edge) == 3
        @test twice_area(v_edge) ≈ twice_area(tri)

        # Corner vertex query
        v_corn = visibility_polygon(tri, point(0.0, 0.0))
        @test length(v_corn) == 3
        @test twice_area(v_corn) ≈ twice_area(tri)
    end

    # 2. Convex Quadrilaterals & Regular Polygons
    @testset "Convex Polygons" begin
        rect = [point(0.0, 0.0), point(5.0, 0.0), point(5.0, 3.0), point(0.0, 3.0)]
        v_rect = visibility_polygon(rect, point(2.0, 1.0))
        @test length(v_rect) == 4
        @test twice_area(v_rect) ≈ twice_area(rect)

        # Regular octagon
        theta = range(0, 2pi, length=9)[1:8]
        oct = [point(3.0 * cos(t), 3.0 * sin(t)) for t in theta]
        v_oct = visibility_polygon(oct, point(0.2, -0.1))
        @test length(v_oct) == 8
        @test twice_area(v_oct) ≈ twice_area(oct)
    end

    # 3. U-shaped Corridor Polygon
    @testset "U-shaped Polygon" begin
        p_u = [
            point(0.0, 0.0),
            point(4.0, 0.0),
            point(4.0, 4.0),
            point(3.0, 4.0),
            point(3.0, 1.0),
            point(1.0, 1.0),
            point(1.0, 4.0),
            point(0.0, 4.0),
        ]
        # Query from bottom corridor interior
        q_u = point(2.0, 0.5)
        v_u = visibility_polygon(p_u, q_u)
        @test length(v_u) == 6
        @test twice_area(v_u) ≈ 9.0
        @test is_inside(q_u, v_u)
        for p in v_u
            @test is_inside(p, p_u)
        end

        # Query from bottom boundary edge
        v_u_bot = visibility_polygon(p_u, point(2.0, 0.0))
        @test length(v_u_bot) == 6
        @test twice_area(v_u_bot) ≈ 10.0

        # Query from corner vertex
        v_u_corn = visibility_polygon(p_u, point(0.0, 0.0))
        @test length(v_u_corn) == 7
        @test twice_area(v_u_corn) ≈ 14.333333333333332

        # Query from reflex vertex
        v_u_refl = visibility_polygon(p_u, point(1.0, 1.0))
        @test length(v_u_refl) == 6
        @test twice_area(v_u_refl) ≈ 14.0
    end

    # 4. L-shaped Polygon
    @testset "L-shaped Polygon" begin
        p_l = [
            point(0.0, 0.0),
            point(3.0, 0.0),
            point(3.0, 1.0),
            point(1.0, 1.0),
            point(1.0, 3.0),
            point(0.0, 3.0),
        ]
        # From corner (0, 0), the full L-shape is visible
        v_corn = visibility_polygon(p_l, point(0.0, 0.0))
        @test length(v_corn) == 6
        @test twice_area(v_corn) ≈ 10.0

        # From reflex corner (1, 1), the full L-shape is visible
        v_refl = visibility_polygon(p_l, point(1.0, 1.0))
        @test length(v_refl) == 6
        @test twice_area(v_refl) ≈ 10.0

        # From horizontal arm (2.5, 0.5), partial occlusion of vertical arm
        q_arm = point(2.5, 0.5)
        v_arm = visibility_polygon(p_l, q_arm)
        @test length(v_arm) == 5
        @test twice_area(v_arm) ≈ 19.0 / 3.0
    end

    # 5. C-shaped Polygon with reflex occlusion
    @testset "C-shaped Polygon" begin
        p_c = [
            point(0.0, 0.0),
            point(6.0, 0.0),
            point(6.0, 6.0),
            point(0.0, 6.0),
            point(0.0, 4.0),
            point(4.0, 4.0),
            point(4.0, 2.0),
            point(0.0, 2.0),
        ]
        q_c = point(2.0, 1.0)
        v_c = visibility_polygon(p_c, q_c)
        @test length(v_c) == 5
        @test twice_area(v_c) ≈ 26.0
    end

    # 6. Comb / Multi-slot Polygon
    @testset "Comb Polygon" begin
        p_comb = [
            point(0.0, 0.0),
            point(6.0, 0.0),
            point(6.0, 4.0),
            point(4.5, 4.0),
            point(4.5, 1.0),
            point(3.5, 1.0),
            point(3.5, 4.0),
            point(2.5, 4.0),
            point(2.5, 1.0),
            point(1.5, 1.0),
            point(1.5, 4.0),
            point(0.0, 4.0),
        ]
        q_comb = point(3.0, 0.5)
        v_comb = visibility_polygon(p_comb, q_comb)
        @test twice_area(v_comb) ≈ 19.5
        @test twice_area(v_comb) < twice_area(p_comb)
        for p in v_comb
            @test is_inside(p, p_comb)
        end
    end

    # 7. ExactKernel Support
    @testset "ExactKernel" begin
        exact = ExactKernel()

        # Triangle
        tri_ex = [point(0//1, 0//1), point(4//1, 0//1), point(1//1, 3//1)]
        v_tri_ex = visibility_polygon(tri_ex, point(1//1, 1//1); kernel = exact)
        @test eltype(eltype(v_tri_ex)) <: Rational
        @test twice_area(v_tri_ex; kernel = exact) == 12//1

        # U-shape exact
        p_u_ex = [
            point(0//1, 0//1), point(4//1, 0//1), point(4//1, 4//1), point(3//1, 4//1),
            point(3//1, 1//1), point(1//1, 1//1), point(1//1, 4//1), point(0//1, 4//1),
        ]
        v_u_ex = visibility_polygon(p_u_ex, point(2//1, 1//2); kernel = exact)
        @test length(v_u_ex) == 6
        @test twice_area(v_u_ex; kernel = exact) == 9//1
        @test eltype(eltype(v_u_ex)) <: Rational

        # L-shape arm exact (19/3)
        p_l_ex = [
            point(0//1, 0//1), point(3//1, 0//1), point(3//1, 1//1),
            point(1//1, 1//1), point(1//1, 3//1), point(0//1, 3//1),
        ]
        v_l_ex = visibility_polygon(p_l_ex, point(5//2, 1//2); kernel = exact)
        @test length(v_l_ex) == 5
        @test twice_area(v_l_ex; kernel = exact) == 19//3
    end

    # 8. Container and Wrapper types
    @testset "Input Wrappers (PntSeq and coordinate vector)" begin
        tri = [point(0.0, 0.0), point(4.0, 0.0), point(1.0, 3.0)]
        ps = PntSeq(tri)
        v_ps = visibility_polygon(ps, point(1.0, 1.0))
        @test length(v_ps) == 3

        v_vec = visibility_polygon(tri, [1.0, 1.0])
        @test length(v_vec) == 3

        @test is_inside(point(1.0, 1.0), tri)
        @test is_inside(point(1.0, 1.0), ps)
        @test !is_inside(point(10.0, 10.0), tri)
    end

    # 9. Errors & Preconditions
    @testset "Errors and Preconditions" begin
        tri = [point(0.0, 0.0), point(4.0, 0.0), point(1.0, 3.0)]

        # Exterior query point
        @test_throws ArgumentError visibility_polygon(tri, point(10.0, 10.0))

        # Clockwise polygon (negative area)
        cw = [point(0.0, 0.0), point(1.0, 3.0), point(4.0, 0.0)]
        @test_throws ArgumentError visibility_polygon(cw, point(1.0, 1.0))

        # Less than 3 vertices
        @test_throws ArgumentError visibility_polygon([point(0.0, 0.0), point(1.0, 1.0)], point(0.5, 0.5))

        # 3D query point vector
        @test_throws ArgumentError visibility_polygon(tri, [1.0, 1.0, 1.0])
    end

    # 10. Polygon Kernel and Star-Shapedness
    @testset "Polygon Kernel & is_star_shaped" begin
        # Convex polygons (kernel is the polygon itself)
        tri = [point(0.0, 0.0), point(4.0, 0.0), point(1.0, 3.0)]
        k_tri = polygon_kernel(tri)
        @test length(k_tri) == 3
        @test twice_area(k_tri) ≈ twice_area(tri)
        @test is_star_shaped(tri)

        sq = [point(0.0, 0.0), point(2.0, 0.0), point(2.0, 2.0), point(0.0, 2.0)]
        k_sq = polygon_kernel(sq)
        @test length(k_sq) == 4
        @test twice_area(k_sq) ≈ 8.0
        @test is_star_shaped(sq)

        # L-shaped polygon (kernel is corner rectangle [0, 1] x [0, 1])
        p_l = [
            point(0.0, 0.0), point(3.0, 0.0), point(3.0, 1.0),
            point(1.0, 1.0), point(1.0, 3.0), point(0.0, 3.0),
        ]
        k_l = polygon_kernel(p_l)
        @test length(k_l) == 4
        @test twice_area(k_l) ≈ 2.0
        @test is_star_shaped(p_l)

        # A point inside the kernel sees the entire polygon
        q_ker = point(0.5, 0.5)
        @test is_inside(q_ker, k_l)
        v_from_ker = visibility_polygon(p_l, q_ker)
        @test twice_area(v_from_ker) ≈ twice_area(p_l)

        # A point outside the kernel sees strictly less than the full polygon
        q_out = point(2.5, 0.5)
        @test !is_inside(q_out, k_l)
        v_from_out = visibility_polygon(p_l, q_out)
        @test twice_area(v_from_out) < twice_area(p_l)

        # Non-star-shaped polygons (kernel is empty)
        p_u = [
            point(0.0, 0.0), point(4.0, 0.0), point(4.0, 4.0), point(3.0, 4.0),
            point(3.0, 1.0), point(1.0, 1.0), point(1.0, 4.0), point(0.0, 4.0),
        ]
        @test isempty(polygon_kernel(p_u))
        @test !is_star_shaped(p_u)

        p_comb = [
            point(0.0, 0.0), point(6.0, 0.0), point(6.0, 4.0), point(4.5, 4.0),
            point(4.5, 1.0), point(3.5, 1.0), point(3.5, 4.0), point(2.5, 4.0),
            point(2.5, 1.0), point(1.5, 1.0), point(1.5, 4.0), point(0.0, 4.0),
        ]
        @test isempty(polygon_kernel(p_comb))
        @test !is_star_shaped(p_comb)

        # 4-pointed star (kernel is an 8-gon)
        p_star = [
            point(2.0, 0.0), point(0.5, 0.5), point(0.0, 2.0), point(-0.5, 0.5),
            point(-2.0, 0.0), point(-0.5, -0.5), point(0.0, -2.0), point(0.5, -0.5),
        ]
        k_star = polygon_kernel(p_star)
        @test length(k_star) == 8
        @test twice_area(k_star) ≈ 8.0 / 3.0
        @test is_star_shaped(p_star)

        # ExactKernel support
        exact = ExactKernel()
        p_l_ex = [point(Rational{BigInt}(p.x), Rational{BigInt}(p.y)) for p in p_l]
        k_l_ex = polygon_kernel(p_l_ex; kernel = exact)
        @test length(k_l_ex) == 4
        @test twice_area(k_l_ex; kernel = exact) == 2//1
        @test is_star_shaped(p_l_ex; kernel = exact)
        @test eltype(eltype(k_l_ex)) <: Rational

        # PntSeq wrapper
        @test is_star_shaped(PntSeq(p_star))
        @test length(polygon_kernel(PntSeq(p_l))) == 4

        # Error cases
        cw_poly = [point(0.0, 0.0), point(1.0, 3.0), point(4.0, 0.0)]
        @test_throws ArgumentError polygon_kernel(cw_poly)
        @test_throws ArgumentError polygon_kernel([point(0.0, 0.0), point(1.0, 1.0)])
    end
end
