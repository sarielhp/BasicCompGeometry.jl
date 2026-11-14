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
        theta = range(0, 2pi, length = 9)[1:8]
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
            point(0//1, 0//1),
            point(4//1, 0//1),
            point(4//1, 4//1),
            point(3//1, 4//1),
            point(3//1, 1//1),
            point(1//1, 1//1),
            point(1//1, 4//1),
            point(0//1, 4//1),
        ]
        v_u_ex = visibility_polygon(p_u_ex, point(2//1, 1//2); kernel = exact)
        @test length(v_u_ex) == 6
        @test twice_area(v_u_ex; kernel = exact) == 9//1
        @test eltype(eltype(v_u_ex)) <: Rational

        # L-shape arm exact (19/3)
        p_l_ex = [
            point(0//1, 0//1),
            point(3//1, 0//1),
            point(3//1, 1//1),
            point(1//1, 1//1),
            point(1//1, 3//1),
            point(0//1, 3//1),
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
        @test_throws ArgumentError visibility_polygon(
            [point(0.0, 0.0), point(1.0, 1.0)],
            point(0.5, 0.5),
        )

        # 3D query point vector
        @test_throws ArgumentError visibility_polygon(tri, [1.0, 1.0, 1.0])
    end

    # 9b. Vantage points within rounding distance of a vertex or an edge
    @testset "Near-vertex and near-edge vantage points" begin
        near_bd(v, P) =
            minimum(
                BasicCompGeometry.dist_point_segment(v, P[i], P[mod1(i + 1, length(P))]) for
                i in eachindex(P)
            ) <= 1e-9
        inside_p(V, P) = all(is_inside(v, P) || near_bd(v, P) for v in V)

        # Regression: a kernel vertex one ulp off the reflex vertex P[7] used to
        # see past the edges through P[7] (twice-area 2.61 > 2.13 of P).
        P = [
            point(-0.14906678641215848, -0.21386971228658577),
            point(0.1933292193742291, -0.780057608580121),
            point(0.21417178328047168, -0.3277410291682997),
            point(0.834060972879189, -0.2748949833825546),
            point(0.4489924475442709, 0.008136599455512365),
            point(0.5546880821518115, 0.8058826354690696),
            point(0.007679344380308622, 0.3045892312763587),
            point(-0.5526334484435932, 0.6710596152929873),
            point(-0.6096767452069599, 0.359235976593314),
            point(-0.9109573256065367, -0.21459331489368902),
        ]
        A = twice_area(P)
        for q in (point(0.00767934438030855, 0.30458923127635873), P[7])
            V = visibility_polygon(P, q)
            @test twice_area(V) ≈ A rtol = 1e-12
            @test inside_p(V, P)
        end
        for q in polygon_kernel(P)
            @test twice_area(visibility_polygon(P, q)) ≈ A rtol = 1e-9
        end

        # Seeded random star-shaped polygons (with reflex vertices) probed at
        # tiny offsets from every vertex and every edge, on both sides.
        for k in 1:15
            rng = Xoshiro(2024 + k)
            n = rand(rng, 5:14)
            angs = sort(2π .* rand(rng, n))
            Q = [(r = 0.25 + 0.75 * rand(rng); point(r * cos(a), r * sin(a))) for a in angs]
            twice_area(Q) > 0 || continue
            AQ = twice_area(Q)
            for i in 1:n, s in (1e-16, 1e-14, 1e-12)
                a, b = Q[i], Q[mod1(i + 1, n)]
                nrm = point(a.y - b.y, b.x - a.x)
                for q in (
                    a + s * point(randn(rng), randn(rng)),
                    a + rand(rng) * (b - a) + (s * randn(rng)) * nrm,
                )
                    V = visibility_polygon(Q, q)
                    @test 0 < twice_area(V) <= AQ * (1 + 1e-12)
                    @test inside_p(V, Q)
                end
            end
        end
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
            point(0.0, 0.0),
            point(3.0, 0.0),
            point(3.0, 1.0),
            point(1.0, 1.0),
            point(1.0, 3.0),
            point(0.0, 3.0),
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
            point(0.0, 0.0),
            point(4.0, 0.0),
            point(4.0, 4.0),
            point(3.0, 4.0),
            point(3.0, 1.0),
            point(1.0, 1.0),
            point(1.0, 4.0),
            point(0.0, 4.0),
        ]
        @test isempty(polygon_kernel(p_u))
        @test !is_star_shaped(p_u)

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
        @test isempty(polygon_kernel(p_comb))
        @test !is_star_shaped(p_comb)

        # 4-pointed star (kernel is an 8-gon)
        p_star = [
            point(2.0, 0.0),
            point(0.5, 0.5),
            point(0.0, 2.0),
            point(-0.5, 0.5),
            point(-2.0, 0.0),
            point(-0.5, -0.5),
            point(0.0, -2.0),
            point(0.5, -0.5),
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

    # 11. Type Stability, Coordinate Promotion, and Allocation Regression Tests
    @testset "Type stability & coordinate types" begin
        exact = ExactKernel()

        # Integer coordinates under FloatingKernel promote cleanly to Float64 without InexactError
        p_int = [point(0, 0), point(10, 0), point(10, 10), point(5, 3), point(0, 10)]
        v_int = @inferred visibility_polygon(p_int, point(2, 1))
        @test v_int isa Vector{Point{2,Float64}}
        @test length(v_int) >= 3
        @test twice_area(v_int) > 0

        k_int = @inferred polygon_kernel(p_int)
        @test k_int isa Vector{Point{2,Float64}}
        @test twice_area(k_int) > 0

        # Float32 homogeneous return types across _ray_edge_intersect, visibility_polygon, polygon_kernel
        hit_res = @inferred BasicCompGeometry._ray_edge_intersect(
            FloatingKernel(),
            point(0.0f0, 0.0f0),
            point(1.0f0, 0.0f0),
            point(1.0f0, -1.0f0),
            point(1.0f0, 1.0f0),
            Float32,
        )
        @test hit_res isa Tuple{Bool,Float32,Point{2,Float32}}
        miss_res = @inferred BasicCompGeometry._ray_edge_intersect(
            FloatingKernel(),
            point(0.0f0, 0.0f0),
            point(1.0f0, 0.0f0),
            point(-1.0f0, -1.0f0),
            point(-1.0f0, 1.0f0),
            Float32,
        )
        @test miss_res isa Tuple{Bool,Float32,Point{2,Float32}}

        p_u32 = Point{2,Float32}[
            point(0.0f0, 0.0f0),
            point(4.0f0, 0.0f0),
            point(4.0f0, 4.0f0),
            point(3.0f0, 4.0f0),
            point(3.0f0, 1.0f0),
            point(1.0f0, 1.0f0),
            point(1.0f0, 4.0f0),
            point(0.0f0, 4.0f0),
        ]
        v_u32 = @inferred visibility_polygon(p_u32, point(2.0f0, 0.5f0))
        @test v_u32 isa Vector{Point{2,Float32}}
        @test twice_area(v_u32) ≈ 9.0

        p_l32 = Point{2,Float32}[
            point(0.0f0, 0.0f0),
            point(3.0f0, 0.0f0),
            point(3.0f0, 1.0f0),
            point(1.0f0, 1.0f0),
            point(1.0f0, 3.0f0),
            point(0.0f0, 3.0f0),
        ]
        k_l32 = @inferred polygon_kernel(p_l32)
        @test k_l32 isa Vector{Point{2,Float32}}
        @test twice_area(k_l32) ≈ 2.0

        # Rational{Int64} preserved under ExactKernel without Union{Rational{BigInt}, Rational{Int64}}
        p_rat = [
            point(0//1, 0//1),
            point(10//1, 0//1),
            point(10//1, 10//1),
            point(5//1, 3//1),
            point(0//1, 10//1),
        ]
        v_rat = @inferred visibility_polygon(p_rat, point(2//1, 1//1); kernel = exact)
        @test v_rat isa Vector{Point{2,Rational{Int64}}}
        k_rat = @inferred polygon_kernel(p_rat; kernel = exact)
        @test k_rat isa Vector{Point{2,Rational{Int64}}}

        # Allocation regression guard on 40-vertex star polygon (previously >650KB / >380KB due to Core.Box)
        star40 = [
            point(
                (iseven(i) ? 2.0 : 1.0) * cos(2pi * i / 40),
                (iseven(i) ? 2.0 : 1.0) * sin(2pi * i / 40),
            ) for i in 0:39
        ]
        q0 = point(0.0, 0.0)
        @inferred visibility_polygon(star40, q0)
        @inferred polygon_kernel(star40)
        alloc_vis = let p = star40, q = q0
            @allocated visibility_polygon(p, q)
        end
        alloc_ker = let p = star40
            @allocated polygon_kernel(p)
        end
        @test alloc_vis < 32_768
        @test alloc_ker < 32_768
    end
end
