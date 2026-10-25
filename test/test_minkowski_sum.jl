using Test
using BasicCompGeometry
using Random

@testset "Minkowski Sum & is_convex" begin
    @testset "is_convex predicate" begin
        # Standard convex polygons
        sq = [point(0.0, 0.0), point(1.0, 0.0), point(1.0, 1.0), point(0.0, 1.0)]
        @test is_convex(sq)
        @test is_convex(reverse(sq)) # CW orientation is recognized as convex

        tri = [point(0.0, 0.0), point(2.0, 0.0), point(1.0, 2.0)]
        @test is_convex(tri)

        hex = [point(cos(2π * k / 6), sin(2π * k / 6)) for k in 0:5]
        @test is_convex(hex)

        # Collinear intermediate vertices along an edge
        sq_collinear = [point(0.0, 0.0), point(0.5, 0.0), point(1.0, 0.0), point(1.0, 1.0), point(0.0, 1.0)]
        @test is_convex(sq_collinear)

        # Duplicate consecutive vertices
        sq_dup = [point(0.0, 0.0), point(1.0, 0.0), point(1.0, 0.0), point(1.0, 1.0), point(0.0, 1.0)]
        @test is_convex(sq_dup)

        # Non-convex shapes
        l_shape = [
            point(0.0, 0.0), point(3.0, 0.0), point(3.0, 1.0),
            point(1.0, 1.0), point(1.0, 3.0), point(0.0, 3.0),
        ]
        @test !is_convex(l_shape)

        u_shape = [
            point(0.0, 0.0), point(3.0, 0.0), point(3.0, 3.0),
            point(2.0, 3.0), point(2.0, 1.0), point(1.0, 1.0),
            point(1.0, 3.0), point(0.0, 3.0),
        ]
        @test !is_convex(u_shape)

        # Self-intersecting star pentagram (all turns left, but wraps around twice: 4π)
        star5_angles = [2π * (2k % 5) / 5 for k in 0:4]
        star5 = [point(cos(a), sin(a)) for a in star5_angles]
        @test !is_convex(star5)

        # Degenerate shapes
        p_pt = point(1.0, 2.0)
        p_seg = Segment(point(0.0, 0.0), point(1.0, 1.0))
        @test is_convex(p_pt)
        @test is_convex(p_seg)
        @test is_convex([point(0.0, 0.0), point(1.0, 1.0)])

        # ExactKernel
        exact = ExactKernel()
        sq_ex = [point(0//1, 0//1), point(1//1, 0//1), point(1//1, 1//1), point(0//1, 1//1)]
        @test is_convex(sq_ex; kernel = exact)
        l_ex = [point(Rational{BigInt}(p.x), Rational{BigInt}(p.y)) for p in l_shape]
        @test !is_convex(l_ex; kernel = exact)

        # PntSeq wrapper
        @test is_convex(PntSeq(sq))
        @test !is_convex(PntSeq(l_shape))
    end

    @testset "minkowski_sum on convex polygons" begin
        # 1. Square + Square
        sq1 = [point(0.0, 0.0), point(1.0, 0.0), point(1.0, 1.0), point(0.0, 1.0)]
        sq2 = [point(0.0, 0.0), point(2.0, 0.0), point(2.0, 2.0), point(0.0, 2.0)]
        m_sq = minkowski_sum(sq1, sq2)
        @test length(m_sq) == 4
        @test twice_area(m_sq) ≈ 18.0 # [0, 3] x [0, 3] -> area 9, twice_area 18
        @test m_sq[1] ≈ point(0.0, 0.0)
        @test m_sq[3] ≈ point(3.0, 3.0)

        # 2. Square + Triangle
        tri = [point(0.0, 0.0), point(1.0, 0.0), point(0.0, 1.0)]
        m_tri = minkowski_sum(sq1, tri)
        @test length(m_tri) == 5
        # Vertices: (0,0), (2,0), (2,1), (1,2), (0,2)
        # Area = 2*2 - 0.5 = 3.5, twice_area = 7.0
        @test twice_area(m_tri) ≈ 7.0

        # 3. Swept polygon: Square + diagonal segment
        seg = Segment(point(0.0, 0.0), point(1.0, 1.0))
        m_seg = minkowski_sum(sq1, seg)
        @test length(m_seg) == 6 # Hexagon
        @test twice_area(m_seg) ≈ 6.0 # Base square area 1 + swept area 2 -> area 3, twice_area 6.0

        # 4. Segment + Segment (parallelogram)
        seg1 = Segment(point(0.0, 0.0), point(2.0, 0.0))
        seg2 = Segment(point(0.0, 0.0), point(0.0, 3.0))
        m_segs = minkowski_sum(seg1, seg2)
        @test length(m_segs) == 4
        @test twice_area(m_segs) ≈ 12.0 # rectangle 2 x 3

        # 5. Point + Polygon (pure translation)
        pt = point(10.0, 20.0)
        m_pt = minkowski_sum(sq1, pt)
        @test length(m_pt) == 4
        @test m_pt[1] ≈ point(10.0, 20.0)
        @test m_pt[3] ≈ point(11.0, 21.0)
        @test twice_area(m_pt) ≈ twice_area(sq1)

        m_pt_rev = minkowski_sum(pt, sq1)
        @test m_pt_rev[1] ≈ point(10.0, 20.0)

        # 6. Point + Segment
        m_p_seg = minkowski_sum(pt, seg1)
        @test m_p_seg isa Segment
        @test m_p_seg.p ≈ point(10.0, 20.0)
        @test m_p_seg.q ≈ point(12.0, 20.0)

        # 7. Point + Point
        m_pp = minkowski_sum(point(1.0, 2.0), point(3.0, 4.0))
        @test m_pp ≈ point(4.0, 6.0)

        # 8. Clockwise input handled seamlessly
        sq1_cw = reverse(sq1)
        m_cw = minkowski_sum(sq1_cw, sq2)
        @test length(m_cw) == 4
        @test twice_area(m_cw) ≈ 18.0

        # 9. PntSeq wrapper
        m_pntseq = minkowski_sum(PntSeq(sq1), PntSeq(tri))
        @test length(m_pntseq) == 5
        @test twice_area(m_pntseq) ≈ 7.0
    end

    @testset "minkowski_sum with ExactKernel" begin
        exact = ExactKernel()
        sq1_ex = [point(0//1, 0//1), point(1//1, 0//1), point(1//1, 1//1), point(0//1, 1//1)]
        tri_ex = [point(0//1, 0//1), point(1//2, 0//1), point(0//1, 1//2)]

        m_ex = minkowski_sum(sq1_ex, tri_ex; kernel = exact)
        @test length(m_ex) == 5
        @test eltype(eltype(m_ex)) <: Rational
        # Vertices: (0, 0), (3//2, 0), (3//2, 1), (1, 3//2), (0, 3//2)
        # Area = (3/2 * 3/2) - (1/2 * 1/2)/2 = 9/4 - 1/8 = 17/8, twice_area = 17//4
        @test twice_area(m_ex; kernel = exact) == 17//4
    end

    @testset "minkowski_sum convex hull and non-convex error handling" begin
        l_shape = [
            point(0.0, 0.0), point(3.0, 0.0), point(3.0, 1.0),
            point(1.0, 1.0), point(1.0, 3.0), point(0.0, 3.0),
        ]
        sq = [point(0.0, 0.0), point(1.0, 0.0), point(1.0, 1.0), point(0.0, 1.0)]

        # Default throws ArgumentError on non-convex inputs
        @test_throws ArgumentError minkowski_sum(l_shape, sq)
        @test_throws ArgumentError minkowski_sum(sq, l_shape)

        # With hull = true, computes CH(P) ⊕ CH(Q)
        m_hull = minkowski_sum(l_shape, sq; hull = true)
        @test length(m_hull) == 5
        @test twice_area(m_hull) ≈ 28.0
        @test is_convex(m_hull)
    end

    @testset "minkowski_sum_decomposition" begin
        l_shape = [
            point(0.0, 0.0), point(3.0, 0.0), point(3.0, 1.0),
            point(1.0, 1.0), point(1.0, 3.0), point(0.0, 3.0),
        ]
        sq = [point(0.0, 0.0), point(1.0, 0.0), point(1.0, 1.0), point(0.0, 1.0)]

        # L-shape triangulates into 4 triangles; each summed with convex square gives 4 convex pieces
        decomp = minkowski_sum_decomposition(l_shape, sq)
        @test length(decomp) == 4
        for piece in decomp
            @test is_convex(piece)
            @test twice_area(piece) > 0
        end

        # Two convex polygons return a single piece
        decomp_conv = minkowski_sum_decomposition(sq, sq)
        @test length(decomp_conv) == 1
        @test length(decomp_conv[1]) == 4
        @test twice_area(decomp_conv[1]) ≈ 8.0 # [0, 2] x [0, 2]

        # ExactKernel support
        exact = ExactKernel()
        l_ex = [point(Rational{BigInt}(p.x), Rational{BigInt}(p.y)) for p in l_shape]
        sq_ex = [point(Rational{BigInt}(p.x), Rational{BigInt}(p.y)) for p in sq]
        decomp_ex = minkowski_sum_decomposition(l_ex, sq_ex; kernel = exact)
        @test length(decomp_ex) == 4
        @test eltype(eltype(decomp_ex[1])) <: Rational
    end

    @testset "Randomized validation against convex_hull of pairwise sums" begin
        rng = Xoshiro(12345)
        for _ in 1:100
            n = rand(rng, 3:15)
            m = rand(rng, 3:15)
            pts1 = [point(rand(rng) * 50.0, rand(rng) * 50.0) for _ in 1:n]
            pts2 = [point(rand(rng) * 50.0, rand(rng) * 50.0) for _ in 1:m]

            poly1 = collect(convex_hull(pts1))
            poly2 = collect(convex_hull(pts2))

            res_fast = minkowski_sum(poly1, poly2)

            all_sums = [p + q for p in poly1 for q in poly2]
            res_hull = collect(convex_hull(all_sums))

            a_fast = twice_area(res_fast)
            a_hull = twice_area(res_hull)

            @test abs(a_fast - a_hull) / max(abs(a_hull), 1.0) < 1e-9
            @test is_convex(res_fast)
        end
    end
end
