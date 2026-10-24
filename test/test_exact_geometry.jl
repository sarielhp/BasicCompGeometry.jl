using Test
using BasicCompGeometry

@testset "Exact geometric primitives" begin
    exact = ExactKernel()

    @testset "Point equality and lexicographic comparison" begin
        p1 = point(big(1) // 3, big(2) // 5)
        p2 = point(big(1) // 3, big(2) // 5)
        p3 = point(big(1) // 3, big(3) // 5)
        p4 = point(big(2) // 5, big(1) // 5)

        @test p1 == p2
        @test p1 != p3
        @test compare_xy(p1, p2; kernel = exact) == 0
        @test compare_xy(p1, p3; kernel = exact) == -1
        @test compare_xy(p3, p1; kernel = exact) == 1
        @test compare_xy(p1, p4; kernel = exact) == -1
        @test compare_xy(p4, p1; kernel = exact) == 1
        @test isless(p1, p3)
        @test !isless(p1, p2)

        # Huge coordinate exact cancellation where floats lose precision
        n = big(2)^120
        pa = point(n, n - 1)
        pb = point(n, n)
        @test compare_xy(pa, pb; kernel = exact) == -1
        @test pa != pb
    end

    @testset "Orientation determinant and sign" begin
        a = point(big(0) // 1, big(0) // 1)
        b = point(big(1) // 1, big(0) // 1)
        c = point(big(0) // 1, big(1) // 1)

        # Left turn (CCW)
        @test orient(a, b, c; kernel = exact) == big(1) // 1
        @test orient_sign(a, b, c; kernel = exact) == 1
        @test (@inferred orient_sign(a, b, c)) == 1

        # Right turn (CW)
        @test orient(a, c, b; kernel = exact) == -big(1) // 1
        @test orient_sign(a, c, b; kernel = exact) == -1

        # Collinear
        collinear_c = point(big(2) // 1, big(0) // 1)
        @test orient(a, b, collinear_c; kernel = exact) == big(0) // 1
        @test orient_sign(a, b, collinear_c; kernel = exact) == 0

        # Determinant value preservation vs sign
        scaled_c = point(big(0) // 1, big(5) // 1)
        @test orient(a, b, scaled_c; kernel = exact) == big(5) // 1
        @test orient_sign(a, b, scaled_c; kernel = exact) == 1

        # Extreme magnitude exact cancellation
        n = big(2)^60
        p = point(big(0), big(0))
        q = point(n, n - 1)
        r = point(n + 1, n)
        @test orient(p, q, r; kernel = exact) == big(1)
        @test orient_sign(p, q, r; kernel = exact) == 1
    end

    @testset "Closed segment membership" begin
        a = point(big(0) // 1, big(0) // 1)
        b = point(big(2) // 1, big(2) // 1)
        seg = Segment(a, b)

        # Endpoints
        @test on_segment(a, a, b; kernel = exact)
        @test on_segment(b, a, b; kernel = exact)
        @test on_segment(a, seg; kernel = exact)
        @test on_segment(b, seg; kernel = exact)

        # Interior exact point
        mid = point(big(1) // 1, big(1) // 1)
        @test on_segment(mid, a, b; kernel = exact)
        @test on_segment(mid, seg; kernel = exact)

        # Collinear outside segment (before a and after b)
        @test !on_segment(point(-big(1) // 1, -big(1) // 1), a, b; kernel = exact)
        @test !on_segment(point(big(3) // 1, big(3) // 1), a, b; kernel = exact)

        # Non-collinear points
        @test !on_segment(point(big(1) // 1, big(0) // 1), a, b; kernel = exact)
        @test !on_segment(point(big(0) // 1, big(1) // 1), a, b; kernel = exact)

        # Vertical segment
        v1 = point(big(3) // 7, big(1) // 2)
        v2 = point(big(3) // 7, big(9) // 2)
        v_seg = Segment(v1, v2)
        @test on_segment(v1, v_seg; kernel = exact)
        @test on_segment(v2, v_seg; kernel = exact)
        @test on_segment(point(big(3) // 7, big(5) // 2), v_seg; kernel = exact)
        @test !on_segment(point(big(3) // 7, big(10) // 2), v_seg; kernel = exact)
        @test !on_segment(point(big(4) // 7, big(5) // 2), v_seg; kernel = exact)

        # Degenerate single-point segment
        pt = point(big(1) // 1, big(1) // 1)
        @test on_segment(pt, pt, pt; kernel = exact)
        @test !on_segment(point(big(1) // 1, big(2) // 1), pt, pt; kernel = exact)
    end

    @testset "Proper segment crossing" begin
        # Classic interior crossing
        s1 = Segment(point(big(0) // 1, big(0) // 1), point(big(2) // 1, big(2) // 1))
        s2 = Segment(point(big(0) // 1, big(2) // 1), point(big(2) // 1, big(0) // 1))
        @test proper_crossing(s1, s2; kernel = exact)
        @test proper_crossing(s1.p, s1.q, s2.p, s2.q; kernel = exact)

        # Vertical segment crossing nonvertical segment
        v_seg = Segment(point(big(1) // 1, -big(1) // 1), point(big(1) // 1, big(1) // 1))
        h_seg = Segment(point(big(0) // 1, big(0) // 1), point(big(2) // 1, big(0) // 1))
        @test proper_crossing(v_seg, h_seg; kernel = exact)

        # Endpoint contact (T-junction) is NOT a proper crossing
        t_seg = Segment(point(big(1) // 1, big(1) // 1), point(big(1) // 1, big(3) // 1))
        @test !proper_crossing(s1, t_seg; kernel = exact)

        # Shared endpoint (corner) is NOT a proper crossing
        corner_seg = Segment(point(big(0) // 1, big(0) // 1), point(big(0) // 1, big(2) // 1))
        @test !proper_crossing(s1, corner_seg; kernel = exact)

        # Collinear overlapping segments are NOT proper crossings
        col1 = Segment(point(big(0) // 1, big(0) // 1), point(big(2) // 1, big(2) // 1))
        col2 = Segment(point(big(1) // 1, big(1) // 1), point(big(3) // 1, big(3) // 1))
        @test !proper_crossing(col1, col2; kernel = exact)

        # Parallel disjoint segments
        par1 = Segment(point(big(0) // 1, big(0) // 1), point(big(2) // 1, big(0) // 1))
        par2 = Segment(point(big(0) // 1, big(1) // 1), point(big(2) // 1, big(1) // 1))
        @test !proper_crossing(par1, par2; kernel = exact)
    end

    @testset "Midpoint and affine interpolation" begin
        p = point(big(1) // 3, big(2) // 7)
        q = point(big(5) // 9, big(4) // 11)
        seg = Segment(p, q)

        # Midpoint
        m = midpoint(p, q)
        @test m == point((big(1) // 3 + big(5) // 9) // 2, (big(2) // 7 + big(4) // 11) // 2)
        @test typeof(m.x) === Rational{BigInt}
        @test midpoint(seg) == m

        # Integer inputs promote to exact Rational
        mi = midpoint(point(1, 2), point(4, 5))
        @test mi == point(5 // 2, 7 // 2)

        # Affine interpolation
        t = big(1) // 4
        interp = affine_interpolate(p, q, t)
        @test interp == p + (q - p) * t
        @test typeof(interp.x) === Rational{BigInt}
        @test affine_interpolate(seg, t) == interp
        @test at(seg, t) == interp
        @test at(seg, big(0) // 1) == p
        @test at(seg, big(1) // 1) == q
    end

    @testset "Nonvertical segment evaluation at x" begin
        a = point(big(1) // 3, big(2) // 5)
        b = point(big(4) // 3, big(7) // 5)
        seg = Segment(a, b)

        # At endpoints
        @test segment_at_x(a, b, big(1) // 3) == a
        @test segment_at_x(a, b, big(4) // 3) == b
        @test segment_at_x(seg, big(1) // 3) == a
        @test segment_at_x(seg, big(4) // 3) == b

        # In interior
        x_mid = big(5) // 6
        res = segment_at_x(seg, x_mid)
        @test res !== nothing
        @test res.x == x_mid
        # Verify collinearity with segment endpoints
        @test orient_sign(a, b, res; kernel = exact) == 0
        @test on_segment(res, seg; kernel = exact)

        # Outside projection returns nothing
        @test segment_at_x(seg, big(0) // 1) === nothing
        @test segment_at_x(seg, big(2) // 1) === nothing

        # Vertical segment violates precondition and must throw ArgumentError
        v1 = point(big(1) // 2, big(0) // 1)
        v2 = point(big(1) // 2, big(3) // 1)
        @test_throws ArgumentError segment_at_x(v1, v2, big(1) // 2)
        @test_throws ArgumentError segment_at_x(Segment(v1, v2), big(1) // 2)
    end

    @testset "Angular comparison of nonzero vectors" begin
        # Precondition enforcement: zero vector throws
        z = point(0, 0)
        u = point(1, 0)
        @test_throws ArgumentError compare_angle(z, u; kernel = exact)
        @test_throws ArgumentError compare_angle(u, z; kernel = exact)
        @test_throws ArgumentError angle_less(z, u; kernel = exact)
        @test_throws ArgumentError angle_less(u, z; kernel = exact)

        # Counterclockwise sequence in [0, 2π)
        # Quadrant 1, boundary, Quadrant 2, boundary, Quadrant 3, boundary, Quadrant 4
        v1 = point(big(1), big(0))        # 0
        v2 = point(big(2), big(1))        # ~26.5 deg
        v3 = point(big(1), big(1))        # 45 deg
        v4 = point(big(0), big(1))        # 90 deg
        v5 = point(big(-1), big(1))       # 135 deg
        v6 = point(big(-1), big(0))       # 180 deg
        v7 = point(big(-1), big(-1))      # 225 deg
        v8 = point(big(0), big(-1))       # 270 deg
        v9 = point(big(1), big(-1))       # 315 deg

        seq = [v1, v2, v3, v4, v5, v6, v7, v8, v9]
        for i in 1:length(seq)-1
            for j in i+1:length(seq)
                @test compare_angle(seq[i], seq[j]; kernel = exact) == -1
                @test compare_angle(seq[j], seq[i]; kernel = exact) == 1
                @test angle_less(seq[i], seq[j]; kernel = exact)
                @test !angle_less(seq[j], seq[i]; kernel = exact)
            end
        end

        # Same-direction equivalence (scalar multiples)
        v_scaled = point(big(3), big(3))
        @test compare_angle(v3, v_scaled; kernel = exact) == 0
        @test compare_angle(v_scaled, v3; kernel = exact) == 0
        @test !angle_less(v3, v_scaled; kernel = exact)
        @test !angle_less(v_scaled, v3; kernel = exact)

        # Angular wraparound: v9 (315 deg) is greater than v1 (0 deg) in [0, 2π) order
        @test !angle_less(v9, v1; kernel = exact)
        @test compare_angle(v9, v1; kernel = exact) == 1
    end

    @testset "Canonical undirected segment key" begin
        a = point(big(1) // 2, big(3) // 4)
        b = point(big(1) // 4, big(5) // 6)

        k1 = segment_key(a, b; kernel = exact)
        k2 = segment_key(b, a; kernel = exact)
        @test k1 == k2
        @test k1.p == b && k1.q == a

        # From Segment
        s = Segment(a, b)
        @test segment_key(s; kernel = exact) == k1

        # Iteration and unpacking
        u, v = k1
        @test u == b && v == a
        @test k1[1] == b && k1[2] == a
        @test length(k1) == 2

        # Tuple version
        @test segment_key((a, b); kernel = exact) == (b, a)
        @test segment_key((b, a); kernel = exact) == (b, a)

        # Usable as dictionary key
        dict = Dict{Segment{2,Rational{BigInt}},Int}()
        dict[segment_key(a, b; kernel = exact)] = 42
        @test dict[segment_key(b, a; kernel = exact)] == 42

        # id_pair
        @test id_pair(5, 2) == (2, 5)
        @test id_pair(1, 4) == (1, 4)
    end

    @testset "Polygon twice_area and segments" begin
        # CCW triangle with area 1/2 -> twice_area == 1
        tri = [point(big(0), big(0)), point(big(1), big(0)), point(big(0), big(1))]
        @test twice_area(tri; kernel = exact) == big(1)

        # CW triangle has negative twice_area
        tri_cw = [point(big(0), big(0)), point(big(0), big(1)), point(big(1), big(0))]
        @test twice_area(tri_cw; kernel = exact) == big(-1)

        # Square of side 3
        sq = [
            point(big(0), big(0)),
            point(big(3), big(0)),
            point(big(3), big(3)),
            point(big(0), big(3)),
        ]
        @test twice_area(sq; kernel = exact) == big(18)

        # Polygon segments
        segs = polygon_segments(tri)
        @test length(segs) == 3
        @test segs[1] == Segment(tri[1], tri[2])
        @test segs[2] == Segment(tri[2], tri[3])
        @test segs[3] == Segment(tri[3], tri[1])
    end
end
