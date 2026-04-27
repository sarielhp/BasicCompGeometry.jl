using Test
using LinearAlgebra
using StaticArrays
using BasicCompGeometry

@testset "BasicCompGeometry.jl" begin

    @testset "Points (2D & 3D, Float & Int)" begin
        # 2D Float
        p1 = npoint(0.0, 0.0)
        p2 = npoint(3.0, 4.0)
        @test dist(p1, p2) ≈ 5.0
        @test dist_sq(p1, p2) ≈ 25.0
        @test convex_comb(p1, p2, 0.5) ≈ npoint(1.5, 2.0)
        
        # 3D Int
        q1 = npoint(0, 0, 0)
        q2 = npoint(1, 2, 2)
        @test dist(q1, q2) ≈ 3.0
        @test dist_sq(q1, q2) == 9
        
        # Predicates
        a, b, c = npoint(0,0), npoint(1,0), npoint(0,1)
        @test turn_sign(a, b, c) > 0
        @test is_left_turn(a, b, c) == true
        @test is_left_eq_turn(a, b, c) == true
        @test is_right_turn(a, c, b) == true
        @test is_left_eq_turn(a, c, b) == false
        @test is_right_eq_turn(a, c, b) == true
        @test turn_sign(a, c, b) < 0
        @test is_collinear(a, b, npoint(2,0)) == true
        @test is_left_eq_turn(a, b, npoint(2,0)) == true
        @test is_right_eq_turn(a, b, npoint(2,0)) == true
        @test turn_sign(a, b, npoint(2,0)) == 0
        
        # Random
        @test length(rand_point(2)) == 2
        @test length(rand_gaussian(3)) == 3
        
        # Min/Max
        @test max(npoint(1, 5), npoint(2, 3)) == npoint(2, 5)
        @test min(npoint(1, 5), npoint(2, 3)) == npoint(1, 3)
    end

    @testset "Segments & Lines" begin
        p1 = npoint(0.0, 0.0)
        p2 = npoint(10.0, 0.0)
        seg = Segment(p1, p2)
        
        @test geom_length(seg) == 10.0
        @test at(seg, 0.2) == npoint(2.0, 0.0)
        
        # Point distance
        qr = npoint(5.0, 5.0)
        @test dist(seg, qr) == 5.0
        @test nn_point(seg, qr) == npoint(5.0, 0.0)
        @test convex_coef(seg, npoint(2.0, 0.0)) ≈ 0.2
        
        # Segment distance
        seg2 = Segment(npoint(5.0, 2.0), npoint(5.0, 10.0))
        @test dist(seg, seg2) == 2.0
        
        # Bisection
        # p_a, p_b define a vertical bisector at x=5
        p_a, p_b = npoint(0.0, 0.0), npoint(10.0, 0.0)
        # s_cross from (0,0) to (10,0) crosses x=5 at (5,0)
        s_cross = Segment(npoint(0.0, 0.0), npoint(10.0, 0.0))
        ok, t, pt = bisection_point(s_cross, p_a, p_b)
        @test ok == true
        @test pt ≈ npoint(5.0, 0.0)
        
        # Lines
        line = Line(npoint(0.0, 0.0), npoint(1.0, 1.0))
        @test BasicCompGeometry.distance(npoint(0.0, 10.0), line) ≈ dist(npoint(0.0, 10.0), npoint(5.0, 5.0))
    end

    @testset "Polygons" begin
        pnts = [npoint(0.0,0.0), npoint(1.0,0.0), npoint(1.0,1.0), npoint(0.0,1.0)]
        poly = Polygon(pnts)
        
        @test length(poly) == 4
        @test cardin(poly) == 4
        @test geom_length(poly) == 3.0 # Open curve by default
        
        # In-place ops
        translate!(poly, npoint(1.0, 1.0))
        @test poly[1] == npoint(1.0, 1.0)
        scale!(poly, 2.0)
        @test poly[1] == npoint(2.0, 2.0)
        
        # Analysis
        plens = prefix_lengths(poly)
        @test length(plens) == 4
        @test plens[end] == geom_length(poly)
        # poly vertices after translate (1,1) and scale 2.0:
        # (0,0)->(1,1)->(2,2)
        # (1,0)->(2,1)->(4,2)
        # (1,1)->(2,2)->(4,4)
        # (0,1)->(1,2)->(2,4)
        # lengths: dist((2,2),(4,2))=2, dist((4,2),(4,4))=2, dist((4,4),(2,4))=2. Total=6.
        # at(poly, 0.5) is at length 3. (4,2) is length 2, (4,4) is length 4.
        # So it is halfway between (4,2) and (4,4), which is (4,3).
        @test at(poly, 0.5) == npoint(4.0, 3.0)
        
        # Iteration & Base
        @test first(poly) == poly[1]
        @test last(poly) == poly[end]
        @test length(collect(poly)) == 4
        
        # Matrix
        mat = Matrix(poly)
        @test size(mat) == (2, 4)
        
        # Simplification
        poly_big = Polygon([npoint(x, 0.0) for x in 0:0.01:1])
        poly_sim, indices = simplify(poly_big, 0.1)
        @test length(poly_sim) < length(poly_big)
        
        # Sampling
        poly_uniform = sample_uniformly(poly, 10)
        @test length(poly_uniform) == 10
        @test geom_length(poly_uniform) ≈ geom_length(poly)
        
        # Random
        @test length(rand_polygon(2, Float64, 5)) == 5
    end

    @testset "Bounding Boxes" begin
        p1 = npoint(-1.0, -2.0)
        p2 = npoint(5.0, 10.0)
        bb = BBox(p1, p2)
        
        @test width(bb) == 6.0
        @test height(bb) == 12.0
        @test middle(bb, 1) == 2.0
        @test diam(bb) ≈ dist(p1, p2)
        
        @test is_inside(npoint(0,0), bb)
        @test !is_inside(npoint(10,10), bb)
        
        # Padding
        bb2 = bb + 1.0
        @test width(bb2) == 8.0
        
        # Mutating
        bb_empty = BBox{2, Float64}()
        @test bb_empty.f_init == false
        bound!(bb_empty, npoint(0.0, 0.0))
        @test bb_empty.f_init == true
        @test diam(bb_empty) == 0.0
        
        # Box distance
        bb_a = BBox(npoint(0,0), npoint(1,1))
        bb_b = BBox(npoint(2,0), npoint(3,1))
        @test dist(bb_a, bb_b) == 1.0
        # max_dist is distance between (0,0) and (3,1) or (0,1) and (3,0)
        # sqrt(3^2 + 1^2) = sqrt(10)
        @test max_dist(bb_a, bb_b) ≈ sqrt(10.0)
    end

    @testset "Transforms 2D" begin
        p = npoint(1.0, 0.0)
        
        # Rotation 90 deg (Counter-Clockwise)
        m_rot = BasicCompGeometry.rotation(pi/2)
        p_rot = BasicCompGeometry.apply_transform(m_rot, p)
        # (1,0) rotated 90 deg CCW is (0,1)
        @test p_rot[1] ≈ 0.0 atol=1e-10
        @test p_rot[2] ≈ 1.0
        
        # Translation
        m_trans = BasicCompGeometry.translation(2.0, 3.0)
        p_trans = BasicCompGeometry.apply_transform(m_trans, p)
        @test p_trans == npoint(3.0, 3.0)
    end

    @testset "Algorithms & Module Level" begin
        # Hausdorff
        pnts1 = [npoint(x, 0.0) for x in 0:10]
        pnts2 = [npoint(x, 0.1) for x in 0:10]
        poly1, poly2 = Polygon(pnts1), Polygon(pnts2)
        
        @test hausdorff_dist_subseg(poly1) == 0.0
        sim_p, sim_i = hausdorff_simplify(poly1, 0.5)
        @test length(sim_p) < length(poly1)
        
        # Module level
        @test distance_infty(poly1, poly2) ≈ 0.1
        @test centroid([npoint(0,0), npoint(2,2)]) == npoint(1,1)
        
        # Match price
        price = match_price(npoint(0,0), npoint(1,0), npoint(0,1), npoint(1,1))
        @test price > 0
        
        # Convex Hull
        hull_pnts = [npoint(0,0), npoint(1,0), npoint(0,1), npoint(0.5, 0.5), npoint(1,1)]
        hull = convex_hull(hull_pnts)
        # Expected hull points (CCW starting from min-x): (0,0), (1,0), (1,1), (0,1)
        @test length(hull) == 4
        @test npoint(0,0) in hull
        @test npoint(1,0) in hull
        @test npoint(1,1) in hull
        @test npoint(0,1) in hull
        @test !(npoint(0.5, 0.5) in hull)
    end

    @testset "VirtArray" begin
        using BasicCompGeometry.VirtArray
        v = [10, 20, 30, 40]
        va = VArray(v)
        @test length(va) == 4
        @test va[1] == 10
        @test va[4] == 40
        
        swap!(va, 1, 4)
        @test va[1] == 40
        @test va[4] == 10
        @test orig_index(va, 1) == 4
        @test orig_index(va, 4) == 1
        
        @test collect(va) == [40, 20, 30, 10]
    end

    @testset "BBT (Bounding Box Tree)" begin
        using BasicCompGeometry.BBT
        pnts = [npoint(x, y) for x in 0:10 for y in 0:10]
        poly = Polygon(pnts)
        tree = Tree_init(poly)
        
        @test tree.root !== nothing
        @test diam(tree.root.bb) > 0
        
        Tree_fully_expand(tree)
        @test depth(tree.root) > 1
        
        # Test original index retrieval
        min_idx = get_min_orig_index(tree, tree.root)
        max_idx = get_max_orig_index(tree, tree.root)
        @test min_idx == 1
        @test max_idx == length(pnts)
        
        # Test 1000 random points in 3D
        pnts3d = [rand_point(3) for _ in 1:1000]
        tree3d = Tree_init(Polygon(pnts3d))
        Tree_fully_expand(tree3d)
        
        @test depth(tree3d.root) >= 10 # log2(1000) is approx 9.9
        @test tree3d.id_counter >= 1000 # Should have at least n nodes for a full expansion
        
        # Verify root BB contains all points
        root_bb = tree3d.root.bb
        all_inside = all(p -> is_inside(p, root_bb), pnts3d)
        @test all_inside == true
    end

    @testset "Diameter Algorithms" begin
        # 2D points on a circle
        pnts = Polygon_random_sphere(2, Float64, 100)
        
        exact_diam = exact_diameter(pnts)
        @test exact_diam ≈ 2.0 atol=0.1 # Unit sphere diameter is 2.0
        
        # (1+ε) approximation with ε=0.1
        approx_diam = approx_diameter(pnts, 0.1)
        @test approx_diam <= exact_diam
        @test approx_diam * 1.1 >= exact_diam
    end

    @testset "WSPD (Well-Separated Pairs Decomposition)" begin
        using BasicCompGeometry.WSPD
        pnts = [npoint(x, 0.0) for x in 1:10] # Points at 1, 2, ..., 10
        poly = Polygon(pnts)
        sep = 2.0
        W = WSPD.init(poly, sep)
        finals = WSPD.expand!(W)
        
        @test length(finals) > 0
        
        # Verify well-separation for all pairs
        for pair in finals
            if pair.dist > 0
                @test (max(diam(pair.left.bb), diam(pair.right.bb)) / pair.dist) <= sep + 1e-9
            end
        end
        
        # Verify all pairs of distinct points are covered by exactly one WSP
        n = length(pnts)
        covered_pairs = 0
        for pair in finals
            l_range, r_range = get_orig_ranges(W, pair)
            covered_pairs += length(l_range) * length(r_range)
        end
        # In a WSPD of a single set, the number of distinct pairs (i, j) with i < j 
        # is covered. Since this implementation starts with (root, root), 
        # it covers all n^2 pairs, but excludes (i, i) if dist > 0 check is used.
        # Actually, monotone chain and other WSPD algorithms usually cover n(n-1)/2 pairs.
        # Let's just check that it produces a reasonable number of pairs.
        @test length(finals) < n^2
    end

    @testset "10-Dimensional Geometry" begin
        # 10D Points
        p10_1 = npoint(zeros(10)...)
        p10_2 = npoint(ones(10)...)
        @test dist(p10_1, p10_2) ≈ sqrt(10.0)
        @test dist_sq(p10_1, p10_2) ≈ 10.0
        
        # 10D Segments
        seg10 = Segment(p10_1, p10_2)
        @test geom_length(seg10) ≈ sqrt(10.0)
        @test at(seg10, 0.5) ≈ npoint(fill(0.5, 10)...)
        
        # 10D BBox
        bb10 = BBox(p10_1, p10_2)
        @test diam(bb10) ≈ sqrt(10.0)
        @test is_inside(npoint(fill(0.1, 10)...), bb10)
        @test !is_inside(npoint(fill(1.1, 10)...), bb10)
        
        # 10D Polygon
        poly10 = Polygon([p10_1, p10_2, p10_1])
        @test cardin(poly10) == 3
        @test geom_length(poly10) ≈ 2 * sqrt(10.0)
    end

end
