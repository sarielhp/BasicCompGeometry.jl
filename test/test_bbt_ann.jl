using Test
using BasicCompGeometry
using BasicCompGeometry.BBT
using StaticArrays

@testset "Virtual array endpoints and bichromatic WSPD" begin
    values = [10, 20, 30]
    va = BasicCompGeometry.VirtArray.VArray(values)
    BasicCompGeometry.VirtArray.swap!(va, 1, 3)
    @test first(va) == 30
    @test last(va) == 10
    @test collect(va) == [30, 20, 10]
    @test values == [10, 20, 30]
    @test collect(eachindex(va)) == [1, 2, 3]

    points = [
        point(0.0, 0.0),
        point(2.0, 1.0),
        point(4.0, 0.0),
        point(1.0, 3.0),
        point(3.0, 4.0),
        point(5.0, 3.0),
    ]
    state = WSPD.init(PntSeq(points), 0.0)
    top = WSPD.get_top(state)
    @test WSPD.top_diam_ub(state) >= maximum(dist(p, q) for p in points, q in points)
    @test WSPD.min_orig_index(state, top) == 1
    @test WSPD.max_orig_index(state, top) == length(points)
    pairs = WSPD.expand_bichromatic!(state, 3)
    represented = Set{Tuple{Int,Int}}()
    for pair in pairs
        i, j = WSPD.reps_orig_indexes(state, pair)
        @test WSPD.get_reps(state, pair) == (points[i], points[j])
        left, right = WSPD.get_orig_ranges(state, pair)
        @test i in left && j in right
        @test min(i, j) <= 3 < max(i, j)
        push!(represented, minmax(i, j))
    end
    @test represented == Set((i, j) for i = 1:3 for j = 4:6)
    @test WSPD.top_diam_ub(state) == 0.0
    @test WSPD.top_refine!(state) === nothing
    @test BBT.exact_naive_scan(state.tree, points[4]) == (0.0, points[4], 4)
end

@testset "BBT Approximate Nearest Neighbor" begin
    # Create a grid of points
    pnts = [point(Float64(x), Float64(y)) for x = 1:10 for y = 1:10]
    tree = Tree_init(PntSeq(pnts))
    Tree_fully_expand(tree)

    # Query point near (5.1, 5.1)
    q = point(5.1, 5.1)

    # Exact NN (c=1.0)
    d, p, idx = approx_nn(tree, q, 1.0)
    @test p ≈ point(5.0, 5.0)
    @test d ≈ dist(q, point(5.0, 5.0))
    # Original index for (5,5) in a 1:10 grid is (5-1)*10 + 5 = 45
    @test idx == 45

    # Query point far away
    q2 = point(100.0, 100.0)
    d2, p2, idx2 = approx_nn(tree, q2, 1.0)
    @test p2 == point(10.0, 10.0)
    @test idx2 == 100

    # Test with approximation parameter c > 1.0
    # For a grid, c=2.0 should still find something reasonable.
    q3 = point(5.5, 5.5)
    d3, p3, idx3 = approx_nn(tree, q3, 2.0)
    @test d3 >= dist(q3, point(5.0, 5.0)) # At least as good as some point
    # Exact NN for (5.5, 5.5) could be (5,5), (6,5), (5,6), (6,6)
    @test p3 in [point(5.0, 5.0), point(6.0, 5.0), point(5.0, 6.0), point(6.0, 6.0)]

    # Test with many points
    n = 1000
    pnts_rand = [rand_point(3) for _ = 1:n]
    tree_rand = Tree_init(PntSeq(pnts_rand))
    Tree_fully_expand(tree_rand)

    q_rand = rand_point(3)
    d_approx, p_approx, idx_approx = approx_nn(tree_rand, q_rand, 1.5)

    # Brute force search
    best_d = Inf
    best_idx = -1
    for i = 1:n
        cur_d = dist(q_rand, pnts_rand[i])
        if cur_d < best_d
            best_d = cur_d
            best_idx = i
        end
    end

    @test d_approx <= 1.5 * best_d
end

@testset "BBT Split Recording" begin
    pnts = [point(0.0, 0.0), point(10.0, 10.0)]
    tree = Tree_init(PntSeq(pnts))
    Tree_fully_expand(tree)

    root = tree.root
    @test root.split_dim !== nothing
    @test root.split_val == 5.0
    @test root.left !== nothing
    @test root.right !== nothing

    # Check that points are correctly partitioned
    # Pnts are at indices 1 and 2 in the original sequence
    # Partition might have swapped them in tree.PS
    left_pnt = tree.PS[first(root.left.r)]
    right_pnt = tree.PS[first(root.right.r)]

    dim = root.split_dim
    @test left_pnt[dim] <= root.split_val
    @test right_pnt[dim] > root.split_val
end

@testset "BBT Silly NN" begin
    pnts = [point(x, y) for x = 1:10 for y = 1:10]
    tree = Tree_init(PntSeq(pnts))
    Tree_fully_expand(tree)

    q = point(5.1, 5.1)
    d, p, idx = silly_nn(tree, q)

    # Silly NN might not find the exact closest point, 
    # but it should return a point from the tree and its correct distance.
    @test p in pnts
    @test d ≈ dist(q, p)
    @test idx >= 1 && idx <= 100

    # Check that silly_nn distance is >= exact_nn distance
    d_exact, _, _ = approx_nn(tree, q, 1.0)
    @test d >= d_exact
end

@testset "BBT Hybrid NN" begin
    pnts = [point(x, y) for x = 1:10 for y = 1:10]
    tree = Tree_init(PntSeq(pnts))
    Tree_fully_expand(tree)

    q = point(5.1, 5.1)

    # Exact Hybrid (c=1.0)
    d, p, idx = hybrid_nn(tree, q, 1.0)
    @test p ≈ point(5.0, 5.0)
    @test d ≈ dist(q, point(5.0, 5.0))
    @test idx == 45

    # Approx Hybrid (c=1.5)
    d_ap, p_ap, idx_ap = hybrid_nn(tree, q, 1.5)
    @test d_ap <= 1.5 * dist(q, point(5.0, 5.0))
end
