using BasicCompGeometry
import BasicCompGeometry: Points
using BasicCompGeometry.BBT
using BasicCompGeometry.WSPD
using Test
using LinearAlgebra

@testset "MatPntSeq" begin
    D, N = 3, 10
    M = rand(D, N)
    mp = MatPntSeq(M)

    @testset "Basics" begin
        @test mp isa AbsPntSeq{3,Float64}
        @test dim(mp) == 3
        @test length(mp) == N
        @test size(mp) == N
        @test mp[1] == M[:, 1]
        @test mp[N] == M[:, N]
        @test Points(mp) === mp.pnts
    end

    @testset "Zero-copy" begin
        original_val = M[1, 1]
        M[1, 1] = 999.0
        @test mp[1][1] == 999.0

        # Test that updating the pnt seq (if possible) updates the matrix
        # Since Point is SVector (immutable), we have to set the whole point
        # But setindex! is defined for AbsPntSeq
        new_p = Point{3,Float64}(1.0, 2.0, 3.0)
        mp[2] = new_p
        @test M[:, 2] == [1.0, 2.0, 3.0]
    end

    @testset "Sharing: translate!, scale!, views, adjoints, fixed size" begin
        A = [0.0 1.0 2.0; 0.0 0.0 1.0]
        P = MatPntSeq(A)
        BasicCompGeometry.translate!(P, Point(1.0, 1.0))
        @test A == [1.0 2.0 3.0; 1.0 1.0 2.0]
        BasicCompGeometry.scale!(P, 2.0)
        @test A == [2.0 4.0 6.0; 2.0 2.0 4.0]

        B = rand(3, 6)
        for S in (view(B, :, 2:4), view(B, 1:2, :), adjoint(rand(5, 3)))
            Q = MatPntSeq(S)
            @test all(Q[j] == S[:, j] for j in axes(S, 2))
            S[1, 1] += 1.0
            @test Q[1][1] == S[1, 1]
            Q[2] = Q[1]
            @test S[:, 2] == S[:, 1]
        end

        p = P[1]
        @test_throws Exception (p[1] = 0.0)
        @test_throws MethodError push!(P, P[1])
        @test_throws MethodError pop!(P)
        @test_throws MethodError empty!(P)
    end

    @testset "copy = true" begin
        A = [0.0 1.0 2.0; 0.0 0.0 1.0]
        before = Base.copy(A)
        P = MatPntSeq(A; copy = true)
        @test P.mat !== A
        @test P.mat isa Matrix{Float64}
        @test all(P[j] == A[:, j] for j in 1:3)
        A[1, 1] = 99.0
        @test P[1][1] == 0.0
        A[1, 1] = 0.0
        P[2] = Point(7.0, 7.0)
        BasicCompGeometry.translate!(P, Point(1.0, 1.0))
        BasicCompGeometry.scale!(P, 3.0)
        @test A == before
        @test Matrix(P) == [3.0 24.0 9.0; 3.0 24.0 6.0]

        V = view(rand(4, 5), 2:3, :)
        Q = MatPntSeq(V; copy = true)
        @test Q.mat isa Matrix{Float64}
        @test all(Q[j] == V[:, j] for j in axes(V, 2))
        @test MatPntSeq(A; copy = false).mat === A
    end

    @testset "Algorithms compatibility" begin
        @test exact_diameter(mp) > 0
        @test approx_diameter(mp, 0.1) > 0

        tree = Tree_init(mp)
        @test depth(tree.root) >= 1

        W = WSPD.init(mp, 2.0)
        finals = WSPD.expand!(W)
        @test length(finals) > 0
    end
end

@testset "MatVecPntSeq" begin
    A = [0.0 1.0 2.0; 0.0 0.0 1.0; 5.0 5.0 5.0]
    P = MatVecPntSeq(A)

    @testset "Basics and type stability" begin
        @test P isa AbsPntSeq{AnyDim,Float64}
        @test dim(P) == 3
        @test length(P) == 3
        @test size(P) == 3
        @test (@inferred MatVecPntSeq(A)) isa MatVecPntSeq{Float64,Matrix{Float64}}
        @test isconcretetype(eltype(P))
        @test eltype(P) === typeof(P[1])
        @test all(P[j] == A[:, j] for j in 1:3)
        @test [collect(p) for p in P] == [A[:, j] for j in 1:3]
        @test first(P) == A[:, 1] && last(P) == A[:, 3]
        @test collect(P[2]) isa Vector{Float64}
        @test P[3] - P[1] == [2.0, 1.0, 0.0]
        @test dist(P, 1, 3) ≈ sqrt(5.0)
        @test dist(P, 2, 2) == 0.0
        @test Matrix(P) == A && Matrix(P) !== A
        @test occursin("3 points in dimension 3", sprint(show, P))
        @test_throws BoundsError P[4]
        @test_throws BoundsError P[1][4]
    end

    @testset "Same dimension-free type for every D" begin
        for D in (1, 2, 7, 500)
            Q = MatVecPntSeq(rand(D, 4))
            @test typeof(Q) === MatVecPntSeq{Float64,Matrix{Float64}}
            @test dim(Q) == D
        end
    end

    @testset "Sharing, read-only points, fixed size" begin
        B = Base.copy(A)
        Q = MatVecPntSeq(B)
        p = Q[2]
        B[1, 2] = 42.0
        @test p[1] == 42.0 && Q[2][1] == 42.0
        @test_throws Exception (p[1] = 0.0)
        Q[3] = [7.0, 8.0, 9.0]
        @test B[:, 3] == [7.0, 8.0, 9.0]
        BasicCompGeometry.translate!(Q, [1.0, 1.0, 1.0])
        @test B[:, 1] == [1.0, 1.0, 6.0]
        BasicCompGeometry.scale!(Q, 2.0)
        @test B[:, 1] == [2.0, 2.0, 12.0]
        @test_throws DimensionMismatch BasicCompGeometry.translate!(Q, [1.0, 1.0])
        @test_throws MethodError push!(Q, [1.0, 2.0, 3.0])
        @test_throws MethodError pop!(Q)
        @test_throws MethodError empty!(Q)

        V = view(rand(4, 5), 2:3, :)
        R = MatVecPntSeq(V)
        V[1, 1] = -1.0
        @test R[1][1] == -1.0
    end

    @testset "copy = true" begin
        B = Base.copy(A)
        Q = MatVecPntSeq(B; copy = true)
        @test Q.mat !== B && Q.mat isa Matrix{Float64}
        B[1, 1] = 99.0
        @test Q[1][1] == 0.0
        Q[2] = [7.0, 7.0, 7.0]
        BasicCompGeometry.translate!(Q, [1.0, 1.0, 1.0])
        @test B == [99.0 1.0 2.0; 0.0 0.0 1.0; 5.0 5.0 5.0]
        R = MatVecPntSeq(view(rand(4, 5), 2:3, :); copy = true)
        @test R.mat isa Matrix{Float64}
        @test MatVecPntSeq(B; copy = false).mat === B
    end

    @testset "Built from a matrix or from the matrix of points: the same" begin
        pts = [[0.0, 1.0, 2.0], [3.0, 5.0, 8.0], [1.0, 1.0, 1.0], [2.0, 0.0, 4.0]]
        from_points = MatVecPntSeq(reduce(hcat, pts))
        from_matrix = MatVecPntSeq([0.0 3.0 1.0 2.0; 1.0 5.0 1.0 0.0; 2.0 8.0 1.0 4.0])
        @test typeof(from_points) === typeof(from_matrix)
        @test [collect(p) for p in from_points] == [collect(p) for p in from_matrix]
        @test all(dist(from_points, i, j) == dist(from_matrix, i, j) for i in 1:4, j in 1:4)
        @test Matrix(from_points) == Matrix(from_matrix)
    end

    @testset "ColumnPoint checks its column" begin
        X = rand(3, 4)
        @test BasicCompGeometry.ColumnPoint(X, 2) == X[:, 2]
        @test_throws BoundsError BasicCompGeometry.ColumnPoint(X, 5)
        @test_throws BoundsError BasicCompGeometry.ColumnPoint(X, 0)
    end

    @testset "Agrees with MatPntSeq" begin
        X = rand(5, 30)
        S, Vq = MatPntSeq(X), MatVecPntSeq(X)
        @test all(collect(Vq[j]) == collect(S[j]) for j in 1:30)
        @test all(
            isapprox(dist(Vq, i, j), dist(S, i, j); rtol = 1.0e-14) for i in 1:30, j in 1:30
        )
    end
end

@testset "GPntSeq with Vector points: length checks and dim" begin
    G = GPntSeq{3,Float64,Vector{Float64}}([[0.0, 0.0, 0.0]])
    push!(G, [1.0, 2.0, 3.0])
    @test length(G) == 2 && dim(G) == 3
    @test_throws DimensionMismatch push!(G, [1.0, 2.0])
    @test_throws DimensionMismatch GPntSeq{3,Float64,Vector{Float64}}([
        [0.0, 0.0, 0.0],
        [1.0, 2.0],
    ])
    A = GPntSeq{AnyDim,Float64,Vector{Float64}}([[1.0, 2.0], [3.0, 4.0]])
    @test dim(A) == 2
    @test_throws DimensionMismatch push!(A, [1.0])
    @test_throws DimensionMismatch GPntSeq{AnyDim,Float64,Vector{Float64}}([
        [1.0],
        [1.0, 2.0],
    ])
    E = GPntSeq{AnyDim,Float64,Vector{Float64}}(Vector{Float64}[])
    @test_throws ArgumentError dim(E)
    push!(E, [1.0, 2.0, 3.0])
    @test dim(E) == 3
    S = PntSeq([point(1.0, 2.0)])
    push!(S, [3.0, 4.0])
    @test S[2] == point(3.0, 4.0)
    @test_throws DimensionMismatch push!(S, [1.0, 2.0, 3.0])
end

@testset "Sequence collection, transformations, and simplification" begin
    points = [point(0.0, 0.0), point(1.0, 1.0), point(2.0, 0.0)]
    sequence = PntSeq(copy(points))
    @test first(sequence) == points[1]
    @test last(sequence) == points[3]
    @test size(sequence, 1) == 3
    @test size(sequence, 2) == 1
    @test dist(sequence, 2, 2) == 0.0
    @test pop!(sequence) == points[3]
    empty!(sequence)
    @test isempty(sequence)
    @test occursin("0 points", sprint(show, sequence))
    matrix_sequence = MatPntSeq(Matrix(PntSeq(points)))
    @test occursin("3 points", sprint(show, matrix_sequence))
    @test apply_transform(translation(point(2.0, -1.0)), points[2]) ≈ point(3.0, 0.0)
    @test Line(point(0, 0), point(1.0, 2.0)) isa Line{2,Float64}
    @test geom_length(Segment{2,Float64}()) == 0.0
    @test occursin("..", sprint(show, Segment(points[1], points[2])))
    @test isapprox(BBox2F(), BBox2F())
    @test BBox2F() == BBox2F()
    @test !isapprox(BBox2F(), BBox(points[1], points[3]))
    box = BBox(points)
    @test bottom_left(box) == point(0.0, 0.0)
    @test top_right(box) == point(2.0, 1.0)
    @test BasicCompGeometry.d_min(box, 1) == 0.0
    @test BasicCompGeometry.d_max(box, 2) == 1.0
    @test occursin("uninitialized", sprint(show, BBox2F()))
    @test union_bbox([box, bbox(Segment(points[1], points[3]))]) == box
    let large_boxes =
            [BBox(point(Float64(i), 0.0), point(Float64(i) + 1.0, 1.0)) for i in 1:10_000]
        expected_large = BBox(point(1.0, 0.0), point(10_001.0, 1.0))
        @test (@inferred union_bbox(large_boxes)) == expected_large
        @test (@allocated union_bbox(large_boxes)) <= 1024
    end
    let mixed_boxes = Any[
            BBox(point(0, 0), point(1, 1)),
            BBox{2,Float64}(),
            Segment(point(0.0, 0.0), point(3.5, 2.5)),
        ]
        ub = union_bbox(mixed_boxes)
        @test ub isa BBox2F
        @test ub == BBox(point(0.0, 0.0), point(3.5, 2.5))
    end
    @test_throws ArgumentError union_bbox(BBox2F[])
    @test_throws DimensionMismatch union_bbox([
        BBox(point(0, 0), point(1, 1)),
        BBox(point(0, 0, 0), point(1, 1, 1)),
    ])

    mktempdir(mkpath(joinpath(@__DIR__, "..", "output", "tmp"))) do dir
        file = joinpath(dir, "points.plt")
        BasicCompGeometry.write_plt(file, PntSeq(points))
        @test length(readlines(file)) == 3
        @test split(first(readlines(file)), ',') == ["0.0", "0.0"]
    end

    @test BasicCompGeometry.hausdorff_dist_subseg(PntSeq(points[1:2])) == 0.0
    empty_curve, empty_indexes = hausdorff_simplify(PntSeq2F(), 0.1)
    @test isempty(empty_curve) && isempty(empty_indexes)
    singleton, indexes = hausdorff_simplify(PntSeq(points[1:1]), 0.1)
    @test collect(singleton) == points[1:1]
    @test indexes == [1]
    zigzag = PntSeq([point(Float64(i), isodd(i) ? 1.0 : 0.0) for i in 0:30])
    for tolerance in (0.1, 0.6, 1.0)
        simplified, indexes = hausdorff_simplify(zigzag, tolerance)
        @test first(indexes) == 1 && last(indexes) == length(zigzag)
        @test collect(simplified) == [zigzag[i] for i in indexes]
        @test all(
            r -> BasicCompGeometry.hausdorff_dist_subseg(zigzag, r) <= tolerance,
            (indexes[i]:indexes[i+1] for i in 1:(length(indexes)-1)),
        )
    end
end

@testset "Points of run-time length: dist, convex_comb, push_distinct!, split_edges" begin
    M = [0.0 3.0 3.0 6.0; 0.0 4.0 4.0 8.0; 1.0 1.0 1.0 1.0]
    S = PntSeq([Point{3,Float64}(c) for c in eachcol(M)])
    V = GPntSeq{AnyDim,Float64,Vector{Float64}}([Vector(c) for c in eachcol(M)])
    C = MatVecPntSeq(M)
    @test dist(V[1], V[2]) == dist(S[1], S[2]) == dist(C[1], C[2]) == 5.0
    @test dist_sq(V[1], C[2]) == 25.0
    @test @inferred(dist(V[1], V[2])) isa Float64
    # Measured in a specialized function: at top level Julia 1.10 boxes the result.
    dist_allocs(p, q) = @allocated dist(p, q)
    dist_allocs(V[1], V[2])
    @test dist_allocs(V[1], V[2]) == 0
    @test_throws DimensionMismatch dist([1.0, 2.0], [1.0, 2.0, 3.0])
    @test convex_comb(V[1], V[2], 0.5) ==
          [1.5, 2.0, 1.0] ==
          collect(convex_comb(S[1], S[2], 0.5))
    @test convex_comb(C[1], C[2], 0.5) isa Vector{Float64}
    @test geom_length(V) == geom_length(S) == geom_length(C) == 10.0
    @test prefix_lengths(V) ==
          prefix_lengths(S) ==
          prefix_lengths(C) ==
          [0.0, 5.0, 5.0, 10.0]

    W = GPntSeq{AnyDim,Float64,Vector{Float64}}()
    @test push_distinct!(W, [1.0, 2.0])
    @test !push_distinct!(W, [1.0, 2.0])
    @test_throws DimensionMismatch push_distinct!(W, [1.0, 2.0, 3.0])
    @test length(W) == 1

    for P in (V, C)
        Q = split_edges(P)
        @test Q isa GPntSeq{AnyDim,Float64,Vector{Float64}}
        @test collect(Q) == collect(split_edges(S))
    end
    @test size(M) == (3, 4)   # split_edges of the matrix sequence did not touch it
end

@testset "convex_comb of static points with different coordinate types stays static" begin
    c = convex_comb(Point{2,Int}(0, 0), Point{2,Float64}(1.0, 2.0), 0.5)
    @test c isa Point{2,Float64}
    @test c == Point{2,Float64}(0.5, 1.0)
end
