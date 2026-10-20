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
    zigzag = PntSeq([point(Float64(i), isodd(i) ? 1.0 : 0.0) for i = 0:30])
    for tolerance in (0.1, 0.6, 1.0)
        simplified, indexes = hausdorff_simplify(zigzag, tolerance)
        @test first(indexes) == 1 && last(indexes) == length(zigzag)
        @test collect(simplified) == [zigzag[i] for i in indexes]
        @test all(
            r -> BasicCompGeometry.hausdorff_dist_subseg(zigzag, r) <= tolerance,
            (indexes[i]:indexes[i+1] for i = 1:(length(indexes)-1)),
        )
    end
end
