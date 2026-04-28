using BasicCompGeometry
using BasicCompGeometry.BBT
using BasicCompGeometry.WSPD
using Test
using LinearAlgebra

@testset "MatPolygon" begin
    D, N = 3, 10
    M = rand(D, N)
    mp = MatPolygon(M)
    
    @testset "Basics" begin
        @test mp isa AbsPolygon{3, Float64}
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
        
        # Test that updating the polygon (if possible) updates the matrix
        # Since Point is SVector (immutable), we have to set the whole point
        # But setindex! is defined for AbsPolygon
        new_p = Point{3, Float64}(1.0, 2.0, 3.0)
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
