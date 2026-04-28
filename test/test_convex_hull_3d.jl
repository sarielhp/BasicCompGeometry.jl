using Test
using BasicCompGeometry
using BasicCompGeometry.ConvexHull3D
using LinearAlgebra

@testset "ConvexHull3D" begin
    @testset "Cube" begin
        # 8 corners of a cube
        pnts = [
            point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0),
            point(1.0, 1.0, 0.0), point(0.0, 1.0, 0.0),
            point(0.0, 0.0, 1.0), point(1.0, 0.0, 1.0),
            point(1.0, 1.0, 1.0), point(0.0, 1.0, 1.0)
        ]
        # Add some interior points
        push!(pnts, point(0.5, 0.5, 0.5))
        push!(pnts, point(0.2, 0.8, 0.3))

        verts, faces = convex_hull_3d(pnts)
        
        # All 8 corners should be in the hull
        @test length(verts) == 8
        @test all(i <= 8 for i in verts)
        
        # Number of faces for a cube (triangulated)
        # 6 faces * 2 triangles/face = 12 triangles
        @test length(faces) == 12
        
        # Euler's formula: V - E + F = 2
        # For cube: 8 - 12 + 6 = 2
        # Our triangles: V=8, F=12, so E must be 18?
        # V=8, F=12 -> 8 - E + 12 = 2 -> E = 18.
        # A triangulated cube has 12 original edges + 6 diagonal edges = 18. Correct.
    end

    @testset "Tetrahedron" begin
        pnts = [
            point(0.0, 0.0, 0.0),
            point(1.0, 0.0, 0.0),
            point(0.0, 1.0, 0.0),
            point(0.0, 0.0, 1.0)
        ]
        verts, faces = convex_hull_3d(pnts)
        @test length(verts) == 4
        @test length(faces) == 4
    end

    @testset "Random Sphere" begin
        # Points on a sphere should all be on the hull
        n = 100
        pnts = [point(normalize(randn(3))...) for _ in 1:n]
        verts, faces = convex_hull_3d(pnts)
        
        # For points on a sphere, most/all should be on the hull
        @test length(verts) >= n * 0.9
    end
end
