using Test
using BasicCompGeometry
using BasicCompGeometry.MVBB: rotating_calipers_min_area, OBBox, volume, approx_diam, approx_mvbb
using StaticArrays
using LinearAlgebra

@testset "MVBB Tests" begin
    @testset "2D Rotating Calipers" begin
        # Rectangle points
        pnts = [
            point(0.0, 0.0),
            point(2.0, 0.0),
            point(2.0, 1.0),
            point(0.0, 1.0)
        ]
        hull = convex_hull(PntSeq(pnts))
        bbox2d = rotating_calipers_min_area(hull)
        @test bbox2d.area ≈ 2.0 atol=1e-6
        
        # Rotated rectangle
        θ = pi/6
        rot = [cos(θ) -sin(θ); sin(θ) cos(θ)]
        rpnts = [Point{2,Float64}(rot * p) for p in pnts]
        rhull = convex_hull(PntSeq(rpnts))
        rbbox2d = rotating_calipers_min_area(rhull)
        @test rbbox2d.area ≈ 2.0 atol=1e-6
    end

    @testset "3D approx_diam" begin
        # Points on a line
        pnts = [point(Float64(i), 0.0, 0.0) for i in 1:100]
        pp = approx_diam(PntSeq(pnts), 0.01)
        @test pp.distance ≈ 99.0 atol=1e-6
        
        # Random points in a cube
        n = 1000
        pnts_cube = [rand_point(3) for _ in 1:n]
        pp_cube = approx_diam(PntSeq(pnts_cube), 0.1)
        exact_d = exact_diameter(PntSeq(pnts_cube))
        @test pp_cube.distance <= exact_d
        @test pp_cube.distance * (1 + 0.1) >= exact_d
    end

    @testset "3D approx_mvbb" begin
        # Points in a rotated box
        # Box: [0, 10] x [0, 2] x [0, 1]
        pnts = [point(rand()*10, rand()*2, rand()) for _ in 1:1000]
        
        # Add corners
        push!(pnts, point(0.0, 0.0, 0.0))
        push!(pnts, point(10.0, 0.0, 0.0))
        push!(pnts, point(0.0, 2.0, 0.0))
        push!(pnts, point(0.0, 0.0, 1.0))
        push!(pnts, point(10.0, 2.0, 1.0))
        
        # Rotate the points
        # Use some arbitrary rotation
        v = normalize(point(1.0, 2.0, 3.0))
        θ = pi/4
        # Rodrigues' rotation formula
        function rotate_vec(p, k, theta)
            return p*cos(theta) + (k × p)*sin(theta) + k*dot(k, p)*(1 - cos(theta))
        end
        
        rpnts = [rotate_vec(p, v, θ) for p in pnts]
        
        ε = 0.1
        obb = approx_mvbb(PntSeq(rpnts), ε)
        
        vol = volume(obb)
        # Ideal volume is 10 * 2 * 1 = 20
        @test vol >= 20.0 - 1e-6
        # Approximate MVBB might be slightly larger, but should be reasonable
        @test vol < 30.0 # Loose bound
    end
end
