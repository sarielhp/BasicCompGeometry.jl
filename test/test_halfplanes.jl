using Test
using LinearAlgebra
using Random
using BasicCompGeometry

@testset "Halfplanes" begin
    # 1. Constructors and basic properties
    p = point(0.0, 0.0)
    u = point(1.0, 0.0) # Horizontal line along x-axis
    l = Line(p, u)
    h_default = Halfplane(l)

    @test boundary(h_default) == l
    @test boundary_line(h_default) == l
    @test h_default.boundary == l

    # By right-halfplane convention:
    # Points with y < 0 are right turn (inside)
    # Points with y = 0 are on boundary (inside)
    # Points with y > 0 are left turn (outside)
    @test is_inside(point(2.0, -1.0), h_default) == true
    @test is_inside(point(5.0, 0.0), h_default) == true
    @test is_inside(point(2.0, 1.0), h_default) == false
    @test (point(2.0, -1.0) in h_default) == true
    @test (point(2.0, 1.0) in h_default) == false
    @test is_inside(h_default, point(2.0, -1.0)) == true

    # Interior & Boundary
    @test in_interior(point(2.0, -1.0), h_default) == true
    @test in_interior(point(2.0, 0.0), h_default) == false
    @test on_boundary(point(2.0, 0.0), h_default) == true
    @test on_boundary(point(2.0, 1.0), h_default) == false

    # Distance
    @test distance(point(2.0, -3.0), h_default) == 0.0
    @test distance(point(2.0, 3.0), h_default) ≈ 3.0
    @test distance(h_default, point(2.0, 3.0)) ≈ 3.0

    # 2. Defining halfplane using point z
    # Point z = (0, 1) is in upper halfplane (y >= 0)
    z_upper = point(0.0, 1.0)
    h_upper = Halfplane(l, z_upper)
    @test is_inside(z_upper, h_upper) == true
    @test is_inside(point(10.0, 2.0), h_upper) == true
    @test is_inside(point(10.0, -2.0), h_upper) == false

    # Point z = (0, -1) is in lower halfplane (y <= 0)
    z_lower = point(0.0, -1.0)
    h_lower = Halfplane(l, z_lower)
    @test is_inside(z_lower, h_lower) == true
    @test is_inside(point(10.0, -2.0), h_lower) == true
    @test is_inside(point(10.0, 2.0), h_lower) == false

    # Constructing with p, u, z directly
    h_upper_direct = Halfplane(p, u, z_upper)
    @test is_inside(z_upper, h_upper_direct) == true
    @test is_inside(point(0.0, -1.0), h_upper_direct) == false

    # Point z on the line must throw ArgumentError
    z_on_line = point(5.0, 0.0)
    @test_throws ArgumentError Halfplane(l, z_on_line)
    @test_throws ArgumentError Halfplane(p, u, z_on_line)

    # 3. Complement / Negation
    h_comp = complement(h_default)
    @test is_inside(point(0.0, 1.0), h_comp) == true
    @test is_inside(point(0.0, -1.0), h_comp) == false
    @test (-h_default) == h_comp

    # Boundary points belong to both
    @test is_inside(point(3.0, 0.0), h_default) == true
    @test is_inside(point(3.0, 0.0), h_comp) == true

    # 4. Depth operation
    # Create 4 halfplanes bounding [0, 1] x [0, 1]:
    # x >= 0  => Line((0,0), (0,1)) with z=(1,0)
    # x <= 1  => Line((1,0), (0,1)) with z=(0,0)
    # y >= 0  => Line((0,0), (1,0)) with z=(0,1)
    # y <= 1  => Line((0,1), (1,0)) with z=(0,0)
    h_left   = Halfplane(Line(point(0.0, 0.0), point(0.0, 1.0)), point(1.0, 0.0))
    h_right  = Halfplane(Line(point(1.0, 0.0), point(0.0, 1.0)), point(0.0, 0.0))
    h_bottom = Halfplane(Line(point(0.0, 0.0), point(1.0, 0.0)), point(0.0, 1.0))
    h_top    = Halfplane(Line(point(0.0, 1.0), point(1.0, 0.0)), point(0.0, 0.0))

    H = [h_left, h_right, h_bottom, h_top]

    # Point inside unit square: contained in all 4 halfplanes
    q_in = point(0.5, 0.5)
    @test depth(H, q_in) == 4
    @test depth(q_in, H) == 4

    # Point with x=1.5, y=0.5: contained in left, bottom, top (3 halfplanes)
    q_out1 = point(1.5, 0.5)
    @test depth(H, q_out1) == 3

    # Point with x=1.5, y=1.5: contained in left, bottom (2 halfplanes)
    q_out2 = point(1.5, 1.5)
    @test depth(H, q_out2) == 2

    # Point with x=-1.5, y=-1.5: contained in right, top (2 halfplanes)
    q_out3 = point(-1.5, -1.5)
    @test depth(H, q_out3) == 2

    # 5. Random halfplane generation
    rng = Random.Xoshiro(42)
    rh1 = rand_halfplane(rng)
    @test rh1 isa Halfplane{Float64}
    # Point p should be in unit square [0, 1]^2
    @test 0.0 <= rh1.boundary.p[1] <= 1.0
    @test 0.0 <= rh1.boundary.p[2] <= 1.0
    # Direction u should be approximately unit length
    @test norm(rh1.boundary.u) ≈ 1.0

    # Test random_halfplane alias and Base.rand
    rh2 = random_halfplane(rng)
    @test rh2 isa Halfplane{Float64}

    rh_vec = rand(rng, Halfplane, 20)
    @test length(rh_vec) == 20
    @test all(h -> h isa Halfplane{Float64}, rh_vec)

    rh2f_vec = rand(rng, Halfplane2F, 10)
    @test length(rh2f_vec) == 10
    @test all(h -> h isa Halfplane2F, rh2f_vec)

    # Verify depth with random halfplanes
    q_test = point(0.5, 0.5)
    d = depth(rh_vec, q_test)
    @test 0 <= d <= 20
    @test d == count(h -> is_inside(q_test, h), rh_vec)

    # 6. Polygon clipping by halfplane
    # Clip unit square [0,1]x[0,1] by halfplane y <= 0.5
    square = [point(0.0, 0.0), point(1.0, 0.0), point(1.0, 1.0), point(0.0, 1.0)]
    h_half = Halfplane(Line(point(0.0, 0.5), point(1.0, 0.0))) # Right halfplane is y <= 0.5
    clipped = clip(square, h_half)
    @test length(clipped) == 4
    # All vertices of clipped polygon must be inside h_half
    @test all(pt -> is_inside(pt, h_half), clipped)
    # Highest y-coordinate in clipped polygon should be 0.5
    @test maximum(pt[2] for pt in clipped) ≈ 0.5

    # Clip PntSeq
    ps_square = PntSeq(square)
    ps_clipped = clip(ps_square, h_half)
    @test ps_clipped isa PntSeq{2, Float64}
    @test length(ps_clipped) == 4

    # 7. File I/O: write_halfplanes and read_halfplanes
    tmp_path = joinpath(mktempdir(), "test_hps.txt")
    test_hps = [h_left, h_right, h_bottom, h_top]
    write_halfplanes(tmp_path, test_hps)
    @test isfile(tmp_path)

    loaded_hps = read_halfplanes(tmp_path)
    @test length(loaded_hps) == 4
    for i in 1:4
        @test loaded_hps[i].boundary.p ≈ test_hps[i].boundary.p
        @test loaded_hps[i].boundary.u ≈ test_hps[i].boundary.u
    end
    rm(tmp_path; force = true)
end
