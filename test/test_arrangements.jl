using Test
using LinearAlgebra
using Random
using BasicCompGeometry

@testset "Arrangement Vertical Decomposition" begin
    # 1. Zero lines in view [0, 1] x [0, 1]
    view_unit = BBox(point(0.0, 0.0), point(1.0, 1.0))
    decomp0 = vertical_decomposition(Line{2,Float64}[], view_unit)
    @test length(decomp0) == 1
    @test area(decomp0) ≈ 1.0
    @test isempty(decomp0.vertices)

    # 2. Single horizontal line y = 0.5
    l_horiz = Line(point(0.0, 0.5), point(1.0, 0.0))
    decomp1 = vertical_decomposition([l_horiz], view_unit)
    @test length(decomp1) == 2
    @test area(decomp1) ≈ 1.0
    @test isempty(decomp1.vertices)
    @test locate(point(0.5, 0.2), decomp1) !== nothing
    @test locate(point(0.5, 0.8), decomp1) !== nothing

    # 3. Two crossing lines: y = x and y = 1 - x
    l_diag1 = Line(point(0.0, 0.0), point(1.0, 1.0))
    l_diag2 = Line(point(0.0, 1.0), point(1.0, -1.0))
    lines_x = [l_diag1, l_diag2]

    decomp_x = vertical_decomposition(lines_x, view_unit)
    # 2 crossing lines inside unit square yield 6 vertical trapezoids
    @test length(decomp_x) == 6
    @test area(decomp_x) ≈ 1.0
    @test length(decomp_x.vertices) == 1
    @test decomp_x.vertices[1] ≈ point(0.5, 0.5)

    # Every trapezoid has positive or non-negative area and valid coordinates
    for trap in trapezoids(decomp_x)
        @test trap.x_left <= trap.x_right
        @test area(trap) > 0.0
        @test length(vertices(trap)) == 4
        @test polygon(trap) isa PntSeq{2,Float64}
    end

    # Point location tests
    t_left = locate(point(0.1, 0.5), decomp_x)
    @test t_left !== nothing
    @test is_inside(point(0.1, 0.5), t_left)

    t_top = locate(point(0.5, 0.9), decomp_x)
    @test t_top !== nothing
    @test is_inside(point(0.5, 0.9), t_top)

    @test locate(point(2.0, 2.0), decomp_x) === nothing

    # 4. Three lines forming a triangle
    l1 = Line(point(0.0, 0.2), point(1.0, 0.0)) # y = 0.2
    l2 = Line(point(0.2, 0.0), point(0.0, 1.0)) # x = 0.2
    l3 = Line(point(0.0, 1.0), point(1.0, -1.0)) # x + y = 1

    lines3 = [l1, l2, l3]
    decomp3 = vertical_decomposition(lines3, view_unit)
    @test area(decomp3) ≈ 1.0
    @test length(decomp3.vertices) == 3

    # Halfplane bundle test
    hps3 = [Halfplane(l, point(0.3, 0.3)) for l in lines3]
    decomp_hp = vertical_decomposition(hps3, view_unit)
    @test length(decomp_hp) == length(decomp3)
    @test area(decomp_hp) ≈ 1.0

    # 5. Auto-view constructor
    decomp_auto = ArrVertDecomp(hps3; factor = 1.2)
    @test decomp_auto isa ArrVertDecomp{Float64}
    view_area = width(decomp_auto.view) * height(decomp_auto.view)
    @test area(decomp_auto) ≈ view_area

    # Accessors and interface
    @test length(decomp_auto) > 0
    @test decomp_auto[1] isa VerticalTrapezoid{Float64}
    @test view_box(decomp_auto) == decomp_auto.view
    @test lines(decomp_auto) == decomp_auto.lines

    # 6. Random halfplanes test
    rng = Random.Xoshiro(12345)
    rand_hps = [rand_halfplane(rng) for _ = 1:6]
    rand_view = default_view(rand_hps; factor = 1.1)
    rand_decomp = ArrVertDecomp(rand_hps, rand_view)

    # Invariant: sum of trapezoid areas must equal view rectangle area
    target_area = width(rand_view) * height(rand_view)
    @test isapprox(area(rand_decomp), target_area, atol = 1e-8)

    # Every trapezoid corners form counter-clockwise polygon
    for trap in rand_decomp
        @test trap.x_left <= trap.x_right
        @test area(trap) >= -1e-12
        # corners CCW orientation test:
        # bottom-left -> bottom-right -> top-right -> top-left
        c = trap.corners
        @test c[1][1] ≈ trap.x_left atol = 1e-9
        @test c[2][1] ≈ trap.x_right atol = 1e-9
        @test c[3][1] ≈ trap.x_right atol = 1e-9
        @test c[4][1] ≈ trap.x_left atol = 1e-9
        @test c[1][2] <= c[4][2] + 1e-9
        @test c[2][2] <= c[3][2] + 1e-9
    end

    # 7. Line & Halfplane clipping to BBox
    bb_clip = BBox(point(0.0, 0.0), point(10.0, 10.0))
    l_diag = Line(point(0.0, 0.0), point(1.0, 1.0))
    seg_diag = BasicCompGeometry.clip(l_diag, bb_clip)
    @test seg_diag isa Segment{2,Float64}
    @test seg_diag.p ≈ point(0.0, 0.0)
    @test seg_diag.q ≈ point(10.0, 10.0)

    l_miss = Line(point(20.0, 0.0), point(0.0, 1.0))
    @test BasicCompGeometry.clip(l_miss, bb_clip) === nothing

    hp_clip = Halfplane(l_diag)
    seg_hp = BasicCompGeometry.clip(hp_clip, bb_clip)
    @test seg_hp isa Segment{2,Float64}
    @test seg_hp.p ≈ point(0.0, 0.0)
    @test seg_hp.q ≈ point(10.0, 10.0)

    # 8. Trapezoid chains and face reconstruction
    chains_x = trapezoid_chains(decomp_x)
    @test length(chains_x) == 4 # 2 crossing lines in box yield 4 faces
    faces_x = faces(decomp_x)
    @test length(faces_x) == 4
    for f in faces_x
        @test f isa PntSeq{2,Float64}
        @test length(f) >= 3
    end

    # Random arrangement face count & partition invariant
    rand_chains = trapezoid_chains(rand_decomp)
    rand_faces = faces(rand_decomp)
    @test length(rand_chains) == length(rand_faces)
    @test length(rand_faces) > 0
    # Every face is non-empty and has >= 3 vertices
    for f in rand_faces
        @test length(f) >= 3
    end
    # The faces function directly from bundle
    faces_direct = faces(rand_hps, rand_view)
    @test length(faces_direct) == length(rand_faces)
end
