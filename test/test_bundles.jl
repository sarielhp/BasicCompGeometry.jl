using Test
using LinearAlgebra
using BasicCompGeometry

@testset "Line and Halfplane Bundles" begin
    # 1. Line-line intersection
    l1 = Line(point(0.0, 0.0), point(1.0, 0.0))  # y = 0
    l2 = Line(point(0.0, 0.0), point(0.0, 1.0))  # x = 0
    l3 = Line(point(1.0, 0.0), point(-1.0, 1.0)) # x + y = 1

    pt12 = intersect(l1, l2)
    @test pt12 !== nothing
    @test pt12 ≈ point(0.0, 0.0)

    pt13 = intersect(l1, l3)
    @test pt13 !== nothing
    @test pt13 ≈ point(1.0, 0.0)

    pt23 = intersect(l2, l3)
    @test pt23 !== nothing
    @test pt23 ≈ point(0.0, 1.0)

    # Parallel lines
    l_par = Line(point(0.0, 1.0), point(1.0, 0.0)) # y = 1
    @test intersect(l1, l_par) === nothing
    @test intersect_lines(l1, l_par) === nothing

    # 2. Halfplane-halfplane intersection
    h1 = Halfplane(l1, point(0.0, -1.0))
    h2 = Halfplane(l2, point(1.0, 0.0))
    @test intersect(h1, h2) ≈ point(0.0, 0.0)

    # 3. Vertices of a line bundle
    lines = [l1, l2, l3]
    verts = vertices(lines)
    @test length(verts) == 3
    expected_verts = [point(0.0, 0.0), point(1.0, 0.0), point(0.0, 1.0)]
    for ev in expected_verts
        @test any(v -> dist(v, ev) < 1e-10, verts)
    end

    # Vertices of a halfplane bundle
    h3 = Halfplane(l3, point(0.0, 0.0))
    hps = [h1, h2, h3]
    verts_hp = vertices(hps)
    @test length(verts_hp) == 3
    for ev in expected_verts
        @test any(v -> dist(v, ev) < 1e-10, verts_hp)
    end

    # Degenerate: concurrent lines through origin
    l_diag = Line(point(0.0, 0.0), point(1.0, 1.0))
    concurrent_lines = [l1, l2, l_diag]
    verts_all = vertices(concurrent_lines; unique = false)
    @test length(verts_all) == 3
    @test all(v -> dist(v, point(0.0, 0.0)) < 1e-10, verts_all)

    verts_uniq = vertices(concurrent_lines; unique = true)
    @test length(verts_uniq) == 1
    @test dist(verts_uniq[1], point(0.0, 0.0)) < 1e-10

    # Bundle with fewer than 2 elements
    @test isempty(vertices(Line{2,Float64}[]))
    @test isempty(vertices([l1]))

    # Type safety: mixed bundles should throw ArgumentError
    @test_throws ArgumentError vertices(Any[l1, h1])

    # 4. Bounding box of bundle
    bb_lines = BBox(lines)
    @test bb_lines.f_init == true
    @test bottom_left(bb_lines) ≈ point(0.0, 0.0)
    @test top_right(bb_lines) ≈ point(1.0, 1.0)
    @test bbox(lines) == bb_lines

    bb_hps = BBox(hps)
    @test bottom_left(bb_hps) ≈ point(0.0, 0.0)
    @test top_right(bb_hps) ≈ point(1.0, 1.0)

    # Empty / parallel bundles
    bb_empty = BBox(Line{2,Float64}[])
    @test bb_empty.f_init == false

    bb_par = BBox([l1, l_par])
    @test bb_par.f_init == false

    # 5. Default view of bundle
    dv = default_view(lines; factor = 1.1)
    @test dv.f_init == true
    # Original extent is [0, 1] x [0, 1], width=1, height=1, center=(0.5, 0.5)
    # Scaled by 1.1: half_size = 0.55 -> [-0.05, 1.05]
    @test bottom_left(dv)[1] ≈ -0.05 atol = 1e-10
    @test bottom_left(dv)[2] ≈ -0.05 atol = 1e-10
    @test top_right(dv)[1] ≈ 1.05 atol = 1e-10
    @test top_right(dv)[2] ≈ 1.05 atol = 1e-10
    @test width(dv) ≈ 1.1 atol = 1e-10
    @test height(dv) ≈ 1.1 atol = 1e-10

    # All bundle vertices must lie strictly inside the default view
    for v in verts
        @test is_inside(v, dv) == true
    end

    # Halfplane bundle default view
    dv_hp = default_view(hps)
    @test dv_hp == dv
    @test defaultview(hps) == dv

    # Two lines intersecting at a single point: width/height = 0 should expand with padding
    dv_two = default_view([l1, l2])
    @test dv_two.f_init == true
    @test width(dv_two) > 0.0
    @test height(dv_two) > 0.0
    @test is_inside(point(0.0, 0.0), dv_two) == true
end
