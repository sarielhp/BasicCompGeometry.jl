using Test
using BasicCompGeometry
import BasicCompGeometry: after_tangent_to_polygon, at, before_tangent_to_polygon
using BasicCompGeometry: direction
using LinearAlgebra

@testset "CurvePolygon2D and point_on tests" begin
    # Test CircleArc geom_length and at
    arc = CircleArc(point(0.0, 0.0), 1.0, 0.0, pi / 2)
    @test geom_length(arc) ≈ pi / 2
    @test at(arc, 0.0) ≈ point(1.0, 0.0)
    @test at(arc, 1.0) ≈ point(0.0, 1.0)
    @test at(arc, 0.5) ≈ point(cos(pi / 4), sin(pi / 4))

    # Test Segment geom_length and at
    seg = Segment(point(0.0, 1.0), point(0.0, 0.0))
    @test geom_length(seg) == 1.0
    @test at(seg, 0.5) ≈ point(0.0, 0.5)

    # Test CurvePolygon2D
    seg2 = Segment(point(0.0, 0.0), point(1.0, 0.0))
    cp = CurvePolygon2D([seg2, arc, seg])
    @test cp isa CurvePolygon2D{Float64}
    @test geom_length(cp) ≈ 1.0 + pi / 2 + 1.0

    # Test point_on for CurvePolygon2D
    @test point_on(cp, 0.0) ≈ point(0.0, 0.0)
    @test point_on(cp, 1.0) ≈ point(0.0, 0.0) # end of seg

    total_l = geom_length(cp)
    t_arc_mid = (1.0 + (pi / 4)) / total_l
    @test point_on(cp, t_arc_mid) ≈ point(cos(pi / 4), sin(pi / 4))

    # Test point_on for regular polygon (PntSeq / AbsPntSeq)
    pnts = [point(0.0, 0.0), point(2.0, 0.0), point(2.0, 2.0), point(0.0, 2.0)]
    poly = PntSeq(pnts)
    # Closed parameterization: total length 8.0
    @test point_on(poly, 0.0, closed = true) ≈ point(0.0, 0.0)
    @test point_on(poly, 0.25, closed = true) ≈ point(2.0, 0.0)
    @test point_on(poly, 0.5, closed = true) ≈ point(2.0, 2.0)
    @test point_on(poly, 0.75, closed = true) ≈ point(0.0, 2.0)
    # Test direction and tangent_line
    v = direction(pi / 4)
    @test v ≈ point(cos(pi / 4), sin(pi / 4))

    line_t, q_t = tangent_line(poly, 0.0)
    @test q_t ≈ point(2.0, 0.0) # max x vertex

    u_t = before_tangent_to_polygon(cp, poly, 0.0)
    v_t = after_tangent_to_polygon(cp, poly, 0.0)
    @test 0.0 <= u_t <= 1.0
    @test 0.0 <= v_t <= 1.0
end

@testset "ParabolaArc tests" begin
    # Directrix y = -1, focus (0, 1) -> y = x^2 / 4, vertex at (0, 0), p = 1
    para = Parabola(point(0.0, 1.0), Line(point(0.0, -1.0), point(1.0, 0.0)))
    arc = ParabolaArc(para, point(-2.0, 1.0), point(2.0, 1.0))
    @test at(arc, 0.0) ≈ point(-2.0, 1.0)
    @test at(arc, 1.0) ≈ point(2.0, 1.0)
    @test at(arc, 0.5) ≈ point(0.0, 0.0)
    @test geom_length(arc) > 4.0 # chord is 4.0, curved arc > 4.0

    # Line intersections
    h_line = Line(point(0.0, 0.5), point(1.0, 0.0)) # y = 0.5
    pts = intersect_line_curve(h_line, arc)
    @test length(pts) == 2
    for (pt, s) in pts
        @test pt[2] ≈ 0.5
        @test -2.0 <= pt[1] <= 2.0
        @test 0.0 <= s <= 1.0
    end
end

@testset "CurvePolygon2D is_inside tests (Circle & Square)" begin
    # Unit Circle CurvePolygon2D
    a1 = CircleArc(point(0.0, 0.0), 1.0, 0.0, pi / 2)
    a2 = CircleArc(point(0.0, 0.0), 1.0, pi / 2, pi)
    a3 = CircleArc(point(0.0, 0.0), 1.0, pi, 3pi / 2)
    a4 = CircleArc(point(0.0, 0.0), 1.0, 3pi / 2, 2pi)
    cp_circle = CurvePolygon2D([a1, a2, a3, a4])

    @test is_inside(point(0.0, 0.0), cp_circle) == true
    @test is_inside(point(0.5, 0.5), cp_circle) == true
    @test is_inside(point(1.5, 0.5), cp_circle) == false
    @test is_inside(point(-1.2, 0.0), cp_circle) == false

    # Unit Square CurvePolygon2D
    s1 = Segment(point(0.0, 0.0), point(1.0, 0.0))
    s2 = Segment(point(1.0, 0.0), point(1.0, 1.0))
    s3 = Segment(point(1.0, 1.0), point(0.0, 1.0))
    s4 = Segment(point(0.0, 1.0), point(0.0, 0.0))
    cp_sq = CurvePolygon2D([s1, s2, s3, s4])

    @test is_inside(point(0.5, 0.5), cp_sq) == true
    @test is_inside(point(0.1, 0.9), cp_sq) == true
    @test is_inside(point(-0.1, 0.5), cp_sq) == false
    @test is_inside(point(1.1, 0.5), cp_sq) == false
end

@testset "Point vs Square Voronoi Curved Polygon (100,000 random points)" begin
    using Random

    # Helper functions as in examples/sq_vs_point.jl
    function square_boundary_lines()
        left = Line(point(0.0, 0.0), point(0.0, 1.0))
        right = Line(point(1.0, 0.0), point(0.0, 1.0))
        bottom = Line(point(0.0, 0.0), point(1.0, 0.0))
        top = Line(point(0.0, 1.0), point(1.0, 0.0))
        return [left, right, bottom, top]
    end

    function compute_voronoi_curved_polygon(p::Point{2,T}) where {T}
        lines = square_boundary_lines()
        n = 4
        parabolas = [Parabola(p, lines[i]) for i in 1:n]
        interior = point(T(0.5), T(0.5))
        normals = Vector{Point{2,T}}(undef, n)
        for i in 1:n
            ni = Point{2,T}(-lines[i].u[2], lines[i].u[1])
            ni = ni / norm(ni)
            if dot(interior - lines[i].p, ni) < 0
                ni = -ni
            end
            normals[i] = ni
        end
        candidates = Point{2,T}[]
        for i in 1:n
            for j in (i+1):n
                nd = normals[i] - normals[j]
                if norm(nd) < 1e-12
                    continue
                end
                ci = dot(lines[i].p, normals[i])
                cj = dot(lines[j].p, normals[j])
                d = ci - cj
                dir = Point{2,T}(-nd[2], nd[1])
                pt0 = Point{2,T}(0.0, 0.0)
                if abs(nd[1]) > 1e-12
                    pt0 = Point{2,T}(d / nd[1], 0.0)
                elseif abs(nd[2]) > 1e-12
                    pt0 = Point{2,T}(0.0, d / nd[2])
                end
                bis = Line(pt0, dir)
                pts = intersect_line_curve(bis, parabolas[i])
                for (pt, _) in pts
                    d_p = dist(pt, p)
                    ok = true
                    for kk in 1:n
                        sd = dot(pt - lines[kk].p, normals[kk])
                        if sd < d_p - 1e-8
                            ok = false
                            break
                        end
                    end
                    if ok && abs(dot(pt - lines[i].p, normals[i]) - d_p) < 1e-8
                        push!(candidates, pt)
                    end
                end
            end
        end
        sort!(candidates, by = pt -> atan(pt[2] - p[2], pt[1] - p[1]))
        unique_pts = Point{2,T}[]
        for pt in candidates
            if isempty(unique_pts) || dist(pt, unique_pts[end]) > 1e-8
                push!(unique_pts, pt)
            end
        end

        m = length(unique_pts)
        arcs = ParabolaArc{T}[]
        for idx in 1:m
            p_cur = unique_pts[idx]
            p_next = unique_pts[mod1(idx + 1, m)]
            mid = (p_cur + p_next) / 2
            best_k = 0
            best_diff = Inf
            for k in 1:n
                d = abs(dist(mid, p) - BasicCompGeometry.distance(mid, lines[k]))
                if d < best_diff
                    best_diff = d
                    best_k = k
                end
            end
            push!(arcs, ParabolaArc(parabolas[best_k], p_cur, p_next))
        end
        return CurvePolygon2D(arcs)
    end

    # Test with site p = (0.3, 0.4) as in sq_vs_point.jl
    p = point(0.3, 0.4)
    cp = compute_voronoi_curved_polygon(p)
    @test length(cp.curves) == 4
    @test geom_length(cp) > 0.0

    Random.seed!(42)
    n_samples = 100_000
    pts = [point(rand(), rand()) for _ in 1:n_samples]

    # Verify that for all 100,000 points, distance test and is_inside test are identical
    mismatches = 0
    for q in pts
        dist_to_site = dist(q, p)
        dist_to_boundary = min(q.x, 1.0 - q.x, q.y, 1.0 - q.y)
        inside_by_dist = dist_to_site < dist_to_boundary
        inside_by_lib = is_inside(q, cp)
        if inside_by_dist != inside_by_lib
            mismatches += 1
        end
    end
    @test mismatches == 0

    # Also test points strictly outside the square
    outside_pts = [point(-0.2, 0.5), point(1.2, 0.5), point(0.5, -0.2), point(0.5, 1.2)]
    for q in outside_pts
        @test is_inside(q, cp) == false
    end
end

@testset "CurvePolygon2D & point_on type stability and zero-allocations" begin
    curves = BasicCompGeometry.Curve2DF[]
    for _ in 1:50
        push!(curves, Segment(Point2F(0.0, 0.0), Point2F(1.0, 0.0)))
        push!(curves, CircleArc(Point2F(0.0, 0.0), 1.0, 0.0, pi / 2))
    end
    cp100 = CurvePolygon2D(curves)
    q = Point2F(0.5, 0.5)
    @inferred is_inside(q, cp100)
    is_inside_alloc(pt, poly) = @allocated is_inside(pt, poly)
    is_inside_alloc(q, cp100)
    @test is_inside_alloc(q, cp100) == 0

    ps = PntSeq([Point2F(0, 0), Point2F(1, 0), Point2F(1, 1), Point2F(0, 1)])
    res_ti = @inferred BasicCompGeometry.tangent_intersections_cp(cp100, ps, 0.4)
    @test length(res_ti) == 5
    ti_alloc(poly, qpoly) =
        @allocated BasicCompGeometry.tangent_intersections_cp(poly, qpoly, 0.4)[3]
    ti_alloc(cp100, ps)
    @test ti_alloc(cp100, ps) == 0

    # Tangent hit (1 point) and 2-hit circle arc intersection
    arc = CircleArc(Point2F(0.0, 0.0), 1.0, 0.0, Float64(pi))
    @test length(intersect_line_curve(Line(Point2F(0.0, 1.0), Point2F(1.0, 0.0)), arc)) == 1
    @test length(intersect_line_curve(Line(Point2F(0.0, 0.5), Point2F(1.0, 0.0)), arc)) == 2
    @test isempty(intersect_line_curve(Line(Point2F(0.0, 2.0), Point2F(1.0, 0.0)), arc))

    # point_on closed type stability and zero allocations
    po_alloc(seq, t) = @allocated point_on(seq, t; closed = true)[1]
    po_alloc(ps, 0.5)
    @test po_alloc(ps, 0.5) == 0

    ps_i = PntSeq([Point2I(0, 0), Point2I(2, 0), Point2I(2, 2), Point2I(0, 2)])
    @test (@inferred point_on(ps_i, 0.0; closed = true)) === Point2F(0.0, 0.0)
    @test (@inferred point_on(ps_i, 0.5; closed = true)) ≈ Point2F(2.0, 2.0)
    @test (@inferred point_on(ps_i, 1.0; closed = true)) === Point2F(0.0, 0.0)

    ps_f32 = PntSeq([
        Point{2,Float32}(0, 0),
        Point{2,Float32}(2, 0),
        Point{2,Float32}(2, 2),
        Point{2,Float32}(0, 2),
    ])
    @test (@inferred point_on(ps_f32, 0.5f0; closed = true)) isa Point{2,Float32}
    @test (@inferred point_on(PntSeq([Point2I(1, 2)]), 0.5; closed = true)) ===
          Point2F(1.0, 2.0)

    # Degenerate zero-length edge and explicitly closed sequence
    ps_deg = PntSeq([Point2F(0, 0), Point2F(0, 0), Point2F(1, 0), Point2F(0, 0)])
    @test point_on(ps_deg, 0.5; closed = true) ≈ Point2F(1.0, 0.0)
    @test point_on(ps_deg, 0.0; closed = true) ≈ Point2F(0.0, 0.0)
end
