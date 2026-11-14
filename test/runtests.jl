using Test
using LinearAlgebra
using Random
using StaticArrays
using BasicCompGeometry
using BasicCompGeometry:
    Polygon_random_sphere,
    at,
    bisection_point,
    bound!,
    cardin,
    convex_coef,
    distance_infty,
    hausdorff_dist_subseg,
    height,
    match_price,
    middle,
    prefix_lengths,
    width

using Aqua

# Bytes allocated by `f(args...)`: the minimum of three measurements inside a
# function after a warm-up call. A single `@allocated` at test-set scope, or the
# first one after a warm-up on Julia 1.10, also counts the measurement's own
# overhead, so zero-allocation tests use this.
function allocs(f, args...)
    f(args...)
    return minimum(@allocated(f(args...)) for _ in 1:3)
end
Random.seed!(0xBADC0DE)
Aqua.test_all(BasicCompGeometry)

@testset "BasicCompGeometry.jl" begin

    @testset "Default exports and explicit generic imports" begin
        generic_names = (
            :area,
            :direction,
            :faces,
            :height,
            :width,
            :boundary,
            :complement,
            :clip,
            :depth,
            :lines,
            :locate,
            :vertices,
            :distance,
        )
        exported = Set(names(BasicCompGeometry))
        default_consumer = Module(:DefaultGeometryConsumer)
        Core.eval(default_consumer, :(using BasicCompGeometry))
        explicit_consumer = Module(:ExplicitGeometryConsumer)
        Core.eval(
            explicit_consumer,
            :(using BasicCompGeometry:
                area,
                direction,
                faces,
                height,
                width,
                boundary,
                complement,
                clip,
                depth,
                lines,
                locate,
                vertices,
                distance),
        )

        for name in generic_names
            @test name ∉ exported
            @test !isdefined(default_consumer, name)
            @test isdefined(BasicCompGeometry, name)
            @test getfield(explicit_consumer, name) === getfield(BasicCompGeometry, name)
        end
        for name in (:Point, :Segment, :BBox, :Polygon, :Circle, :Sphere, :convex_hull)
            @test getfield(default_consumer, name) === getfield(BasicCompGeometry, name)
        end

        # Another geometry package may export these same generic names. Check
        # both import orders without adding that package as a test dependency.
        other_geometry = Module(:OtherGeometry)
        for name in generic_names
            Core.eval(other_geometry, :(function $name end))
            Core.eval(other_geometry, :(export $name))
        end
        for bcg_first in (true, false)
            mixed_consumer = Module(:MixedGeometryConsumer)
            Core.eval(mixed_consumer, :(const OtherGeometry = $other_geometry))
            imports = (:(using BasicCompGeometry), :(using .OtherGeometry))
            for import_expr in (bcg_first ? imports : reverse(imports))
                Core.eval(mixed_consumer, import_expr)
            end
            for name in generic_names
                # Resolve dynamically created imports in their current world.
                @test Core.eval(mixed_consumer, name) === Core.eval(other_geometry, name)
            end
        end

        bb = BBox(point(0.0, 0.0), point(2.0, 3.0))
        @test explicit_consumer.width(bb) == 2.0
        @test BasicCompGeometry.height(bb) == 3.0
        @test BasicCompGeometry.depth === BBT.depth
    end

    @testset "Opt-in Extras public helpers" begin
        extra_names = (
            :AbsPolygon,
            :Canvas,
            :MatPolygon,
            :Points,
            :Polygon_random_sphere,
            :ReadWrite,
            :VecPntSeq2F,
            :VecPoint2I,
            :VecPolygon2F,
            :VertTrapezoid,
            :VirtArray,
            :a_coeff,
            :add_latex_macros!,
            :add_latex_packages!,
            :after_tangent_to_polygon,
            :alpha1,
            :alpha2,
            :append_latex_preamble!,
            :area,
            :at,
            :axis_direction,
            :b_coeff,
            :before_tangent_to_polygon,
            :bisection_point,
            :bound!,
            :boundary,
            :c_coeff,
            :cairo_draw_halfplane,
            :cairo_draw_halfplanes,
            :cairo_draw_latex,
            :cairo_draw_latex_page,
            :cairo_draw_points,
            :cairo_draw_polygon,
            :cairo_draw_setup,
            :cairo_set_line_width,
            :cardin,
            :center,
            :clip,
            :compare_right_arc_to_point,
            :compare_right_arcs,
            :complement,
            :convex_coef,
            :d_max,
            :d_min,
            :defaultview,
            :depth,
            :derivative,
            :description,
            :direction,
            :dist_point_segment,
            :dist_segment_segment,
            :dist_subspace,
            :distance,
            :distance_infty,
            :end_point,
            :expand!,
            :expand_add!,
            :faces,
            :get_file_path,
            :get_latex_preamble,
            :hausdorff_dist_subseg,
            :height,
            :init!,
            :latex_to_pdf,
            :lines,
            :locate,
            :match_price,
            :max_dist_subspace,
            :middle,
            :open_canvas,
            :p_coeff,
            :pdf_merge,
            :pnt_seq_random_sphere,
            :prefix_lengths,
            :push_distinct!,
            :r_major,
            :r_minor,
            :rand_polygon,
            :random_halfplane,
            :read_halfplanes,
            :read_latex_snippet,
            :reset_latex_preamble!,
            :rotation_angle,
            :set_latex_preamble!,
            :start_point,
            :split_edges,
            :vertex,
            :vertices,
            :width,
            :write_halfplanes,
            :write_plt,
        )
        default_names = (
            :AbsCurve2D,
            :AbsPntSeq,
            :AbstractKernel,
            :ArrVertDecomp,
            :ArrVertDecomp2F,
            :BBT,
            :BBox,
            :BBox2F,
            :Circle,
            :Circle2F,
            :Circle2I,
            :CircleArc,
            :CircleArc2F,
            :ConvexHull3D,
            :CubicBezier,
            :CubicBezier2F,
            :CubicBezier3F,
            :CubicSpline,
            :Curve2D,
            :Curve2DF,
            :CurvePolygon2D,
            :CurvePolygon2DF,
            :Ellipse,
            :Ellipse2F,
            :EllipticArc,
            :EllipticArc2F,
            :ExactKernel,
            :FloatingKernel,
            :Halfplane,
            :Halfplane2F,
            :Hyperbola,
            :IpeDraw,
            :Kd_tree,
            :Line,
            :LongestConvexSubset,
            :MVBB,
            :MatPntSeq,
            :MonotoneTriangulation,
            :MetricSpace,
            :Parabola,
            :ParabolaArc,
            :ParabolaArc2F,
            :Plane,
            :Plane2F,
            :PntSeq,
            :PntSeq2F,
            :PntSeq2I,
            :PntSeq3F,
            :Point,
            :Point2F,
            :Point2I,
            :Point3F,
            :Polygon,
            :Polygon2F,
            :Polygon2I,
            :Polygon3F,
            :PolygonTriangulation,
            :Segment,
            :Segment2F,
            :Sphere,
            :Sphere2F,
            :Sphere3F,
            :VerticalTrapezoid,
            :VerticalTrapezoid2F,
            :WSPD,
            :affine_interpolate,
            :angle_less,
            :apply_transform,
            :approx_diameter,
            :approx_diameter_subspace,
            :bbox,
            :bisector,
            :bottom_left,
            :boundary_line,
            :centroid,
            :compare_angle,
            :compare_distance,
            :compare_squared_distance,
            :compare_xy,
            :convex_comb,
            :convex_hull,
            :corner,
            :default_view,
            :diam,
            :dim,
            :dist,
            :dist_sq,
            :exact_diameter,
            :exact_diameter_subspace,
            :flatten,
            :geom_length,
            :hausdorff_simplify,
            :id_pair,
            :in_interior,
            :interpolate_catmull_rom,
            :interpolate_natural_spline,
            :intersect_hyperbolas,
            :intersect_line_curve,
            :intersect_lines,
            :intersections,
            :invert,
            :is_collinear,
            :is_convex,
            :is_inside,
            :is_left_eq_turn,
            :is_left_turn,
            :is_right_eq_turn,
            :is_right_turn,
            :is_star_shaped,
            :is_x_monotone,
            :is_y_monotone,
            :max_dist,
            :merge_chains,
            :midpoint,
            :minkowski_sum,
            :minkowski_sum_decomposition,
            :NoHit,
            :OverlapHit,
            :PointHit,
            :build_chain_index,
            :build_rmq,
            :chain_window,
            :first_positive_vertical_contact,
            :incoming_sector_contains,
            :nedges,
            :open_sector_contains,
            :range_argmin,
            :vertical_ray_segment,
            :AbsLatticeBody,
            :ConvexPolygonBody,
            :FloatChordBody,
            :LatticeDisk,
            :chord_range,
            :discrete_hull,
            :discrete_hull_columns,
            :discrete_hull_start,
            :find_edge_direction,
            :lattice_convex_hull,
            :monotone_chain_indices,
            :monotone_chains,
            :nn_point,
            :on_boundary,
            :on_segment,
            :orient,
            :orient_sign,
            :pairwise_intersections,
            :point,
            :point_max,
            :point_min,
            :point_on,
            :polar,
            :polar_line,
            :polygon,
            :polygon_diagonals,
            :polygon_kernel,
            :polygon_segments,
            :proper_crossing,
            :rand_gaussian,
            :rand_halfplane,
            :rand_pnt_seq,
            :rand_point,
            :rotation,
            :sample_uniformly,
            :scale!,
            :segment_at_x,
            :segment_key,
            :simplify,
            :stack_diagonals,
            :subdivide,
            :tangent_line,
            :top_right,
            :translate!,
            :translation,
            :trapezoid_chains,
            :trapezoids,
            :triangulate_ear_cut,
            :triangulate_ear_cut_kd_tree,
            :triangulate_ear_cut_u_kd_tree,
            :triangulate_ear_cut_z_order,
            :triangulate_monotone,
            :triangulate_polygon,
            :triangulate_seidel,
            :triangulate_sweep,
            :twice_area,
            :turn_sign,
            :union_bbox,
            :vertical_decomposition,
            :view_box,
            :visibility_polygon,
        )
        extras = BasicCompGeometry.Extras
        @test Set(names(extras)) == Set((:Extras, extra_names...))
        @test Set(names(BasicCompGeometry)) == Set((:BasicCompGeometry, default_names...))

        extras_consumer = Module(:ExtrasGeometryConsumer)
        Core.eval(extras_consumer, :(using BasicCompGeometry.Extras))
        for name in extra_names
            @test getfield(extras, name) === getfield(BasicCompGeometry, name)
            @test Core.eval(extras_consumer, name) === getfield(BasicCompGeometry, name)
        end
        @test !isdefined(extras_consumer, :Point)
        @test !isdefined(extras_consumer, :convex_hull)

        # Exercise both the full package path and the documented CG alias.
        for imports in (
            :(using BasicCompGeometry, BasicCompGeometry.Extras),
            quote
                import BasicCompGeometry as CG
                using .CG, .CG.Extras
            end,
        )
            consumer = Module(:CombinedGeometryConsumer)
            Core.eval(consumer, imports)
            for name in (extra_names..., default_names...)
                @test Core.eval(consumer, name) === getfield(BasicCompGeometry, name)
            end
            @test Core.eval(consumer, :(width(BBox(point(0.0, 0.0), point(2.0, 3.0))))) ==
                  2.0
            @test Core.eval(consumer, :(area(Ellipse(point(0.0, 0.0), 2.0, 3.0)))) ≈ 6π
        end
    end

    @testset "Points (2D & 3D, Float & Int)" begin
        # 2D Float
        p1 = point(0.0, 0.0)
        p2 = point(3.0, 4.0)
        @test dist(p1, p2) ≈ 5.0
        @test dist_sq(p1, p2) ≈ 25.0
        @test convex_comb(p1, p2, 0.5) ≈ point(1.5, 2.0)

        # 3D Int
        q1 = point(0, 0, 0)
        q2 = point(1, 2, 2)
        @test dist(q1, q2) ≈ 3.0
        @test dist_sq(q1, q2) == 9

        # Predicates
        a, b, c = point(0, 0), point(1, 0), point(0, 1)
        @test turn_sign(a, b, c) > 0
        @test is_left_turn(a, b, c) == true
        @test is_left_eq_turn(a, b, c) == true
        @test is_right_turn(a, c, b) == true
        @test is_left_eq_turn(a, c, b) == false
        @test is_right_eq_turn(a, c, b) == true
        @test turn_sign(a, c, b) < 0
        @test is_collinear(a, b, point(2, 0)) == true
        @test is_left_eq_turn(a, b, point(2, 0)) == true
        @test is_right_eq_turn(a, b, point(2, 0)) == true
        @test turn_sign(a, b, point(2, 0)) == 0

        # Random & Val(D) inference
        @test length(rand_point(2)) == 2
        @test length(rand_gaussian(3)) == 3
        @test (@inferred rand_point(Val(3))) isa Point{3,Float64}
        @test (@inferred rand_point(Val(2), Float32)) isa Point{2,Float32}
        @test (@inferred rand_gaussian(Val(3))) isa Point{3,Float64}
        @test (@inferred rand_gaussian(Val(2), Float32)) isa Point{2,Float32}
        let sample_points =
                () -> (rand_point(Val(3)), rand_gaussian(Val(3)), rand_gaussian(3))
            @test allocs(sample_points) == 0
        end

        # dist_subspace type stability across Int, Float32, Rational
        @test (@inferred BasicCompGeometry.dist_subspace(
            point(0, 0),
            point(0, 0),
            point(1, 0),
        )) === 0.0
        @test (@inferred BasicCompGeometry.dist_subspace(
            point(0, 3),
            point(0, 0),
            point(1, 0),
        )) === 3.0
        @test (@inferred BasicCompGeometry.dist_subspace(
            point(0.0f0, 0.0f0),
            point(0.0f0, 0.0f0),
            point(1.0f0, 0.0f0),
        )) === 0.0f0
        @test (@inferred BasicCompGeometry.dist_subspace(
            point(0 // 1, 0 // 1),
            point(0 // 1, 0 // 1),
            point(1 // 1, 0 // 1),
        )) === 0.0

        # Min/Max
        @test point_max(point(1, 5), point(2, 3)) == point(2, 5)
        @test point_min(point(1, 5), point(2, 3)) == point(1, 3)
    end

    @testset "Segments & Lines" begin
        p1 = point(0.0, 0.0)
        p2 = point(10.0, 0.0)
        seg = Segment(p1, p2)

        @test geom_length(seg) == 10.0
        @test at(seg, 0.2) == point(2.0, 0.0)

        # Point distance
        qr = point(5.0, 5.0)
        @test dist(seg, qr) == 5.0
        @test nn_point(seg, qr) == point(5.0, 0.0)
        @test convex_coef(seg, point(2.0, 0.0)) ≈ 0.2

        # Segment distance
        seg2 = Segment(point(5.0, 2.0), point(5.0, 10.0))
        @test dist(seg, seg2) == 2.0
        let seg_dist_bench =
                (a, b) -> BasicCompGeometry.dist_segment_segment(a.p, a.q, b.p, b.q)
            @test allocs(seg_dist_bench, seg, seg2) == 0
        end
        @test dist(
            Segment(point(0.0, 0.0), point(1e-4, 0.0)),
            Segment(point(5e-5, -5e-5), point(5e-5, 5e-5)),
        ) == 0.0
        @test dist(Segment(point(0, 0), point(2, 0)), Segment(point(0, 3), point(2, 3))) ==
              3.0
        @test (@inferred BasicCompGeometry.dist_segment_segment(
            point(0.0f0, 0.0f0),
            point(2.0f0, 0.0f0),
            point(1.0f0, 1.0f0),
            point(1.0f0, 3.0f0),
        )) === 1.0f0
        # Crossing at a tiny angle: must not fall back to the endpoint distance.
        let L = 1e4, θ = 1e-7
            @test dist(
                Segment(point(-L, 0.0), point(L, 0.0)),
                Segment(point(-L * cos(θ), -L * sin(θ)), point(L * cos(θ), L * sin(θ))),
            ) < 1e-9
        end
        @test BasicCompGeometry.dist_segment_segment(
            point(0.0, 0.0),
            point(0.0, 0.0),
            point(1.0, 1.0),
            point(1.0, 1.0),
        ) ≈ sqrt(2)

        # nn_point type stability across Int, Float32, Rational
        @test (@inferred nn_point(point(0, 0), point(0, 0), point(1, 1))) ===
              point(0.0, 0.0)
        @test (@inferred nn_point(point(0, 0), point(2, 0), point(1, 1))) ===
              point(1.0, 0.0)
        @test (@inferred nn_point(
            point(0.0f0, 0.0f0),
            point(2.0f0, 0.0f0),
            point(1.0f0, 1.0f0),
        )) === point(1.0f0, 0.0f0)
        @test (@inferred nn_point(
            point(0 // 1, 0 // 1),
            point(0 // 1, 0 // 1),
            point(1 // 1, 1 // 1),
        )) === point(0 // 1, 0 // 1)
        @test (@inferred nn_point(
            point(0 // 1, 0 // 1),
            point(2 // 1, 0 // 1),
            point(1 // 2, 1 // 1),
        )) === point(1 // 2, 0 // 1)

        # Bisection
        # p_a, p_b define a vertical bisector at x=5
        p_a, p_b = point(0.0, 0.0), point(10.0, 0.0)
        # s_cross from (0,0) to (10,0) crosses x=5 at (5,0)
        s_cross = Segment(point(0.0, 0.0), point(10.0, 0.0))
        ok, t, pt = bisection_point(s_cross, p_a, p_b)
        @test ok == true
        @test pt ≈ point(5.0, 0.0)

        # bisection_point type stability across Int, Float32, Rational
        @test (@inferred bisection_point(
            Segment(point(0, 0), point(0, 10)),
            point(0, 0),
            point(10, 0),
        )) === (false, 0.0, point(0.0, 0.0))
        @test (@inferred bisection_point(
            Segment(point(0, 0), point(10, 0)),
            point(0, 0),
            point(10, 0),
        )) === (true, 0.5, point(5.0, 0.0))
        @test (@inferred bisection_point(
            Segment(point(0.0f0, 0.0f0), point(10.0f0, 0.0f0)),
            point(0.0f0, 0.0f0),
            point(10.0f0, 0.0f0),
        )) === (true, 0.5f0, point(5.0f0, 0.0f0))
        @test (@inferred bisection_point(
            Segment(point(0 // 1, 0 // 1), point(0 // 1, 1 // 1)),
            point(0 // 1, 0 // 1),
            point(1 // 1, 0 // 1),
        )) === (false, 0 // 1, point(0 // 1, 0 // 1))
        @test (@inferred bisection_point(
            Segment(point(0 // 1, 0 // 1), point(1 // 1, 0 // 1)),
            point(0 // 1, 0 // 1),
            point(1 // 1, 0 // 1),
        )) === (true, 1 // 2, point(1 // 2, 0 // 1))

        # Lines
        line = Line(point(0.0, 0.0), point(1.0, 1.0))
        @test BasicCompGeometry.distance(point(0.0, 10.0), line) ≈
              dist(point(0.0, 10.0), point(5.0, 5.0))
    end

    @testset "PntSeqs" begin
        pnts = [point(0.0, 0.0), point(1.0, 0.0), point(1.0, 1.0), point(0.0, 1.0)]
        poly = PntSeq(pnts)

        @test length(poly) == 4
        @test dim(poly) == 2
        @test cardin(poly) == 4
        @test geom_length(poly) == 3.0 # Open curve by default

        # In-place ops
        translate!(poly, point(1.0, 1.0))
        @test poly[1] == point(1.0, 1.0)
        scale!(poly, 2.0)
        @test poly[1] == point(2.0, 2.0)

        # Analysis
        plens = prefix_lengths(poly)
        @test length(plens) == 4
        @test plens[end] == geom_length(poly)
        # poly vertices after translate (1,1) and scale 2.0:
        # (0,0)->(1,1)->(2,2)
        # (1,0)->(2,1)->(4,2)
        # (1,1)->(2,2)->(4,4)
        # (0,1)->(1,2)->(2,4)
        # lengths: dist((2,2),(4,2))=2, dist((4,2),(4,4))=2, dist((4,4),(2,4))=2. Total=6.
        # at(poly, 0.5) is at length 3. (4,2) is length 2, (4,4) is length 4.
        # So it is halfway between (4,2) and (4,4), which is (4,3).
        @test at(poly, 0.5) == point(4.0, 3.0)

        # Iteration & Base
        @test first(poly) == poly[1]
        @test last(poly) == poly[end]
        @test length(collect(poly)) == 4

        # Matrix
        mat = Matrix(poly)
        @test size(mat) == (2, 4)

        # Simplification
        poly_big = Polygon([point(x, 0.0) for x in 0:0.01:1])
        poly_sim, indices = simplify(poly_big, 0.1)
        @test length(poly_sim) < length(poly_big)

        # Sampling
        poly_uniform = sample_uniformly(poly, 10)
        @test length(poly_uniform) == 10
        @test geom_length(poly_uniform) ≈ geom_length(poly)

        # Type stability and O(m + n) allocation for dist, at, and sample_uniformly
        p32 = PntSeq([Point{2,Float32}(Float32(i), Float32(i % 3)) for i in 1:20])
        @test (@inferred dist(p32, 1, 1)) === 0.0f0
        @test (@inferred dist(p32, 1, 2)) isa Float32
        @test (@inferred at(p32, 0.0f0)) isa Point{2,Float32}
        @test (@inferred at(p32, 0.5f0)) isa Point{2,Float32}
        @test (@inferred at(p32, 1.0f0)) isa Point{2,Float32}
        p32_u = @inferred sample_uniformly(p32, 12)
        @test p32_u isa PntSeq{2,Float32}
        @test p32_u[1] == p32[1] && p32_u[end] == p32[end]

        pint = PntSeq([Point{2,Int}(0, 0), Point{2,Int}(2, 0), Point{2,Int}(2, 2)])
        @test (@inferred at(pint, 0.0)) === point(0.0, 0.0)
        @test (@inferred at(pint, 0.5)) === point(2.0, 0.0)
        @test (@inferred sample_uniformly(pint, 5)) isa PntSeq{2,Float64}

        p1k = PntSeq([point(Float64(i), sin(i)) for i in 1:1000])
        _run_su(p, n) = sample_uniformly(p, n)
        _run_su(p1k, 1000)
        @test (@allocated _run_su(p1k, 1000)) < 100_000

        # Random
        @test length(rand_pnt_seq(2, Float64, 5)) == 5
    end

    @testset "Bounding Boxes" begin
        p1 = point(-1.0, -2.0)
        p2 = point(5.0, 10.0)
        bb = BBox(p1, p2)

        @test width(bb) == 6.0
        @test height(bb) == 12.0
        @test middle(bb, 1) == 2.0
        @test diam(bb) ≈ dist(p1, p2)

        @test is_inside(point(0, 0), bb)
        @test !is_inside(point(10, 10), bb)

        # Padding
        bb2 = bb + 1.0
        @test width(bb2) == 8.0

        # Mutating
        bb_empty = BBox{2,Float64}()
        @test bb_empty.f_init == false
        bound!(bb_empty, point(0.0, 0.0))
        @test bb_empty.f_init == true
        @test diam(bb_empty) == 0.0

        # Box distance
        bb_a = BBox(point(0, 0), point(1, 1))
        bb_b = BBox(point(2, 0), point(3, 1))
        @test dist(bb_a, bb_b) == 1.0
        # max_dist is distance between (0,0) and (3,1) or (0,1) and (3,0)
        # sqrt(3^2 + 1^2) = sqrt(10)
        @test max_dist(bb_a, bb_b) ≈ sqrt(10.0)
    end

    @testset "Transforms 2D" begin
        p = point(1.0, 0.0)

        # Rotation 90 deg (Counter-Clockwise)
        m_rot = BasicCompGeometry.rotation(pi/2)
        p_rot = BasicCompGeometry.apply_transform(m_rot, p)
        # (1,0) rotated 90 deg CCW is (0,1)
        @test p_rot[1] ≈ 0.0 atol=1e-10
        @test p_rot[2] ≈ 1.0

        # Translation
        m_trans = BasicCompGeometry.translation(2.0, 3.0)
        p_trans = BasicCompGeometry.apply_transform(m_trans, p)
        @test p_trans == point(3.0, 3.0)

        # Zero-allocation transform & AbstractMatrix support
        let tf_bench = (m, pt) -> BasicCompGeometry.apply_transform(m, pt)
            @test allocs(tf_bench, m_trans, p) == 0
        end
        m_static = SMatrix{3,3,Float64}(m_trans)
        @test (@inferred apply_transform(m_static, p)) == point(3.0, 3.0)
        @test_throws DimensionMismatch apply_transform(zeros(2, 2), p)
    end

    @testset "Algorithms & Module Level" begin
        # Hausdorff
        pnts1 = [point(x, 0.0) for x in 0:10]
        pnts2 = [point(x, 0.1) for x in 0:10]
        poly1, poly2 = PntSeq(pnts1), PntSeq(pnts2)

        @test hausdorff_dist_subseg(poly1) == 0.0
        sim_p, sim_i = hausdorff_simplify(poly1, 0.5)
        @test length(sim_p) < length(poly1)

        # Module level
        @test distance_infty(poly1, poly2) ≈ 0.1
        @test centroid([point(0, 0), point(2, 2)]) == point(1, 1)

        # Match price
        price = match_price(point(0, 0), point(1, 0), point(0, 1), point(1, 1))
        @test price > 0

        # Convex Hull
        hull_pnts = [point(0, 0), point(1, 0), point(0, 1), point(0.5, 0.5), point(1, 1)]
        hull = convex_hull(hull_pnts)
        # Expected hull points (CCW starting from min-x): (0,0), (1,0), (1,1), (0,1)
        @test length(hull) == 4
        @test point(0, 0) in hull
        @test point(1, 0) in hull
        @test point(1, 1) in hull
        @test point(0, 1) in hull
        @test !(point(0.5, 0.5) in hull)
    end

    @testset "VirtArray" begin
        using BasicCompGeometry.VirtArray
        v = [10, 20, 30, 40]
        va = VArray(v)
        @test length(va) == 4
        @test va[1] == 10
        @test va[4] == 40

        swap!(va, 1, 4)
        @test va[1] == 40
        @test va[4] == 10
        @test orig_index(va, 1) == 4
        @test orig_index(va, 4) == 1

        @test collect(va) == [40, 20, 30, 10]
    end

    @testset "BBT (Bounding Box Tree)" begin
        using BasicCompGeometry.BBT
        pnts = [point(x, y) for x in 0:10 for y in 0:10]
        poly = PntSeq(pnts)
        tree = Tree_init(poly)

        @test tree.root !== nothing
        @test diam(tree.root.bb) > 0

        Tree_fully_expand(tree)
        @test depth(tree.root) > 1

        # Test original index retrieval (type-stable and zero-allocation)
        min_idx = @inferred get_min_orig_index(tree, tree.root)
        max_idx = @inferred get_max_orig_index(tree, tree.root)
        @test min_idx == 1
        @test max_idx == length(pnts)
        tree_uncached = Tree_init(poly)
        Tree_fully_expand(tree_uncached)
        @test (@allocated get_min_orig_index(tree_uncached, tree_uncached.root)) == 0
        @test (@allocated get_max_orig_index(tree_uncached, tree_uncached.root)) == 0

        # Test 1000 random points in 3D
        pnts3d = [rand_point(3) for _ in 1:1000]
        tree3d = Tree_init(PntSeq(pnts3d))
        Tree_fully_expand(tree3d)

        @test depth(tree3d.root) >= 10 # log2(1000) is approx 9.9
        @test tree3d.id_counter >= 1000 # Should have at least n nodes for a full expansion

        # Verify root BB contains all points
        root_bb = tree3d.root.bb
        all_inside = all(p -> is_inside(p, root_bb), pnts3d)
        @test all_inside == true

        # Kd_tree alias and parent-node upward navigation
        @test Kd_tree === BBT
        @test BasicCompGeometry.Kd_tree === BBT
        kdt = Kd_tree.Tree_init(poly)
        @test kdt.root !== nothing
        @test isnothing(kdt.root.parent)
        @test parent(kdt.root) === nothing
        @test Kd_tree.parent(kdt.root) === nothing
        @test Kd_tree.parent_node(kdt.root) === nothing
        @test Kd_tree.is_root(kdt.root)
        @test !isnothing(kdt.root.left)
        @test !isnothing(kdt.root.right)
        @test kdt.root.left.parent === kdt.root
        @test parent(kdt.root.left) === kdt.root
        @test Kd_tree.parent_node(kdt.root.left) === kdt.root
        @test !Kd_tree.is_root(kdt.root.left)
        @test kdt.root.right.parent === kdt.root

        # Lazy refinement sets parent links on newly split children
        refined = Kd_tree.Tree_refine_node(kdt, kdt.root.left)
        @test refined === kdt.root.left
        @test !isnothing(refined.left)
        @test !isnothing(refined.right)
        @test refined.left.parent === refined
        @test parent(refined.right) === refined
        @test parent(parent(refined.left)) === kdt.root

        # Verify every node after full expansion links back to root
        function verify_parent_links(node, expected_parent, root_node)
            @test node.parent === expected_parent
            @test parent(node) === expected_parent
            @test parent_node(node) === expected_parent
            @test is_root(node) == isnothing(expected_parent)
            if isnothing(node.left) && isnothing(node.right)
                curr = node
                steps = 0
                while !is_root(curr)
                    curr = parent(curr)
                    steps += 1
                    @test steps <= length(root_node.r)
                end
                @test curr === root_node
            else
                !isnothing(node.left) && verify_parent_links(node.left, node, root_node)
                !isnothing(node.right) && verify_parent_links(node.right, node, root_node)
            end
        end
        verify_parent_links(tree.root, nothing, tree.root)
        verify_parent_links(tree3d.root, nothing, tree3d.root)

        # Legacy 11-argument Node constructor compatibility (defaults parent to nothing)
        legacy_node = typeof(kdt.root)(
            kdt.root.bb,
            1:1,
            nothing,
            nothing,
            true,
            0.0,
            999,
            nothing,
            nothing,
            nothing,
            nothing,
        )
        @test isnothing(legacy_node.parent)
        @test is_root(legacy_node)

        # Flat array-backed 2D Kd-tree (FlatNode2D / FlatTree2D / Tree_init_flat)
        @test sizeof(Kd_tree.FlatNode2D{Float64}) == 64
        empty_ft = Kd_tree.Tree_init_flat(Point{2,Float64}[])
        @test isempty(empty_ft)
        @test length(empty_ft) == 0

        ft_pts = [
            Point(
                cos(0.37 * i) * (1.0 + 0.2 * (i % 5)),
                sin(0.37 * i) * (1.0 + 0.1 * (i % 7)),
            ) for i in 1:120
        ]
        ft = Kd_tree.Tree_init_flat(ft_pts; leaf_size = 8)
        ft_ps = Kd_tree.Tree_init_flat(PntSeq(ft_pts); leaf_size = 8)
        @test !isempty(ft)
        @test length(ft) == 120
        @test length(ft_ps.nodes) == length(ft.nodes)
        @test is_root(ft, 1)
        @test parent(ft, 1) == 0
        @test parent_node(ft, 1) == 0

        seen_slots = falses(120)
        for id in eachindex(ft.nodes)
            nd = ft.nodes[id]
            if id > 1
                @test !is_root(ft, id)
                pid = parent(ft, id)
                @test pid == parent_node(ft, id)
                @test 1 <= pid < id
                pnd = ft.nodes[pid]
                @test pnd.left == id || pnd.right == id
                @test pnd.r_start <= nd.r_start <= nd.r_end <= pnd.r_end
            end
            if nd.is_leaf
                @test (nd.r_end - nd.r_start + 1) <= 8 || (nd.x1 == nd.x0 && nd.y1 == nd.y0)
                for slot in nd.r_start:nd.r_end
                    @test !seen_slots[slot]
                    seen_slots[slot] = true
                    orig = ft.perm[slot]
                    @test ft.leaf_of[orig] == id
                    @test ft.slot_of[orig] == slot
                    @test ft.px[slot] == ft_pts[orig].x
                    @test ft.py[slot] == ft_pts[orig].y
                    @test nd.x0 <= ft.px[slot] <= nd.x1
                    @test nd.y0 <= ft.py[slot] <= nd.y1
                end
            else
                lnd = ft.nodes[nd.left]
                rnd = ft.nodes[nd.right]
                @test lnd.r_start == nd.r_start
                @test lnd.r_end + 1 == rnd.r_start
                @test rnd.r_end == nd.r_end
            end
        end
        @test all(seen_slots)

        # Subset build + integer and Rational{BigInt} 10^320 coordinates
        sub_px = Float64[p.x for p in ft_pts[1:2:120]]
        sub_py = Float64[p.y for p in ft_pts[1:2:120]]
        sub_ids = collect(1:2:120)
        ft_sub = Kd_tree.Tree_init_flat(sub_px, sub_py, sub_ids, 120; leaf_size = 6)
        @test length(ft_sub) == 60
        @test count(!iszero, ft_sub.leaf_of) == 60
        @test all(i -> ft_sub.leaf_of[i] == 0 && ft_sub.slot_of[i] == 0, 2:2:120)

        int_pts = [Point(i, (i * 17) % 41) for i in 1:40]
        ft_int = Kd_tree.Tree_init_flat(int_pts; leaf_size = 4)
        @test length(ft_int) == 40
        @test is_root(ft_int, 1)

        huge_scale = big(10)^320
        big_pts = [
            Point{2,Rational{BigInt}}(
                big(i) * huge_scale // 1,
                big((i * 13) % 29) * huge_scale // 1,
            ) for i in 1:40
        ]
        ft_big = Kd_tree.Tree_init_flat(big_pts; leaf_size = 4)
        @test length(ft_big) == 40
        @test is_root(ft_big, 1)

        # Unshrunk-cell 96-byte FlatUTree2D checks
        @test sizeof(Kd_tree.FlatUNode2D{Float64}) == 96
        empty_fut = Kd_tree.Tree_init_flat_u(Point{2,Float64}[])
        @test isempty(empty_fut) && length(empty_fut) == 0
        fut_u = Kd_tree.Tree_init_flat_u(ft_pts; leaf_size = 8)
        @test fut_u isa Kd_tree.FlatUTree2D{Float64}
        @test length(fut_u) == length(ft_pts)
        @test is_root(fut_u, 1)
        @test Kd_tree.parent(fut_u, 1) == 0
        @test Kd_tree.parent_node(fut_u, 1) == 0
        @test fut_u.nodes[1].ux0 == -Inf && fut_u.nodes[1].ux1 == Inf
        @test fut_u.nodes[1].uy0 == -Inf && fut_u.nodes[1].uy1 == Inf
        @test fut_u.nodes[1].bound_mask == 0x00
        for nid in eachindex(fut_u.nodes)
            nd = fut_u.nodes[nid]
            for k in nd.r_start:nd.r_end
                x = fut_u.px[k]
                y = fut_u.py[k]
                @test nd.ux0 < x <= nd.ux1
                @test nd.uy0 < y <= nd.uy1
            end
            if !is_root(fut_u, nid)
                @test Kd_tree.parent(fut_u, nid) > 0
            end
        end

        fut_sub = Kd_tree.Tree_init_flat_u(
            [1.0, 3.0, 2.0, 4.0],
            [4.0, 2.0, 3.0, 1.0],
            Int32[2, 5, 7, 10],
            12;
            leaf_size = 2,
        )
        @test length(fut_sub) == 4
        @test fut_sub.n_total == 12
        @test fut_sub.leaf_of[5] != 0
        @test fut_sub.leaf_of[1] == 0

        fut_seq = Kd_tree.Tree_init_flat_u(PntSeq(ft_pts); leaf_size = 8)
        @test length(fut_seq) == length(ft_pts)

        fut_int = Kd_tree.Tree_init_flat_u(int_pts; leaf_size = 4)
        @test length(fut_int) == 40
        @test is_root(fut_int, 1)
        for nd in fut_int.nodes
            bm = nd.bound_mask
            for k in nd.r_start:nd.r_end
                @test ((bm & 0x01) == 0 || nd.ux0 < fut_int.px[k])
                @test ((bm & 0x02) == 0 || fut_int.px[k] <= nd.ux1)
                @test ((bm & 0x04) == 0 || nd.uy0 < fut_int.py[k])
                @test ((bm & 0x08) == 0 || fut_int.py[k] <= nd.uy1)
            end
        end

        fut_big = Kd_tree.Tree_init_flat_u(big_pts; leaf_size = 4)
        @test length(fut_big) == 40
        @test is_root(fut_big, 1)
        for nd in fut_big.nodes
            bm = nd.bound_mask
            for k in nd.r_start:nd.r_end
                @test ((bm & 0x01) == 0 || nd.ux0 < fut_big.px[k])
                @test ((bm & 0x02) == 0 || fut_big.px[k] <= nd.ux1)
                @test ((bm & 0x04) == 0 || nd.uy0 < fut_big.py[k])
                @test ((bm & 0x08) == 0 || fut_big.py[k] <= nd.uy1)
            end
        end

        # Multi-point leaves: nearest-neighbor queries must scan every point of a leaf.
        nn_rng = Xoshiro(7)
        nn_pts = PntSeq([point(rand(nn_rng), rand(nn_rng)) for _ in 1:200])
        nn_tree = Tree_init(nn_pts; leaf_size = 8)
        Tree_fully_expand(nn_tree)
        nn_ok = true
        for _ in 1:100
            q = point(rand(nn_rng), rand(nn_rng))
            d_tree, _, _ = BBT.approx_nn(nn_tree, q, 1.0)
            d_scan, _, _ = BBT.exact_naive_scan(nn_pts, q)
            d_tree == d_scan || (nn_ok = false)
        end
        @test nn_ok
        @test get_min_orig_index(nn_tree, nn_tree.root) == 1
        @test get_max_orig_index(nn_tree, nn_tree.root) == 200
        @test_throws ArgumentError Tree_init(nn_pts; leaf_size = 0)

        # Parent pointers must not make `show` dump the whole tree.
        @test length(sprint(show, nn_tree.root.left)) < 200
        @test length(sprint(show, nn_tree)) < 200

        # Flat trees reject non-finite coordinates (they used to hang) and bad original ids.
        @test_throws ArgumentError Kd_tree.Tree_init_flat(
            [point(-Inf, 0.0), point(0.0, 0.0), point(Inf, 0.0)];
            leaf_size = 1,
        )
        @test_throws ArgumentError Kd_tree.Tree_init_flat_u(
            [point(0.0, 0.0), point(NaN, 0.0), point(1.0, 0.0)];
            leaf_size = 1,
        )
        @test_throws ArgumentError Kd_tree.Tree_init_flat(
            [1.0, 2.0, 3.0],
            [1.0, 2.0, 3.0],
            [1, 2, 50],
            3;
            leaf_size = 1,
        )
        @test_throws ArgumentError Kd_tree.Tree_init_flat(
            [1.0, 2.0, 3.0],
            [1.0, 2.0, 3.0],
            [2, 2, 1],
            3;
            leaf_size = 1,
        )
        @test_throws ArgumentError Kd_tree.Tree_init_flat_u(
            [1.0, 2.0],
            [1.0, 2.0],
            [0, 1],
            2;
            leaf_size = 1,
        )

        # Integer coordinate spans beyond typemax must not overflow into a single huge leaf.
        wide_rng = Xoshiro(3)
        wide_pts = [
            point(
                rand(wide_rng, (typemin(Int), typemax(Int) - 1, typemax(Int))),
                rand(wide_rng, 1:10),
            ) for _ in 1:50
        ]
        for ft_wide in (
            Kd_tree.Tree_init_flat(wide_pts; leaf_size = 4),
            Kd_tree.Tree_init_flat_u(wide_pts; leaf_size = 4),
        )
            @test all(
                nd ->
                    !nd.is_leaf ||
                    nd.r_end - nd.r_start + 1 <= 4 ||
                    (nd.x0 == nd.x1 && nd.y0 == nd.y1),
                ft_wide.nodes,
            )
        end
    end

    @testset "Diameter Algorithms" begin
        # 2D points on a circle
        pnts = Polygon_random_sphere(2, Float64, 100)

        exact_diam = exact_diameter(pnts)
        @test exact_diam ≈ 2.0 atol=0.1 # Unit sphere diameter is 2.0

        # (1+ε) approximation with ε=0.1
        approx_diam = approx_diameter(pnts, 0.1)
        @test approx_diam <= exact_diam
        @test approx_diam * 1.1 >= exact_diam

        # Float32 point sequences with Float64 ε
        pnts32 = PntSeq([Point{2,Float32}(p.x, p.y) for p in pnts])
        approx_diam32 = approx_diameter(pnts32, 0.1)
        @test approx_diam32 isa Float32
        @test approx_diam32 <= exact_diam + 1.0f-5
        @test approx_diam32 * 1.1f0 >= exact_diam - 1.0f-5
    end

    @testset "WSPD (Well-Separated Pairs Decomposition)" begin
        using BasicCompGeometry.WSPD
        using BasicCompGeometry.VirtArray: orig_index

        function covered_point_pairs(W, pair)
            left = [orig_index(W.tree.PS, i) for i in pair.left.r]
            right = [orig_index(W.tree.PS, i) for i in pair.right.r]
            return [(min(i, j), max(i, j)) for i in left for j in right if i != j]
        end

        function test_wspd(points, sep)
            W = WSPD.init(PntSeq(points), sep)
            finals = WSPD.expand!(W)
            @test !isempty(finals)

            for pair in finals
                @test pair.dist > 0
                @test WSPD.separation(pair) <= sep + 1e-12
            end

            covered = reduce(vcat, (covered_point_pairs(W, pair) for pair in finals))
            expected = [(i, j) for i in 1:(length(points)-1) for j in (i+1):length(points)]
            @test sort(covered) == expected
            @test length(unique(covered)) == length(covered)
            return length(finals)
        end

        line = [point(Float64(x), 0.0) for x in 1:10]
        cloud = [point(randn(), randn(), randn()) for _ in 1:24]
        @test test_wspd(line, 2.0) < length(line)^2
        @test test_wspd(cloud, 0.5) < length(cloud)^2
        @test test_wspd(reverse(cloud), 0.5) < length(cloud)^2

        # Type stability for Float32 points and mismatched sep types (e.g. Int or Float64)
        line32 = [Point{2,Float32}(Float32(x), 0.0f0) for x in 1:10]
        W32 = @inferred WSPD.init(PntSeq(line32), 2)
        @test isconcretetype(typeof(W32))
        @test typeof(W32.pnts) == PntSeq{2,Float32}
        @test W32.sep === 2.0f0
        @test (@inferred WSPD.separation(WSPD.get_top(W32))) isa Float32
        finals32 = @inferred WSPD.expand!(W32)
        @test !isempty(finals32)
        @test (@inferred WSPD.separation(first(finals32))) isa Float32
    end

    @testset "10-Dimensional Geometry" begin
        # 10D Points
        p10_1 = point(zeros(10)...)
        p10_2 = point(ones(10)...)
        @test dist(p10_1, p10_2) ≈ sqrt(10.0)
        @test dist_sq(p10_1, p10_2) ≈ 10.0

        # 10D Segments
        seg10 = Segment(p10_1, p10_2)
        @test geom_length(seg10) ≈ sqrt(10.0)
        @test at(seg10, 0.5) ≈ point(fill(0.5, 10)...)

        # 10D BBox
        bb10 = BBox(p10_1, p10_2)
        @test diam(bb10) ≈ sqrt(10.0)
        @test is_inside(point(fill(0.1, 10)...), bb10)
        @test !is_inside(point(fill(1.1, 10)...), bb10)

        # 10D PntSeq
        poly10 = PntSeq([p10_1, p10_2, p10_1])
        @test cardin(poly10) == 3
        @test geom_length(poly10) ≈ 2 * sqrt(10.0)
    end

    @testset "Subspace Diameter Exactness (ε=0.0)" begin
        using Random
        rng = Random.Xoshiro(42)
        for _ in 1:100
            n = 100
            pts = [point(rand(rng), rand(rng), rand(rng)) for _ in 1:n]
            ps = PntSeq(pts)
            dir = normalize(point(randn(rng), randn(rng), randn(rng)))

            exact = exact_diameter_subspace(ps, dir)
            wspd = approx_diameter_subspace(ps, 0.0, dir)
            fst = BasicCompGeometry.MVBB.approx_diam_subspace(ps, 0.0, dir).distance

            @test wspd ≈ exact
            @test fst ≈ exact
        end
    end

    include("test_kernels.jl")
    include("test_exact_geometry.jl")
    include("test_monotone_triangulation.jl")
    include("test_visibility_polygon.jl")
    include("test_minkowski_sum.jl")
    include("test_discrete_hull.jl")
    include("test_chain_incidence.jl")
    include("test_pnt_seq_hierarchy.jl")
    include("test_curve_ops.jl")
    include("test_metric_space.jl")
    include("test_read_write.jl")
    include("test_longest_convex_subset.jl")
    include("test_convex_hull_3d.jl")
    include("test_mat_pnt_seq.jl")
    include("test_bbt_ann.jl")
    include("test_mvbb.jl")
    include("test_diameter_exactness.jl")
    include("test_spheres.jl")
    include("test_curve_polygon.jl")
    include("test_hyperbolas.jl")
    include("test_parabolas.jl")
    include("test_halfplanes.jl")
    include("test_bundles.jl")
    include("test_arrangements.jl")

    if Base.find_package("Cairo") !== nothing
        include("test_canvas.jl")
    end

    if Base.find_package("LaTeXStrings") !== nothing &&
       Base.find_package("Cairo") !== nothing
        include("test_latex_cairo.jl")
    end

    include("test_ellipses.jl")
    include("test_bezier.jl")
    include("test_ipedraw.jl")
    include("test_golden.jl")
    include("test_docstrings.jl")

end
