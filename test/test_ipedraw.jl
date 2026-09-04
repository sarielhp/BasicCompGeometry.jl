using Test
using BasicCompGeometry
using BasicCompGeometry.IpeDraw

@testset "IpeDraw Vector Figure Interface" begin
    @test isdefined(BasicCompGeometry, :IpeDraw)
    temp_dir = mktempdir()

    try
        # 1. Canvas creation and sizing
        canvas = IpeCanvas(width=500.0, height=300.0)
        @test BasicCompGeometry.width(canvas) == 500.0
        @test BasicCompGeometry.height(canvas) == 300.0
        @test canvas.width == 500.0
        @test canvas.height == 300.0

        # 2. Layers and Views
        set_layer!(canvas, "background")
        add_layer!(canvas, "foreground")
        add_view!(canvas, ["background", "foreground"])
        @test "background" in canvas.layers
        @test "foreground" in canvas.layers

        # 3. Preamble management
        add_preamble!(canvas, raw"\newcommand{\eps}{\varepsilon}")
        @test occursin(raw"\newcommand{\eps}{\varepsilon}", canvas.preamble)

        # 4. Geometric Primitives
        p1 = Point(50.0, 50.0)
        p2 = Point(200.0, 150.0)
        draw_point!(canvas, p1; stroke=:darkred, fill=:red, shape=:disk)
        draw_point!(canvas, 100.0, 80.0; stroke=:blue, shape=:circle)
        draw_segment!(canvas, p1, p2; stroke=:black, pen=:heavier)
        draw_segment!(canvas, Segment(p1, p2); stroke=:gray4)

        bb = BBox(Point(20.0, 20.0), Point(250.0, 200.0))
        draw_box!(canvas, bb; stroke=:blue, fill=:lightblue)
        draw_box!(canvas, 10.0, 10.0, 50.0, 50.0)

        poly = [Point(10.0, 10.0), Point(50.0, 10.0), Point(30.0, 40.0)]
        draw_polygon!(canvas, poly; close=true, fill=:lightgreen)

        c = Circle(Point(100.0, 100.0), 30.0)
        draw_circle!(canvas, c; stroke=:darkgreen)
        draw_circle!(canvas, 80.0, 80.0, 15.0; stroke=:darkblue)

        draw_arc!(canvas, Point(100.0, 100.0), 25.0, 0.0, π / 2; stroke=:darkred)

        # 5. Conceptual & Algorithmic Helpers
        draw_bar!(canvas, 50.0, 450.0, 120.0;
            label_left = "0",
            label_right = raw"M = \frac{1}{\eps}"
        )
        draw_span!(canvas, 150.0, 300.0, 120.0; fill=:lightgreen, stroke=:darkgreen)
        draw_dimension!(canvas, 150.0, 300.0, 150.0;
            label = raw"\Delta \le 2\eps \mu",
            arrow = :both
        )
        draw_arrow!(canvas, Point(100.0, 60.0), Point(150.0, 110.0))
        draw_curved_arrow!(canvas, Point(50.0, 200.0), Point(200.0, 200.0); bend=25.0)

        draw_label!(canvas, 250.0, 50.0, raw"\mathcal{E}_{\le i} \implies |X - \mu| \le \eps \mu";
            size = :large,
            halign = :center
        )

        # LaTeXStrings support test if available
        if Base.find_package("LaTeXStrings") !== nothing
            @eval begin
                using LaTeXStrings
                draw_label!($canvas, 250.0, 80.0, L"\mathcal{E}_{\le i} \implies |X - \mu| \le \eps \mu")
            end
        end

        # 5b. Ellipses, Arcs, Béziers, Splines, Holes, and Groups
        e = Ellipse(Point(300.0, 200.0), 40.0, 20.0, π / 6)
        draw_ellipse!(canvas, e; fill=:lightblue, stroke=:darkblue, tiling=:hatch)
        draw_ellipse!(canvas, Point(350.0, 200.0), 30.0, 15.0; stroke=:red)

        c_arc = CircleArc(Point(100.0, 100.0), 25.0, 0.0, π)
        draw_arc!(canvas, c_arc; stroke=:darkgreen)

        el_arc = EllipticArc(e, 0.0, π / 2)
        draw_arc!(canvas, el_arc; stroke=:magenta, arrow=:forward)

        bez = CubicBezier(Point(50.0, 50.0), Point(50.0, 100.0), Point(100.0, 100.0), Point(100.0, 50.0))
        draw_bezier!(canvas, bez; stroke=:purple, pen=:heavier)

        knots = PntSeq([Point(200.0, 100.0), Point(220.0, 140.0), Point(260.0, 110.0), Point(300.0, 150.0)])
        draw_spline!(canvas, knots; stroke=:darkblue, pen=:fat)
        draw_bspline!(canvas, knots; stroke=:gray3, dash=:dashed)

        outer_box = PntSeq([Point(50.0, 50.0), Point(150.0, 50.0), Point(150.0, 150.0), Point(50.0, 150.0)])
        inner_hole = PntSeq([Point(75.0, 75.0), Point(125.0, 75.0), Point(125.0, 125.0), Point(75.0, 125.0)])
        draw_polygon_with_holes!(canvas, outer_box, inner_hole; fill=:lightgray, stroke=:black)

        ipe_group(canvas; matrix=[1.0, 0.0, 0.0, 1.0, 10.0, 20.0], opacity=Symbol("50%")) do cv
            draw_point!(cv, 0.0, 0.0; stroke=:red)
            draw_segment!(cv, Point(0.0, 0.0), Point(20.0, 20.0))
        end

        # 6. File serialization and PDF compilation
        base_fig = joinpath(temp_dir, "test_fig")
        artifacts = export_figure(canvas, base_fig;
            caption = "Test caption referencing the lemma.",
            label = "fig:test_sample"
        )

        @test isfile(artifacts.ipe)
        @test isfile(artifacts.tex)
        @test occursin("<ipe version=", read(artifacts.ipe, String))
        @test occursin(raw"\figlab{fig:test_sample}", read(artifacts.tex, String))

        # Check PDF compilation via ipetoipe
        if Sys.which("ipetoipe") !== nothing
            @test isfile(artifacts.pdf)
            @test filesize(artifacts.pdf) > 0
        end

        # 7. Block syntax `open_ipe`
        res = open_ipe(joinpath(temp_dir, "block_fig"); caption="Block test", label="fig:block") do cv
            draw_box!(cv, 0.0, 0.0, 100.0, 100.0; stroke=:darkred)
            draw_label!(cv, 50.0, 50.0, raw"x \in \mathcal{S}")
        end
        @test isfile(res.ipe)
        @test isfile(res.tex)
        if Sys.which("ipetoipe") !== nothing
            @test isfile(res.pdf)
        end

    finally
        rm(temp_dir, recursive=true, force=true)
    end
end
