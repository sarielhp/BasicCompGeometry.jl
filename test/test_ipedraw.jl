using Test
using BasicCompGeometry
using BasicCompGeometry.IpeDraw

@testset "IpeDraw Vector Figure Interface" begin
    @test isdefined(BasicCompGeometry, :IpeDraw)
    temp_dir = mktempdir()
    test_external = get(ENV, "BCG_TEST_EXTERNAL_TOOLS", "false") == "true"

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

        # Compact API, reusable styles, fitting, and page-space labels
        compact = IpeCanvas(width=200.0, height=100.0)
        unit_circle = Circle(point(0.0, 0.0), 1.0)
        fit!(compact, unit_circle; margin=30.0)
        circle_style = Style(fill=:lightblue, fill_opacity=0.2, stroke=:blue, pen=:heavier)
        draw!(compact, unit_circle; style=circle_style)
        mark!(compact, point(0.0, 0.0))
        label!(compact, point(0.0, 1.0), "p"; offset=(4.0, 5.0), anchor=:southwest)
        compact_xml = BasicCompGeometry.IpeDraw.to_xml(compact)
        @test compact.viewport.scale == 20.0
        @test occursin("20.000 0 0 20.000 100.000 50.000 e", compact_xml)
        @test occursin("opacity=\"20%\"", compact_xml)
        @test occursin("pos=\"104.000 75.000\"", compact_xml)

        layer(compact, :annotations) do cv
            mark!(cv, point(0.5, 0.0); stroke=:red)
        end
        @test compact.active_layer == "alpha"
        with_style(compact, Style(stroke=:darkred, pen=:heavier)) do cv
            draw!(cv, Segment(point(-0.5, 0.0), point(0.5, 0.0)))
        end
        @test isempty(compact.active_style.values)
        @test occursin("stroke=\"darkred\" pen=\"heavier\"", compact.elements[end])

        # Scoped viewports, clipping, themes, and page-space furniture
        theme = publication_theme()
        @test theme.region === theme[:region]
        @test theme.region.values.fill_opacity == 0.2
        @test_throws ArgumentError PageBox(0, 0, 0, 10)

        previous_viewport = compact.viewport
        detail_box = PageBox(110, 10, 80, 70)
        inset(compact, detail_box; fit=unit_circle, margin=5) do cv
            draw!(cv, unit_circle; style=theme.region, fill=:lightgreen, stroke=:darkgreen)
        end
        @test compact.viewport === previous_viewport

        clipped_box = BBox(point(-0.5, -0.5), point(0.5, 0.5))
        clip_to(compact, clipped_box) do cv
            draw!(cv, unit_circle; stroke=:red)
        end
        legend!(compact, ["disk" => Style(fill=:lightblue, stroke=:blue)];
                at=(5, 55), width=70)
        scale_bar!(compact, 0.5; at=(10, 12), label="0.5")
        furniture_xml = BasicCompGeometry.IpeDraw.to_xml(compact)
        @test length(findall(" clip=", furniture_xml)) == 2
        @test occursin("disk", furniture_xml)
        @test occursin("0.5", furniture_xml)
        @test compact.viewport === previous_viewport
        @test_throws ArgumentError legend!(compact, ["bad" => :blue])

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
            label = "fig:test_sample",
            outputs = (:ipe, :tex),
        )

        @test isfile(artifacts.ipe)
        @test isfile(artifacts.tex)
        @test occursin("<ipe version=", read(artifacts.ipe, String))
        @test occursin(raw"\figlab{fig:test_sample}", read(artifacts.tex, String))

        if test_external
            pdf_artifacts = export_figure(canvas, base_fig; outputs=(:pdf,))
            @test isfile(pdf_artifacts.pdf)
            @test filesize(pdf_artifacts.pdf) > 0
        end

        # 7. Block syntax `open_ipe`
        res = open_ipe(joinpath(temp_dir, "block_fig");
            caption="Block test", label="fig:block", outputs=(:ipe, :tex)) do cv
            draw_box!(cv, 0.0, 0.0, 100.0, 100.0; stroke=:darkred)
            draw_label!(cv, 50.0, 50.0, raw"x \in \mathcal{S}")
        end
        @test isfile(res.ipe)
        @test isfile(res.tex)

        ipe_only = figure(joinpath(temp_dir, "compact.ipe"); fit=unit_circle) do cv
            draw!(cv, unit_circle; stroke=:blue)
        end
        @test isfile(ipe_only.ipe)
        @test ipe_only.pdf === nothing
        @test ipe_only.tex === nothing

        @test_throws ArgumentError open_ipe(joinpath(temp_dir, "no_source");
                                             outputs=(:tex,), preview=true) do _
        end

    finally
        rm(temp_dir, recursive=true, force=true)
    end
end
