using Test
using BasicCompGeometry
using BasicCompGeometry.IpeDraw

@testset "IpeDraw Vector Figure Interface" begin
    @test isdefined(BasicCompGeometry, :IpeDraw)
    temp_dir = mktempdir()
    test_external = get(ENV, "BCG_TEST_EXTERNAL_TOOLS", "false") == "true"

    try
        # 1. Canvas creation and sizing
        canvas = IpeCanvas(width = 500.0, height = 300.0)
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
        draw_point!(canvas, p1; stroke = :darkred, fill = :red, shape = :disk)
        draw_point!(canvas, 100.0, 80.0; stroke = :blue, shape = :circle)
        draw_segment!(canvas, p1, p2; stroke = :black, pen = :heavier)
        draw_segment!(canvas, Segment(p1, p2); stroke = :gray4)

        bb = BBox(Point(20.0, 20.0), Point(250.0, 200.0))
        draw_box!(canvas, bb; stroke = :blue, fill = :lightblue)
        draw_box!(canvas, 10.0, 10.0, 50.0, 50.0)

        poly = [Point(10.0, 10.0), Point(50.0, 10.0), Point(30.0, 40.0)]
        draw_polygon!(canvas, poly; close = true, fill = :lightgreen)

        c = Circle(Point(100.0, 100.0), 30.0)
        draw_circle!(canvas, c; stroke = :darkgreen)
        draw_circle!(canvas, 80.0, 80.0, 15.0; stroke = :darkblue)

        draw_arc!(canvas, Point(100.0, 100.0), 25.0, 0.0, π / 2; stroke = :darkred)

        # Compact API, reusable styles, fitting, and page-space labels
        compact = IpeCanvas(width = 200.0, height = 100.0)
        unit_circle = Circle(point(0.0, 0.0), 1.0)
        fit!(compact, unit_circle; margin = 30.0)
        circle_style =
            Style(fill = :lightblue, fill_opacity = 0.2, stroke = :blue, pen = :heavier)
        draw!(compact, unit_circle; style = circle_style)
        mark!(compact, point(0.0, 0.0))
        label!(compact, point(0.0, 1.0), "p"; offset = (4.0, 5.0), anchor = :southwest)
        compact_xml = BasicCompGeometry.IpeDraw.to_xml(compact)
        @test compact.viewport.scale == 20.0
        @test occursin("20.000 0 0 20.000 100.000 50.000 e", compact_xml)
        @test occursin("opacity=\"20%\"", compact_xml)
        @test occursin("pos=\"104.000 75.000\"", compact_xml)

        layer(compact, :annotations) do cv
            mark!(cv, point(0.5, 0.0); stroke = :red)
        end
        @test compact.active_layer == "alpha"
        with_style(compact, Style(stroke = :darkred, pen = :heavier)) do cv
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
        inset(compact, detail_box; fit = unit_circle, margin = 5) do cv
            draw!(
                cv,
                unit_circle;
                style = theme.region,
                fill = :lightgreen,
                stroke = :darkgreen,
            )
        end
        @test compact.viewport === previous_viewport

        clipped_box = BBox(point(-0.5, -0.5), point(0.5, 0.5))
        clip_to(compact, clipped_box) do cv
            draw!(cv, unit_circle; stroke = :red)
        end
        legend!(
            compact,
            ["disk" => Style(fill = :lightblue, stroke = :blue)];
            at = (5, 55),
            width = 70,
        )
        scale_bar!(compact, 0.5; at = (10, 12), label = "0.5")
        furniture_xml = BasicCompGeometry.IpeDraw.to_xml(compact)
        @test length(findall(" clip=", furniture_xml)) == 2
        @test occursin("disk", furniture_xml)
        @test occursin("0.5", furniture_xml)
        @test compact.viewport === previous_viewport
        @test_throws ArgumentError legend!(compact, ["bad" => :blue])

        # 5. Conceptual & Algorithmic Helpers
        draw_bar!(
            canvas,
            50.0,
            450.0,
            120.0;
            label_left = "0",
            label_right = raw"M = \frac{1}{\eps}",
        )
        draw_span!(canvas, 150.0, 300.0, 120.0; fill = :lightgreen, stroke = :darkgreen)
        draw_dimension!(
            canvas,
            150.0,
            300.0,
            150.0;
            label = raw"\Delta \le 2\eps \mu",
            arrow = :both,
        )
        draw_arrow!(canvas, Point(100.0, 60.0), Point(150.0, 110.0))
        draw_curved_arrow!(canvas, Point(50.0, 200.0), Point(200.0, 200.0); bend = 25.0)

        draw_label!(
            canvas,
            250.0,
            50.0,
            raw"\mathcal{E}_{\le i} \implies |X - \mu| \le \eps \mu";
            size = :large,
            halign = :center,
        )

        # LaTeXStrings support test if available
        if Base.find_package("LaTeXStrings") !== nothing
            @eval begin
                using LaTeXStrings
                draw_label!(
                    $canvas,
                    250.0,
                    80.0,
                    L"\mathcal{E}_{\le i} \implies |X - \mu| \le \eps \mu",
                )
            end
        end

        # 5b. Ellipses, Arcs, Béziers, Splines, Holes, and Groups
        e = Ellipse(Point(300.0, 200.0), 40.0, 20.0, π / 6)
        draw_ellipse!(canvas, e; fill = :lightblue, stroke = :darkblue, tiling = :hatch)
        draw_ellipse!(canvas, Point(350.0, 200.0), 30.0, 15.0; stroke = :red)

        c_arc = CircleArc(Point(100.0, 100.0), 25.0, 0.0, π)
        draw_arc!(canvas, c_arc; stroke = :darkgreen)

        el_arc = EllipticArc(e, 0.0, π / 2)
        draw_arc!(canvas, el_arc; stroke = :magenta, arrow = :forward)

        bez = CubicBezier(
            Point(50.0, 50.0),
            Point(50.0, 100.0),
            Point(100.0, 100.0),
            Point(100.0, 50.0),
        )
        draw_bezier!(canvas, bez; stroke = :purple, pen = :heavier)

        knots = PntSeq([
            Point(200.0, 100.0),
            Point(220.0, 140.0),
            Point(260.0, 110.0),
            Point(300.0, 150.0),
        ])
        draw_spline!(canvas, knots; stroke = :darkblue, pen = :fat)
        draw_bspline!(canvas, knots; stroke = :gray3, dash = :dashed)

        outer_box = PntSeq([
            Point(50.0, 50.0),
            Point(150.0, 50.0),
            Point(150.0, 150.0),
            Point(50.0, 150.0),
        ])
        inner_hole = PntSeq([
            Point(75.0, 75.0),
            Point(125.0, 75.0),
            Point(125.0, 125.0),
            Point(75.0, 125.0),
        ])
        draw_polygon_with_holes!(
            canvas,
            outer_box,
            inner_hole;
            fill = :lightgray,
            stroke = :black,
        )

        ipe_group(
            canvas;
            matrix = [1.0, 0.0, 0.0, 1.0, 10.0, 20.0],
            opacity = Symbol("50%"),
        ) do cv
            draw_point!(cv, 0.0, 0.0; stroke = :red)
            draw_segment!(cv, Point(0.0, 0.0), Point(20.0, 20.0))
        end

        # 6. File serialization and PDF compilation
        base_fig = joinpath(temp_dir, "test_fig")
        artifacts = export_figure(
            canvas,
            base_fig;
            caption = "Test caption referencing the lemma.",
            label = "fig:test_sample",
            outputs = (:ipe, :tex),
        )

        @test isfile(artifacts.ipe)
        @test isfile(artifacts.tex)
        @test occursin("<ipe version=", read(artifacts.ipe, String))
        @test occursin(raw"\figlab{fig:test_sample}", read(artifacts.tex, String))

        if test_external
            pdf_artifacts = export_figure(canvas, base_fig; outputs = (:pdf,))
            @test isfile(pdf_artifacts.pdf)
            @test filesize(pdf_artifacts.pdf) > 0
        end

        # 7. Block syntax `open_ipe`
        res = open_ipe(
            joinpath(temp_dir, "block_fig");
            caption = "Block test",
            label = "fig:block",
            outputs = (:ipe, :tex),
        ) do cv
            draw_box!(cv, 0.0, 0.0, 100.0, 100.0; stroke = :darkred)
            draw_label!(cv, 50.0, 50.0, raw"x \in \mathcal{S}")
        end
        @test isfile(res.ipe)
        @test isfile(res.tex)

        ipe_only = figure(joinpath(temp_dir, "compact.ipe"); fit = unit_circle) do cv
            draw!(cv, unit_circle; stroke = :blue)
        end
        @test isfile(ipe_only.ipe)
        @test ipe_only.pdf === nothing
        @test ipe_only.tex === nothing

        @test_throws ArgumentError open_ipe(
            joinpath(temp_dir, "no_source");
            outputs = (:tex,),
            preview = true,
        ) do _
        end

    finally
        rm(temp_dir, recursive = true, force = true)
    end
end

@testset "Compact Ipe dispatch and missing tools" begin
    canvas = IpeCanvas(width = 200.0, height = 100.0)
    circle_arc = CircleArc(point(0.0, 0.0), 1.0, 0.0, π)
    ellipse = Ellipse(point(2.0, 2.0), 2.0, 1.0, 0.3)
    elliptic_arc = EllipticArc(ellipse, 0.0, π)
    bezier = CubicBezier(point(0.0, 0.0), point(0.0, 1.0), point(1.0, 1.0), point(1.0, 0.0))
    spline = CubicSpline([bezier])
    polygon_points = [point(0.0, 0.0), point(1.0, 0.0), point(0.0, 1.0)]
    objects = Any[
        BBox(point(0.0, 0.0), point(1.0, 1.0)),
        circle_arc,
        ellipse,
        elliptic_arc,
        bezier,
        spline,
        PntSeq(polygon_points),
        polygon_points,
    ]
    draw!(canvas, objects; stroke = :red)
    @test length(canvas.elements) == length(objects)
    @test all(element -> occursin("stroke=\"red\"", element), canvas.elements)
    draw_segment!(canvas, 0.0, 0.0, 1.0, 1.0)
    draw_spline!(canvas, PntSeq(polygon_points); method = :natural)
    @test length(canvas.elements) == length(objects) + 2
    draw_bspline!(canvas, polygon_points; closed = true)
    @test occursin(" u", canvas.elements[end]) || occursin("\nu", canvas.elements[end])
    theme = publication_theme()
    @test :region in propertynames(theme)
    @test :styles in propertynames(theme, true)
    @test_throws ArgumentError figure(_ -> nothing, "invalid.png")
    @test_throws ArgumentError figure(_ -> nothing, "invalid.pdf"; source = :other)

    mktempdir(mkpath(joinpath(@__DIR__, "..", "output", "tmp"))) do dir
        withenv("PATH" => dir) do
            @test_logs (:warn, r"ipetoipe command not found") @test !compile_pdf(
                canvas,
                joinpath(dir, "missing.pdf"),
            )
            @test_logs (:warn, r"ipetoipe command not found") @test !compile_pdf(
                "missing.ipe",
                joinpath(dir, "missing.pdf"),
            )
            @test_logs (:warn, r"ipe executable not found") @test edit_ipe(
                "missing.ipe",
            ) === nothing
            @test_logs (:warn, r"ipetoipe command not found") (
                :warn,
                r"PDF compilation failed",
            ) begin
                @test_throws ErrorException export_figure(
                    canvas,
                    joinpath(dir, "source");
                    outputs = (:ipe, :pdf),
                )
            end
            @test isfile(joinpath(dir, "source.ipe"))
        end
        if get(ENV, "BCG_TEST_EXTERNAL_TOOLS", "false") == "true"
            pdf = joinpath(dir, "canvas.pdf")
            @test compile_pdf(canvas, pdf)
            @test startswith(read(pdf, String), "%PDF-")
            artifact = figure(joinpath(dir, "no_extension"); keep_source = false) do cv
                draw!(cv, circle_arc)
            end
            @test isfile(artifact.pdf)
            @test artifact.ipe === nothing
        end
    end
end

@testset "IpeDraw Review 010 Performance Regressions" begin
    ID = BasicCompGeometry.IpeDraw

    # 1. Type-stable _apply_tf across coordinate types with and without viewport
    canvas = IpeCanvas(width = 200.0, height = 100.0)
    for pt in (Point(1, 2), Point(1.0f0, 2.0f0), Point(1 // 2, 3 // 4))
        tf_plain = @inferred ID._apply_tf(canvas, pt)
        @test tf_plain isa Point{2,Float64}
        @test tf_plain == Point(Float64(pt.x), Float64(pt.y))
    end
    fit!(canvas, BBox(Point(0.0, 0.0), Point(10.0, 10.0)))
    for pt in (Point(1, 2), Point(1.0f0, 2.0f0), Point(1 // 2, 3 // 4))
        tf_vp = @inferred ID._apply_tf(canvas, pt)
        @test tf_vp isa Point{2,Float64}
    end

    # 2. Parameterized Style, inferred _style_options, function barrier, and fill_opacity
    st = Style(stroke = :blue, fill = :red, fill_opacity = 0.2, extra_ignored = 42)
    @test st isa Style{@NamedTuple{stroke::Symbol,fill::Symbol,fill_opacity::Float64,extra_ignored::Int}}
    @test Style((stroke = :blue,)) == Style(stroke = :blue)

    opts = @inferred ID._style_options(
        Style(),
        st,
        (pen = :heavier,),
        Val(ID.PATH_STYLE),
    )
    @test opts === (stroke = :blue, fill = :red, fill_opacity = 0.2, pen = :heavier)
    @test_throws ArgumentError draw!(canvas, point(0.0, 0.0); bogus_option = 1)

    with_style(canvas, Style(stroke = :darkgreen, fill = :lightgreen, fill_opacity = 0.3)) do cv
        n_before = length(cv.elements)
        draw!(cv, BBox(point(0.0, 0.0), point(2.0, 2.0)); pen = :fat)
        @test length(cv.elements) == n_before + 2
        @test occursin("opacity=\"30%\"", cv.elements[end - 1])
        @test occursin("pen=\"fat\"", cv.elements[end])
    end
    @test isempty(canvas.active_style.values)

    segs = [Segment(point(Float64(i), 0.0), point(Float64(i), 1.0)) for i in 1:5]
    n_before = length(canvas.elements)
    draw!(canvas, segs; style = st)
    @test length(canvas.elements) == n_before + 5
    @test all(el -> occursin("stroke=\"blue\"", el), canvas.elements[(end - 4):end])

    # 3. Single-pass IOBuffer streaming for polygon, spline, B-spline, and holes
    empty!(canvas.elements)
    draw_polygon!(canvas, Point{2,Float64}[])
    draw_polygon!(canvas, [Point(1.0, 1.0)])
    draw_polygon!(canvas, (Point(Float64(i), 0.0) for i in 1:1))
    @test isempty(canvas.elements)

    draw_polygon!(canvas, [Point(0, 0), Point(1, 0), Point(1, 1)]; close = false)
    @test !occursin("\nh\n", canvas.elements[end])
    draw_polygon!(canvas, (Point(Float64(i), Float64(i % 2)) for i in 1:4); close = true)
    @test occursin("\nh\n", canvas.elements[end])

    empty_spline = CubicSpline(CubicBezier{2,Float64}[])
    n_before = length(canvas.elements)
    draw_spline!(canvas, empty_spline)
    @test length(canvas.elements) == n_before

    @test_throws ArgumentError draw_bspline!(canvas, [Point(0.0, 0.0), Point(1.0, 1.0)])
    draw_bspline!(canvas, (Point(Float64(i), Float64(i % 2)) for i in 1:4); closed = false)
    @test occursin(" m\n", canvas.elements[end]) && occursin("\ns\n", canvas.elements[end])
    draw_bspline!(canvas, (Point(0.0, 0.0), Point(1.0, 0.0), Point(0.0, 1.0)); closed = true)
    @test !occursin(" m\n", canvas.elements[end]) && occursin("\nu\n", canvas.elements[end])

    large_pts = [Point(Float64(i), Float64(i % 7)) for i in 1:2000]
    function _measure_path_allocs(cv, pts)
        draw_polygon!(cv, pts)
        poly_bytes = sizeof(cv.elements[end])
        empty!(cv.elements)
        poly_alloc = @allocated draw_polygon!(cv, pts)
        empty!(cv.elements)
        draw_bspline!(cv, pts)
        bspline_bytes = sizeof(cv.elements[end])
        empty!(cv.elements)
        bspline_alloc = @allocated draw_bspline!(cv, pts)
        empty!(cv.elements)
        return (poly_alloc, poly_bytes, bspline_alloc, bspline_bytes)
    end
    _measure_path_allocs(canvas, large_pts[1:10])
    poly_alloc, poly_bytes, bspline_alloc, bspline_bytes = _measure_path_allocs(canvas, large_pts)
    @test poly_alloc < 3 * poly_bytes
    @test bspline_alloc < 3 * bspline_bytes

    # 4. Label anchor lookup and single-pass XML escaping
    expected_anchors = (
        :center => ("center", "center"),
        :north => ("center", "top"),
        :south => ("center", "bottom"),
        :east => ("right", "center"),
        :west => ("left", "center"),
        :northeast => ("right", "top"),
        :northwest => ("left", "top"),
        :southeast => ("right", "bottom"),
        :southwest => ("left", "bottom"),
    )
    for (anc, (ha, va)) in expected_anchors
        draw_label!(canvas, 1.0, 2.0, "a < b & c > d"; anchor = anc)
        el = canvas.elements[end]
        @test occursin("halign=\"$ha\"", el)
        @test occursin("valign=\"$va\"", el)
        @test occursin("\$a &lt; b &amp; c &gt; d\$", el)
    end
    draw_label!(canvas, 1.0, 2.0, "x & y"; style = :normal)
    @test occursin(">x &amp; y</text>", canvas.elements[end])
    @test ID._escape_ipe_xml("a & b < c > d") == "a &amp; b &lt; c &gt; d"
    @test ID._escape_ipe_xml("plain text") == "plain text"
    @test_throws ArgumentError draw_label!(canvas, 0.0, 0.0, "x"; anchor = :invalid_anchor)
end
