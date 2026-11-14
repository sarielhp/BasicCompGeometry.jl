# Golden-file (snapshot) tests for the package's file writers and readers.
#
# Each test builds an output deterministically and compares it byte for byte
# with a committed fixture in `test/golden/`. Any change in output format fails
# here until the fixture is deliberately regenerated and the diff reviewed:
#
#     BCG_UPDATE_GOLDEN=1 julia --project=. -e 'using Test; include("test/test_golden.jl")'
#
# See test/golden/README.md.

using Test
using BasicCompGeometry
using BasicCompGeometry.IpeDraw
using BasicCompGeometry.ReadWrite

const GOLDEN_DIR = joinpath(@__DIR__, "golden")
const GOLDEN_TMP = joinpath(@__DIR__, "..", "output", "tmp", "golden_test")
const GOLDEN_UPDATE = get(ENV, "BCG_UPDATE_GOLDEN", "") == "1"

# Describe the first differing line of two byte strings, or `nothing` if equal.
function golden_first_difference(expected::Vector{UInt8}, actual::Vector{UInt8})
    expected == actual && return nothing
    exp_lines = split(String(copy(expected)), '\n')
    act_lines = split(String(copy(actual)), '\n')
    for i in 1:max(length(exp_lines), length(act_lines))
        e = i <= length(exp_lines) ? exp_lines[i] : nothing
        a = i <= length(act_lines) ? act_lines[i] : nothing
        e == a && continue
        show_line(x) = x === nothing ? "<end of file>" : repr(String(x))
        return "first difference at line $i:\n  expected: $(show_line(e))\n  actual:   $(show_line(a))\n" *
               "  (expected $(length(expected)) bytes, got $(length(actual)) bytes)"
    end
    # Same lines but different bytes (cannot happen for valid splits; be safe).
    k = findfirst(i -> i > length(actual) || expected[i] != actual[i], 1:length(expected))
    return "first differing byte at offset $(something(k, length(expected) + 1))"
end

# Compare `actual` with fixture `name`, or write it when updating.
function check_golden(name::AbstractString, actual::Vector{UInt8})
    path = joinpath(GOLDEN_DIR, name)
    if GOLDEN_UPDATE
        mkpath(GOLDEN_DIR)
        write(path, actual)
        @info "Golden fixture regenerated; review the diff before committing" path
        @test_broken false  # nothing was verified
        return
    end
    if !isfile(path)
        @error "Missing golden fixture $name; regenerate with BCG_UPDATE_GOLDEN=1"
        @test isfile(path)
        return
    end
    diff = golden_first_difference(read(path), actual)
    if diff !== nothing
        @error "Golden mismatch for $name (if intended, regenerate with BCG_UPDATE_GOLDEN=1 and review the diff): $diff"
    end
    @test diff === nothing
end

check_golden_file(name::AbstractString, file::AbstractString) =
    check_golden(name, read(file))

tmp_path(name) = joinpath(GOLDEN_TMP, name)

# Coordinates that stress fixed-point formatting: binary values just below or
# above decimal ties, signed zeros, tiny negatives, and large magnitudes.
const TIE_VALUES = [
    0.0025,
    0.1235,
    -0.0025,
    -0.1235,
    -1e-17,
    1e-17,
    -0.0,
    0.0005,
    1.0005,
    2.5e-4,
    0.0015,
    -0.0045,
    123456.7895,
    1.0e12 + 0.25,
    -9876543.2105,
    3.14159265358979,
]

tie_points() = [Point(TIE_VALUES[i], TIE_VALUES[end+1-i]) for i in eachindex(TIE_VALUES)]

ipe_bytes(canvas) = (save_ipe(canvas, tmp_path("canvas.ipe")); read(tmp_path("canvas.ipe")))

function build_primitives_canvas()
    cv = IpeCanvas(width = 400.0, height = 300.0)
    # Points / marks in all shapes and sizes, at tie coordinates
    for (i, p) in enumerate(tie_points())
        shape = (:disk, :circle, :box)[mod1(i, 3)]
        draw_point!(cv, p; shape = shape, stroke = :red, fill = :blue, size = :large)
    end
    draw_point!(cv, 0.0025, -1e-17)
    draw_points!(cv, [Point(1.0, 2.0), Point(0.1235, 0.0015)]; shape = :circle)
    # Segments
    draw_segment!(cv, Point(0.0025, 0.1235), Point(-1e-17, 1.0e12 + 0.25))
    draw_segment!(
        cv,
        Segment(Point(1.0, 1.0), Point(2.0, 3.0));
        arrow = :pointed,
        rarrow = :normal,
        dash = :dashed,
        pen = :heavier,
        stroke = :darkred,
    )
    draw_segment!(cv, 0.0005, 1.0005, 2.5e-4, -0.0045; arrow = :fancy)
    # Boxes
    draw_box!(
        cv,
        BBox(Point(-0.0025, 0.0025), Point(10.1235, 20.0015));
        fill = :lightblue,
        opacity = Symbol("50%"),
    )
    draw_box!(cv, 5.0, 5.0, -1e-17, 0.0005; dash = :dotted, pen = :fat)
    # Circles
    draw_circle!(cv, Point(0.0025, 0.1235), 0.0015; fill = :red, opacity = "30%")
    draw_circle!(cv, Sphere(Point(100.0, 100.0), 50.0))
    draw_circle!(cv, -1e-17, 1.0005, 2.5e-4)
    # Circular arcs
    draw_arc!(cv, Point(50.0, 50.0), 20.0, 0.0, pi / 2; arrow = :pointed)
    draw_arc!(cv, CircleArc(Point(0.0025, -0.0025), 1.0005, -pi / 3, 5pi / 4))
    # Ellipses and elliptic arcs
    draw_ellipse!(
        cv,
        Point(10.0, 20.0),
        30.0,
        10.0;
        angle = pi / 6,
        fill = :gray7,
        tiling = :hatch,
    )
    draw_ellipse!(cv, Ellipse(Point(0.1235, 0.0025), 2.0, 1.0, 0.0))
    draw_elliptic_arc!(cv, Point(5.0, 5.0), 4.0, 2.0, pi / 4, 0.0, pi; rarrow = :pointed)
    draw_arc!(cv, EllipticArc(Ellipse(Point(-1e-17, 0.0015), 3.0, 1.5, -pi / 8), 0.5, 2.5))
    # Bezier and splines
    pts = tie_points()[1:6]
    draw_bezier!(
        cv,
        CubicBezier(pts[1], pts[2], pts[3], pts[4]);
        arrow = :pointed,
        fill = :none,
        dash = :dashed,
    )
    seq = PntSeq([
        Point(0.0025, 0.0),
        Point(1.1235, 2.0015),
        Point(3.0, -1e-17),
        Point(4.0005, 1.0),
        Point(6.0, 0.1235),
    ])
    draw_spline!(cv, seq)
    draw_spline!(cv, seq; closed = true, stroke = :blue)
    draw_spline!(cv, seq; method = :natural)
    draw_spline!(cv, seq; method = :natural, closed = true, rarrow = :pointed)
    draw_spline!(
        cv,
        CubicSpline(
            [
                CubicBezier(pts[1], pts[2], pts[3], pts[4]),
                CubicBezier(pts[4], pts[5], pts[6], pts[1]),
            ];
            closed = true,
        );
        fill = :lightred,
    )
    # B-splines, open and closed
    draw_bspline!(cv, seq)
    draw_bspline!(cv, seq; closed = true, fill = :lightgreen, pen = :thin)
    draw_bspline!(cv, collect(pts[1:3]))
    # Polygons: closed, open, vectors, iterators with unknown length
    draw_polygon!(cv, seq; fill = :yellow, opacity = Symbol("70%"))
    draw_polygon!(cv, tie_points(); close = false, dash = :dashed)
    draw_polygon!(cv, Iterators.filter(p -> true, pts); stroke = :green)
    draw_polygon!(cv, [Point(1.0, 1.0)])  # too short: emits nothing
    # Polygons with holes
    outer = PntSeq([
        Point(0.0, 0.0),
        Point(10.0025, 0.0),
        Point(10.0025, 10.1235),
        Point(-1e-17, 10.1235),
    ])
    hole1 = PntSeq([Point(1.0005, 1.0005), Point(2.0, 1.0005), Point(2.0, 2.0015)])
    hole2 = PntSeq([Point(5.0, 5.0), Point(6.0, 5.0), Point(6.0, 6.0), Point(5.0, 6.0)])
    draw_polygon_with_holes!(
        cv,
        outer,
        hole1,
        hole2;
        opacity = Symbol("50%"),
        tiling = :crosshatch,
    )
    draw_polygon_with_holes!(cv, outer; fill = :lightblue)
    return cv
end

function build_viewport_canvas(flip_y::Bool)
    cv = IpeCanvas(width = 300.0, height = 200.0, layer = "geo")
    world = PntSeq([Point(-1.0, -0.5), Point(2.0, 1.5), Point(0.25, 0.75)])
    fit!(cv, world; margin = 17.0, flip_y = flip_y)
    draw_polygon!(cv, world; fill = :lightgray)
    draw_circle!(cv, Point(0.5, 0.5), 0.3)
    draw_arc!(cv, Point(0.5, 0.5), 0.5, 0.25, 2.0)
    draw_ellipse!(cv, Point(0.0, 0.0), 0.8, 0.3; angle = 0.7)
    draw_elliptic_arc!(cv, Point(1.0, 1.0), 0.6, 0.2, -0.3, 1.0, 4.0)
    draw_bspline!(cv, world; closed = true)
    draw_spline!(cv, world)
    draw_box!(cv, BBox(Point(-1.0, -0.5), Point(0.0025, 0.1235)))
    draw_point!(cv, Point(0.1235, -0.0025))
    draw_label!(cv, Point(0.0025, 0.1235), "p")
    return cv
end

function build_labels_styles_canvas()
    cv = IpeCanvas(width = 320.0, height = 240.0, layer = "base")
    add_preamble!(cv, raw"\newcommand{\eps}{\varepsilon}")
    add_layer!(cv, "overlay")
    add_view!(cv, ["base"])
    add_view!(cv, ["base", "overlay"]; active = "overlay")
    # Labels: math wrapping, escaping, alignment, anchors, offsets, sizes
    draw_label!(cv, Point(0.0025, 0.1235), raw"\eps < 1 & x > 0")
    draw_label!(cv, -1e-17, 1.0005, raw"$a_i$ already math"; halign = :left)
    draw_label!(
        cv,
        Point(10.0, 10.0),
        "plain <text> & more";
        style = :normal,
        size = :small,
        stroke = :red,
        valign = :top,
    )
    draw_label!(cv, Point(20.0, 20.0), raw"\begin{array}{c} x \end{array}")
    draw_label!(cv, Point(30.0, 30.0), 42; offset = (0.0025, -0.0025), anchor = :northeast)
    for anchor in
        (:center, :north, :south, :east, :west, :northwest, :southwest, :southeast)
        draw_label!(cv, Point(40.0, 40.0), string(anchor); anchor = anchor, style = :normal)
    end
    label!(cv, Point(50.0, 50.0), "q_1"; offset = (3, 4))
    label!(cv, 60.0, 60.0, "r"; appearance = Style(stroke = :blue, size = :large))
    # Styles, with_style, themes
    theme = publication_theme()
    custom = Theme(hot = Style(stroke = :red, fill = :yellow), cold = (stroke = :blue,))
    box = BBox(Point(0.0, 0.0), Point(10.0025, 5.1235))
    draw!(cv, box; style = theme.region, fill = :blue)
    draw!(cv, box; style = custom.hot, fill_opacity = 0.3)
    draw!(cv, box; style = custom[:cold], fill = :green, fill_opacity = 1)
    draw!(cv, box; fill = :red, fill_opacity = 0)
    draw!(
        cv,
        Sphere(Point(1.0, 2.0), 0.5);
        style = theme.boundary,
        fill = :gray5,
        fill_opacity = Symbol("40%"),
    )
    draw!(cv, Point(3.0, 4.0); style = theme.point)
    draw!(cv, Segment(Point(0.0, 0.0), Point(1.0, 1.0)); arrow = :pointed)
    draw!(cv, CircleArc(Point(0.0, 0.0), 2.0, 0.0, 1.0))
    draw!(cv, Ellipse(Point(1.0, 1.0), 2.0, 1.0, 0.5); fill = :red, fill_opacity = 0.5)
    draw!(cv, EllipticArc(Ellipse(Point(1.0, 1.0), 2.0, 1.0, 0.5), 0.0, 1.0))
    draw!(
        cv,
        CubicBezier(Point(0.0, 0.0), Point(1.0, 2.0), Point(3.0, 2.0), Point(4.0, 0.0)),
    )
    draw!(
        cv,
        interpolate_catmull_rom(
            PntSeq([Point(0.0, 0.0), Point(1.0, 1.0), Point(2.0, 0.0)]),
        ),
    )
    draw!(cv, PntSeq([Point(0.0, 0.0), Point(1.0, 0.0), Point(1.0, 1.0)]); close = false)
    draw!(
        cv,
        [Point(0.0, 0.0), Point(2.0, 0.0), Point(2.0, 2.0)];
        fill = :lightblue,
        fill_opacity = 0.2,
    )
    draw!(
        cv,
        Any[Point(5.0, 5.0), Segment(Point(5.0, 5.0), Point(6.0, 6.0)), box];
        stroke = :darkblue,
    )
    mark!(cv, Point(7.0, 7.0); shape = :box)
    mark!(cv, [Point(8.0, 8.0), Point(9.0, 9.0)]; style = theme.point, size = :tiny)
    with_style(cv, Style(stroke = :darkgreen, pen = :fat)) do c
        draw!(c, box)
        draw!(c, Segment(Point(0.0, 0.0), Point(3.0, 0.0)))
        with_style(c; dash = :dashed) do c2
            draw!(c2, box; fill = :lightgreen)
            label!(c2, Point(1.0, 1.0), "styled"; stroke = :red)
        end
    end
    draw!(cv, box)  # style restored
    layer(cv, "overlay") do c
        draw!(c, box; stroke = :red)
        label!(c, Point(2.0, 2.0), "on overlay")
    end
    layer(cv, :extra) do c
        mark!(c, Point(0.0025, 0.0025))
    end
    draw!(cv, Point(0.0, 0.0))  # back on "base"
    return cv
end

function build_diagrams_canvas()
    cv = IpeCanvas(width = 400.0, height = 300.0, paper = "400 300", layer = "main")
    world = PntSeq([Point(0.0, 0.0), Point(4.0, 0.0), Point(4.0, 3.0), Point(0.0, 3.0)])
    fit!(cv, world; margin = 25.0)
    draw!(cv, world; fill = :lightgray)
    # Groups: matrix, opacity, clip
    ipe_group(
        cv;
        matrix = [1.0, 0.0, 0.0, 1.0, 0.0025, -0.1235],
        opacity = Symbol("50%"),
    ) do c
        draw!(c, Point(1.0, 1.0))
        ipe_group(c) do c2
            draw!(c2, Segment(Point(0.0, 0.0), Point(1.0, 1.0)))
        end
    end
    ipe_group(cv; clip = "0 0 m 10 0 l 10 10 l h") do c
        draw!(c, Point(2.0, 2.0))
    end
    clip_to(cv, BBox(Point(1.0, 1.0), Point(3.0, 2.0))) do c
        draw!(c, Sphere(Point(2.0, 1.5), 1.0); fill = :red)
    end
    clip_to(cv, PageBox(10.0, 10.0, 50.0025, 40.1235)) do c
        draw!(c, Point(0.5, 0.5))
    end
    # Page space and insets
    page_space(cv) do c
        draw!(c, Point(5.0, 5.0))
        draw!(c, BBox(Point(0.0, 0.0), Point(400.0, 300.0)); stroke = :gray5)
    end
    inset(
        cv,
        PageBox(280.0, 200.0, 100.0, 80.0);
        fit = PntSeq([Point(-1.0, -1.0), Point(1.0, 1.0)]),
    ) do c
        draw!(c, Sphere(Point(0.0, 0.0), 1.0))
        label!(c, Point(0.0, 0.0), "inset")
    end
    inset(
        cv,
        PageBox(20.0, 200.0, 60.0, 60.0);
        fit = (BBox(Point(0.0, 0.0), Point(2.0, 1.0)),),
        margin = 3,
        flip_y = true,
        clip = false,
    ) do c
        draw!(c, Segment(Point(0.0, 0.0), Point(2.0, 1.0)))
        draw!(c, CircleArc(Point(1.0, 0.5), 0.4, 0.0, 2.0))
    end
    draw!(cv, Point(4.0, 3.0))  # world viewport restored after inset
    # Legend and scale bar
    theme = publication_theme()
    legend!(
        cv,
        [
            "region" => theme.region,
            raw"a & <b>" => Style(fill = :red),
            "line" => Style(stroke = :blue, pen = :fat),
        ],
    )
    legend!(cv, ["sw" => Style(fill = :green)]; position = :southwest, width = 80)
    legend!(cv, ["at" => Style(fill = :yellow)]; at = (100.0025, 50.1235))
    scale_bar!(cv, 1.0)
    scale_bar!(cv, 0.5; position = :southeast, label = raw"$0.5$")
    scale_bar!(cv, 2.0; at = (150.0, 20.0), label = nothing)
    # Bars, spans, dimensions, arrows
    page_space(cv) do c
        draw_bar!(c, 10.0, 110.0, 150.0; label_left = 0, label_right = "n")
        draw_bar!(c, 10.0, 110.0, 120.0; height = 10.0025)
        draw_span!(c, 30.0, 60.0, 150.0)
        draw_span!(c, 70.0, 80.0, 150.0; dash = :dashed, fill = :lightred)
        draw_dimension!(c, Point(10.0, 100.0), Point(110.0, 130.0); label = raw"\ell")
        draw_dimension!(c, 10.0, 110.0, 90.0; label = "100", arrow = :pointed)
        draw_dimension!(c, Point(10.0, 80.0), Point(10.0, 80.0))
        draw_arrow!(c, Point(120.0, 150.0), Point(200.0, 150.0))
        draw_arrow!(c, 120.0, 140.0, 200.0, 100.0; dash = :solid, arrow = :normal)
        draw_curved_arrow!(c, Point(120.0, 90.0), Point(160.0, 120.0), Point(200.0, 90.0))
        draw_curved_arrow!(c, Point(120.0, 80.0), Point(200.0, 80.0); bend = -15.0)
        draw_curved_arrow!(c, Point(0.0025, 0.0025), Point(0.0025, 0.0025))
    end
    return cv
end

const FIG_TEX_BASE = "golden_fig"

@testset "Golden outputs" begin
    mkpath(GOLDEN_TMP)
    GOLDEN_UPDATE &&
        @warn "BCG_UPDATE_GOLDEN=1: regenerating golden fixtures instead of comparing"

    @testset "Ipe primitives" begin
        check_golden("ipe_primitives.ipe", ipe_bytes(build_primitives_canvas()))
    end
    @testset "Ipe viewport" begin
        check_golden("ipe_viewport.ipe", ipe_bytes(build_viewport_canvas(false)))
        check_golden("ipe_viewport_flip.ipe", ipe_bytes(build_viewport_canvas(true)))
    end
    @testset "Ipe labels and styles" begin
        check_golden("ipe_labels_styles.ipe", ipe_bytes(build_labels_styles_canvas()))
    end
    @testset "Ipe diagrams" begin
        check_golden("ipe_diagrams.ipe", ipe_bytes(build_diagrams_canvas()))
    end
    @testset "Ipe figure LaTeX fragment" begin
        tex = save_figure_tex(
            tmp_path(FIG_TEX_BASE * "_fig.tex");
            caption = raw"A $\eps$ caption",
        )
        check_golden_file("figure_fig.tex", tex)
        base = tmp_path("exported")
        artifacts = export_figure(
            IpeCanvas(width = 50.0, height = 40.0),
            base;
            caption = "Exported",
            outputs = (:tex,),
        )
        check_golden_file("export_fig.tex", artifacts.tex)
    end

    @testset "ReadWrite binary writers" begin
        fm = Float32[0.0025 -1.5 1.0f-30; 3.25 -0.0 Inf32; 1.0f6 2.0 -7.125]
        im = Int32[1 -2 typemax(Int32); typemin(Int32) 0 42]
        write_fbin(tmp_path("m.fbin"), fm)
        write_ibin(tmp_path("m.ibin"), im)
        write_fbin(tmp_path("empty.fbin"), Matrix{Float32}(undef, 0, 3))
        check_golden_file("matrix.fbin", tmp_path("m.fbin"))
        check_golden_file("matrix.ibin", tmp_path("m.ibin"))
        check_golden_file("empty.fbin", tmp_path("empty.fbin"))
        # The golden files must also read back to the same values.
        if !GOLDEN_UPDATE
            @test read_fbin_n(joinpath(GOLDEN_DIR, "matrix.fbin")) == permutedims(fm)
            @test read_ibin(joinpath(GOLDEN_DIR, "matrix.ibin")) == im
        end
    end

    @testset "ReadWrite text readers" begin
        # Canonical dump of the parsed points; repr gives exact (round-trip) floats.
        io = IOBuffer()
        dump_seq(label, P) = begin
            println(io, label, " n=", length(P))
            for p in P
                println(io, join((repr(Float64(c)) for c in p), ' '))
            end
        end
        dump_seq("read_plt dim=2", read_plt(joinpath(GOLDEN_DIR, "input_points.plt")))
        dump_seq(
            "read_plt dim=3",
            read_plt(joinpath(GOLDEN_DIR, "input_points.plt"); dim = 3),
        )
        dump_seq(
            "read_geolife_plt",
            read_geolife_plt(joinpath(GOLDEN_DIR, "input_geolife.plt")),
        )
        fv = read_fvecs(joinpath(GOLDEN_DIR, "input.fvecs"))
        println(io, "read_fvecs size=", size(fv))
        for j in axes(fv, 2)
            println(io, join((repr(x) for x in fv[:, j]), ' '))
        end
        check_golden("readers_dump.txt", take!(io))
    end

    @testset "Point-sequence and halfplane text writers" begin
        P = PntSeq([
            Point(0.0025, -1e-17),
            Point(0.1235, 1.0e12 + 0.25),
            Point(-0.0, 1 / 3),
        ])
        BasicCompGeometry.write_plt(tmp_path("seq.plt"), P)
        check_golden_file("write_plt.plt", tmp_path("seq.plt"))
        P3 = PntSeq([Point(1.0, 2.0, 3.0), Point(-4.5, 1e-300, 6.0e22)])
        BasicCompGeometry.write_plt(tmp_path("seq3.plt"), P3)
        check_golden_file("write_plt_3d.plt", tmp_path("seq3.plt"))
        hps = [
            Halfplane(Line(Point(0.0025, -1e-17), Point(1.0, 0.0))),
            Halfplane(Line(Point(1 / 3, 2 / 3), Point(-0.6, 0.8))),
            Halfplane(Line(Point(1.0e12 + 0.25, -0.0), Point(0.0, -1.0))),
        ]
        BasicCompGeometry.write_halfplanes(tmp_path("hps.txt"), hps)
        check_golden_file("halfplanes.txt", tmp_path("hps.txt"))
        if !GOLDEN_UPDATE
            back = BasicCompGeometry.read_halfplanes(joinpath(GOLDEN_DIR, "halfplanes.txt"))
            @test [(h.boundary.p, h.boundary.u) for h in back] == [(h.boundary.p, h.boundary.u) for h in hps]
        end
    end

    rm(GOLDEN_TMP; recursive = true, force = true)
end
