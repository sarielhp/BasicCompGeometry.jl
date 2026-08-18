using Test
using BasicCompGeometry
using Cairo

@testset "Unified Canvas (PDF, SVG, PNG, GIF)" begin
    temp_dir = mktempdir()

    try
        pts = [Point(1.0, 1.0), Point(5.0, 2.0), Point(4.0, 5.0), Point(2.0, 4.0)]
        bb = BBox(Point(0.0, 0.0), Point(6.0, 6.0))

        # Test 1: Multi-page PDF
        pdf_file = joinpath(temp_dir, "test.pdf")
        open_canvas(pdf_file, 400, 400) do canvas
            cairo_draw_setup(canvas, bb, 400, 400)
            cairo_draw_points(canvas, pts, 3.0)
            Cairo.show_page(canvas)

            cairo_draw_setup(canvas, bb, 400, 400)
            cairo_draw_polygon(canvas, pts)
            Cairo.show_page(canvas)
        end
        @test isfile(pdf_file)
        @test filesize(pdf_file) > 0

        # Test 2: Single-page SVG (automatic renaming from _001 to exact filename)
        svg_single = joinpath(temp_dir, "single.svg")
        open_canvas(svg_single, 400, 400) do canvas
            cairo_draw_setup(canvas, bb, 400, 400)
            cairo_draw_points(canvas, pts, 3.0)
        end
        @test isfile(svg_single)
        @test filesize(svg_single) > 0

        # Test 3: Multi-page SVG sequence
        svg_multi = joinpath(temp_dir, "deck.svg")
        open_canvas(svg_multi, 400, 400) do canvas
            # Page 1
            cairo_draw_setup(canvas, bb, 400, 400)
            cairo_draw_points(canvas, pts, 3.0)
            Cairo.show_page(canvas)

            # Page 2
            cairo_draw_setup(canvas, bb, 400, 400)
            cairo_draw_polygon(canvas, pts)
            Cairo.show_page(canvas)
        end
        p1 = joinpath(temp_dir, "deck_001.svg")
        p2 = joinpath(temp_dir, "deck_002.svg")
        @test isfile(p1) && filesize(p1) > 0
        @test isfile(p2) && filesize(p2) > 0

        # Test 4: Formatted pattern SVG (e.g. slide_%02d.svg)
        svg_pat = joinpath(temp_dir, "slide_%02d.svg")
        open_canvas(svg_pat, 400, 400) do canvas
            cairo_draw_setup(canvas, bb, 400, 400)
            cairo_draw_points(canvas, pts, 3.0)
            Cairo.show_page(canvas)

            cairo_draw_setup(canvas, bb, 400, 400)
            cairo_draw_polygon(canvas, pts)
            Cairo.show_page(canvas)
        end
        s1 = joinpath(temp_dir, "slide_01.svg")
        s2 = joinpath(temp_dir, "slide_02.svg")
        @test isfile(s1) && filesize(s1) > 0
        @test isfile(s2) && filesize(s2) > 0

        # Test 5: PNG Sequence
        png_pat = joinpath(temp_dir, "frame_%03d.png")
        open_canvas(png_pat, 200, 200) do canvas
            for _ in 1:3
                cairo_draw_setup(canvas, bb, 200, 200)
                cairo_draw_points(canvas, pts, 2.0)
                Cairo.show_page(canvas)
            end
        end
        f1 = joinpath(temp_dir, "frame_001.png")
        f2 = joinpath(temp_dir, "frame_002.png")
        f3 = joinpath(temp_dir, "frame_003.png")
        @test isfile(f1) && filesize(f1) > 0
        @test isfile(f2) && filesize(f2) > 0
        @test isfile(f3) && filesize(f3) > 0

        # Test 6: Animated GIF (if ffmpeg is available)
        if Sys.which("ffmpeg") !== nothing
            gif_file = joinpath(temp_dir, "anim.gif")
            open_canvas(gif_file, 200, 200; fps=5) do canvas
                for _ in 1:4
                    cairo_draw_setup(canvas, bb, 200, 200)
                    cairo_draw_points(canvas, pts, 2.0)
                    Cairo.show_page(canvas)
                end
            end
            @test isfile(gif_file)
            @test filesize(gif_file) > 0
            @test get_file_path(gif_file) == abspath(gif_file)
        end

        # Test 7: HTML Presentation Canvas
        html_out = joinpath(temp_dir, "my_slides.html")
        c_res = open_canvas(html_out, 400, 400) do canvas
            # Page 1
            cairo_draw_setup(canvas, bb, 400, 400)
            cairo_draw_points(canvas, pts, 3.0)
            description(canvas, "Initial point configuration.")
            Cairo.show_page(canvas)

            # Page 2
            cairo_draw_setup(canvas, bb, 400, 400)
            cairo_draw_polygon(canvas, pts)
            description(canvas, "Convex hull / polygon boundary.")
            Cairo.show_page(canvas)
        end

        expected_html = joinpath(temp_dir, "my_slides", "index.html")
        @test get_file_path(c_res) == abspath(expected_html)
        @test isfile(expected_html)
        @test isfile(joinpath(temp_dir, "my_slides", "page_001.svg"))
        @test isfile(joinpath(temp_dir, "my_slides", "page_002.svg"))

        html_txt = read(expected_html, String)
        @test occursin("my_slides", html_txt)
        @test occursin("Initial point configuration.", html_txt)
        @test occursin("Convex hull / polygon boundary.", html_txt)

    finally
        rm(temp_dir, recursive=true, force=true)
    end
end
