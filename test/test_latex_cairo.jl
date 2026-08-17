using Test
using BasicCompGeometry
using Cairo
using LaTeXStrings

@testset "LaTeX Cairo Extension (LaTeXCairoExt)" begin
    temp_dir = mktempdir()

    try
        # Test 1: latex_to_pdf
        p1 = "\\section*{Test Page 1}\nThis is a mathematical formula: \$E = mc^2\$."
        p2 = "\\section*{Test Page 2}\nAnother formula: \$\\int_0^\\infty e^{-x^2} dx = \\frac{\\sqrt{\\pi}}{2}\$."
        pdf_path = joinpath(temp_dir, "test_doc.pdf")

        latex_to_pdf([p1, p2], pdf_path; paperwidth=600, paperheight=600, margin=40)
        @test isfile(pdf_path)
        @test filesize(pdf_path) > 0

        # Test 2: cairo_draw_latex onto Cairo surface
        cairo_pdf = joinpath(temp_dir, "cairo_math.pdf")
        surface = CairoPDFSurface(cairo_pdf, 500, 500)
        cr = CairoContext(surface)

        set_source_rgb(cr, 0.95, 0.95, 0.95)
        Cairo.rectangle(cr, 0, 0, 500, 500)
        Cairo.fill(cr)

        cairo_draw_latex(cr, 50, 50, L"\sum_{i=1}^n x_i = \Theta(n)"; fontsize=14, halign=:left)
        cairo_draw_latex(cr, 250, 250, L"\text{Prophet Value } p^*"; fontsize=16, halign=:center)
        cairo_draw_latex(cr, 450, 450, L"\alpha \le 2.0"; fontsize=12, halign=:right)

        Cairo.finish(surface)
        @test isfile(cairo_pdf)
        @test filesize(cairo_pdf) > 0

        # Test 3: pdf_merge
        merged_pdf = joinpath(temp_dir, "merged.pdf")
        pdf_merge(merged_pdf, (pdf_path, 1), (cairo_pdf, 1), (pdf_path, 2))
        @test isfile(merged_pdf)
        @test filesize(merged_pdf) > 0
    finally
        rm(temp_dir, recursive=true, force=true)
    end
end
