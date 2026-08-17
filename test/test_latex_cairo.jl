using Test
using BasicCompGeometry
using Cairo
using LaTeXStrings

@testset "LaTeX Cairo Extension (LaTeXCairoExt)" begin
    temp_dir = mktempdir()

    try
        # Test 1: Snippet extraction (IPE Prelim & custom markers)
        tex_sample = """
        \\documentclass{article}
        %%% IPE Prelim start
        \\newcommand{\\R}{\\mathbb{R}}
        \\DeclareMathOperator{\\Vor}{Vor}
        %%% IPE Prelim end
        \\begin{document}
        %%% MACROS BEG
        \\newcommand{\\eps}{\\varepsilon}
        %%% MACROS END
        Hello
        \\end{document}
        """
        tex_file = joinpath(temp_dir, "sample.tex")
        write(tex_file, tex_sample)

        # Auto-detect IPE prelims
        snip_ipe = read_latex_snippet(tex_file)
        @test occursin("\\newcommand{\\R}{\\mathbb{R}}", snip_ipe)
        @test occursin("\\DeclareMathOperator{\\Vor}{Vor}", snip_ipe)
        @test !occursin("\\newcommand{\\eps}{\\varepsilon}", snip_ipe)

        # Custom markers
        snip_custom = read_latex_snippet(tex_file; beg_marker="%%% MACROS BEG", end_marker="%%% MACROS END")
        @test occursin("\\newcommand{\\eps}{\\varepsilon}", snip_custom)
        @test !occursin("\\newcommand{\\R}{\\mathbb{R}}", snip_custom)

        # Test 2: Global preamble management
        reset_latex_preamble!()
        add_latex_packages!("bm", "xcolor")
        append_latex_preamble!(tex_file) # appends IPE prelims from file
        append_latex_preamble!(raw"\newcommand{\Prophet}{\mathcal{P}}")

        preamble_str = get_latex_preamble()
        @test occursin("\\usepackage{bm}", preamble_str)
        @test occursin("\\usepackage{xcolor}", preamble_str)
        @test occursin("\\newcommand{\\R}{\\mathbb{R}}", preamble_str)
        @test occursin("\\newcommand{\\Prophet}{\\mathcal{P}}", preamble_str)

        # Test 3: latex_to_pdf using global macros
        p1 = "\\section*{Test Page 1}\nThis uses macro: \$\\R^2\$ and \$\\Prophet\$."
        p2 = "\\section*{Test Page 2}\nAnother formula: \$\\int_0^\\infty e^{-x^2} dx = \\frac{\\sqrt{\\pi}}{2}\$."
        pdf_path = joinpath(temp_dir, "test_doc.pdf")

        latex_to_pdf([p1, p2], pdf_path; paperwidth=600, paperheight=600, margin=40)
        @test isfile(pdf_path)
        @test filesize(pdf_path) > 0

        # Test 4: cairo_draw_latex onto Cairo surface using macros
        cairo_pdf = joinpath(temp_dir, "cairo_math.pdf")
        surface = CairoPDFSurface(cairo_pdf, 500, 500)
        cr = CairoContext(surface)

        set_source_rgb(cr, 0.95, 0.95, 0.95)
        Cairo.rectangle(cr, 0, 0, 500, 500)
        Cairo.fill(cr)

        cairo_draw_latex(cr, 50, 50, L"\sum_{i=1}^n x_i \in \R"; fontsize=14, halign=:left)
        cairo_draw_latex(cr, 250, 250, L"\Vor(p) \subset \Prophet"; fontsize=16, halign=:center)
        cairo_draw_latex(cr, 450, 450, L"\alpha \le 2.0"; fontsize=12, halign=:right)

        Cairo.finish(surface)
        @test isfile(cairo_pdf)
        @test filesize(cairo_pdf) > 0

        # Test 5: pdf_merge
        merged_pdf = joinpath(temp_dir, "merged.pdf")
        pdf_merge(merged_pdf, (pdf_path, 1), (cairo_pdf, 1), (pdf_path, 2))
        @test isfile(merged_pdf)
        @test filesize(merged_pdf) > 0
    finally
        rm(temp_dir, recursive=true, force=true)
    end
end
