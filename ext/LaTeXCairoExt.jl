module LaTeXCairoExt

using BasicCompGeometry
import Cairo
import LaTeXStrings: LaTeXString

const LATEX_SNIPPET_CACHE = Dict{Tuple{String, Float64, Int, String}, Cairo.CairoSurfaceBase{UInt32}}()

"""
    BasicCompGeometry.cairo_draw_latex(cr, x, y, latex_str;
                                       fontsize=12.0,
                                       halign=:left,
                                       valign=:baseline,
                                       dpi=300,
                                       scale=1.0,
                                       compiler="lualatex",
                                       packages=["amsmath", "amssymb", "amsfonts"])

Render a LaTeX snippet or mathematical formula `latex_str` onto the Cairo context `cr` at position `(x, y)`.

# Arguments
- `cr`: A `CairoContext`.
- `x::Real`, `y::Real`: Insertion coordinate.
- `latex_str::Union{AbstractString, LaTeXString}`: The LaTeX formula or text string (e.g. `L"\\sum_{i=1}^n x_i"` or `"\\textbf{Site } p_k"`).
- `fontsize::Real`: Base font size in points (default: `12.0`).
- `halign::Symbol`: Horizontal alignment: `:left` (default), `:center`, or `:right`.
- `valign::Symbol`: Vertical alignment: `:baseline` (default), `:top`, `:center`, or `:bottom`.
- `dpi::Int`: Rasterization DPI (default: `300`).
- `scale::Real`: Additional scaling factor (default: `1.0`).
- `compiler::String`: TeX engine to use (default: `"lualatex"`).
- `packages::Vector{String}`: Additional LaTeX packages to include.
"""
function BasicCompGeometry.cairo_draw_latex(
    cr,
    x::Real,
    y::Real,
    latex_str::Union{AbstractString, LaTeXString};
    fontsize::Real = 12.0,
    halign::Symbol = :left,
    valign::Symbol = :baseline,
    dpi::Int = 300,
    scale::Real = 1.0,
    compiler::String = "lualatex",
    packages::Vector{String} = ["amsmath", "amssymb", "amsfonts"]
)
    str = String(latex_str)
    cache_key = (str, Float64(fontsize), dpi, compiler)

    img = get!(LATEX_SNIPPET_CACHE, cache_key) do
        temp_dir = mktempdir()
        tex_path = joinpath(temp_dir, "snippet.tex")
        pdf_path = joinpath(temp_dir, "snippet.pdf")
        png_base = joinpath(temp_dir, "snippet")
        png_path = joinpath(temp_dir, "snippet.png")

        pkg_lines = join(["\\usepackage{$p}" for p in packages], "\n")

        # If math mode not explicitly enclosed, wrap math if contains math symbols
        content = strip(str)
        is_math = startswith(content, "\$") || startswith(content, "\\[") || startswith(content, "\\begin{equation")
        body = is_math ? content : "\\textbf{\\textnormal{$content}}"

        tex_doc = """
\\documentclass[border=0pt,varwidth,12pt]{standalone}
$pkg_lines
\\begin{document}
$body
\\end{document}
"""
        open(tex_path, "w") do io
            write(io, tex_doc)
        end

        try
            run(Cmd(`$compiler -interaction=nonstopmode snippet.tex`, dir=temp_dir))
            run(`pdftocairo -png -r $dpi -singlefile $pdf_path $png_base`)
            return Cairo.read_from_png(png_path)
        finally
            rm(temp_dir, recursive=true, force=true)
        end
    end

    raw_w = Float64(Cairo.width(img))
    raw_h = Float64(Cairo.height(img))

    pt_scale = (72.0 / dpi) * (fontsize / 12.0) * scale
    draw_w = raw_w * pt_scale
    draw_h = raw_h * pt_scale

    x_offset = if halign == :center
        -draw_w / 2.0
    elseif halign == :right
        -draw_w
    else # :left
        0.0
    end

    y_offset = if valign == :center
        -draw_h / 2.0
    elseif valign == :top
        0.0
    elseif valign == :bottom
        -draw_h
    else # :baseline
        -draw_h * 0.80
    end

    Cairo.save(cr)
    Cairo.translate(cr, x + x_offset, y + y_offset)
    Cairo.scale(cr, pt_scale, pt_scale)
    Cairo.set_source_surface(cr, img, 0, 0)
    Cairo.paint(cr)
    Cairo.restore(cr)
    return
end

"""
    BasicCompGeometry.latex_to_pdf(pages, output_path;
                                  paperwidth=800,
                                  paperheight=800,
                                  margin=50,
                                  packages=["amsmath", "amssymb", "amsfonts", "microtype"],
                                  compiler="lualatex",
                                  passes=2,
                                  preamble="")

Compile one or more LaTeX pages into a multi-page PDF document.
Page dimensions are specified in big points (bp) to match Cairo points.

# Arguments
- `pages`: Either a single LaTeX string or a `Vector{<:AbstractString}` where each entry represents one page.
- `output_path::String`: Destination file path for the output PDF.
- `paperwidth::Real`: Canvas/page width in bp (default: `800`).
- `paperheight::Real`: Canvas/page height in bp (default: `800`).
- `margin::Real`: Page margins in bp (default: `50`).
- `packages::Vector{String}`: Standard packages to include in preamble.
- `compiler::String`: LaTeX engine to invoke (`"lualatex"` or `"pdflatex"`, default: `"lualatex"`).
- `passes::Int`: Number of compilation passes (default: `2`).
- `preamble::String`: Additional LaTeX preamble definitions (e.g. macros).
"""
function BasicCompGeometry.latex_to_pdf(
    pages::Union{AbstractVector{<:AbstractString}, AbstractString},
    output_path::String;
    paperwidth::Real = 800,
    paperheight::Real = 800,
    margin::Real = 50,
    packages::Vector{String} = ["amsmath", "amssymb", "amsfonts", "microtype"],
    compiler::String = "lualatex",
    passes::Int = 2,
    preamble::String = ""
)
    temp_dir = mktempdir()
    tex_path = joinpath(temp_dir, "document.tex")
    pdf_path = joinpath(temp_dir, "document.pdf")

    body = if pages isa AbstractVector
        join(pages, "\n\\clearpage\n")
    else
        String(pages)
    end

    pkg_lines = join(["\\usepackage{$p}" for p in packages], "\n")

    tex_content = """
\\documentclass[12pt]{article}
\\usepackage[paperwidth=$(paperwidth)bp,paperheight=$(paperheight)bp,margin=$(margin)bp]{geometry}
$pkg_lines
$preamble
\\pagestyle{empty}
\\begin{document}
$body
\\end{document}
"""

    open(tex_path, "w") do io
        write(io, tex_content)
    end

    for _ in 1:max(1, passes)
        run(Cmd(`$compiler -interaction=nonstopmode -output-directory=$temp_dir $tex_path`, dir=temp_dir))
    end

    if !isfile(pdf_path)
        error("LaTeX compilation failed to produce $pdf_path using $compiler")
    end

    mkpath(dirname(output_path))
    cp(pdf_path, output_path, force=true)
    rm(temp_dir, recursive=true, force=true)
    return output_path
end

"""
    BasicCompGeometry.pdf_merge(output_pdf::String, sources...)

Merge multiple PDF files or specific pages into `output_pdf` using `qpdf`.

Each element of `sources` can be:
- A filename `String` (all pages of that file).
- A `Tuple` or `Pair`: `(filename, page_spec)` where `page_spec` can be an `Int`, `String` (e.g. `"1-3"`), or `Vector{Int}`.
- A raw command-line argument string for `qpdf`.

# Example
```julia
pdf_merge("final.pdf", (diag_pdf, 1), (text_pdf, 1), (diag_pdf, 2), (text_pdf, "2-3"))
```
"""
function BasicCompGeometry.pdf_merge(output_pdf::String, sources...)
    qpdf_args = String[]
    for src in sources
        if src isa Pair || src isa Tuple
            file, spec = src
            push!(qpdf_args, String(file))
            if spec isa AbstractVector{<:Integer}
                push!(qpdf_args, join(spec, ","))
            else
                push!(qpdf_args, string(spec))
            end
        elseif src isa AbstractString
            # Check if source contains embedded space-separated args or just a filename
            tokens = split(String(src))
            append!(qpdf_args, tokens)
        else
            push!(qpdf_args, string(src))
        end
    end

    mkpath(dirname(output_pdf))
    run(`qpdf --empty --pages $qpdf_args -- $output_pdf`)
    return output_pdf
end

end # module
