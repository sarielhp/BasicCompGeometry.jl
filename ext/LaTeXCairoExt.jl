module LaTeXCairoExt

using BasicCompGeometry
import Cairo
import LaTeXStrings: LaTeXString

# --- Global Session LaTeX State ---
const DEFAULT_PACKAGES = ["amsmath", "amssymb", "amsfonts", "microtype"]
const GLOBAL_PACKAGES = Set{String}(DEFAULT_PACKAGES)
const GLOBAL_MACROS = Ref{String}("")

const LATEX_SNIPPET_CACHE = Dict{Tuple{String, Float64, Int, String, String}, Cairo.CairoSurfaceBase{UInt32}}()

"""
    BasicCompGeometry.read_latex_snippet(file_or_content; beg_marker=nothing, end_marker=nothing)

Read text from a `.tex` file or string, optionally extracting the section between `beg_marker` and `end_marker`.
If no markers are given, the entire content is returned.
"""
function BasicCompGeometry.read_latex_snippet(
    file_or_content::AbstractString;
    beg_marker::Union{AbstractString, Regex, Nothing} = nothing,
    end_marker::Union{AbstractString, Regex, Nothing} = nothing
)
    text = isfile(file_or_content) ? read(file_or_content, String) : String(file_or_content)
    isnothing(beg_marker) && return strip(text)

    r_beg = beg_marker isa Regex ? beg_marker : Regex(escape_string(String(beg_marker)))
    m_beg = match(r_beg, text)
    isnothing(m_beg) && return strip(text)

    sub = text[(m_beg.offset + length(m_beg.match)):end]
    isnothing(end_marker) && return strip(sub)

    r_end = end_marker isa Regex ? end_marker : Regex(escape_string(String(end_marker)))
    m_end = match(r_end, sub)
    return isnothing(m_end) ? strip(sub) : strip(sub[1:m_end.offset - 1])
end

"""
    BasicCompGeometry.append_latex_preamble!(source; beg_marker=nothing, end_marker=nothing)

Append LaTeX packages or macros from a string or `.tex` file to the global session preamble.
If `source` is a file path, automatically extracts delimited snippets (e.g. `%%% IPE Prelim start`).
"""
function BasicCompGeometry.append_latex_preamble!(
    source::Union{AbstractString, LaTeXString};
    beg_marker::Union{AbstractString, Regex, Nothing} = nothing,
    end_marker::Union{AbstractString, Regex, Nothing} = nothing
)
    snippet = BasicCompGeometry.read_latex_snippet(String(source); beg_marker=beg_marker, end_marker=end_marker)
    if !isempty(snippet)
        if isempty(GLOBAL_MACROS[])
            GLOBAL_MACROS[] = snippet
        else
            GLOBAL_MACROS[] = GLOBAL_MACROS[] * "\n" * snippet
        end
    end
    empty!(LATEX_SNIPPET_CACHE)
    return GLOBAL_MACROS[]
end

BasicCompGeometry.add_latex_macros!(source; kwargs...) = BasicCompGeometry.append_latex_preamble!(source; kwargs...)

"""
    BasicCompGeometry.set_latex_preamble!(source; beg_marker=nothing, end_marker=nothing)

Set the global session LaTeX preamble from a string or `.tex` file, replacing existing macros.
"""
function BasicCompGeometry.set_latex_preamble!(
    source::Union{AbstractString, LaTeXString};
    beg_marker::Union{AbstractString, Regex, Nothing} = nothing,
    end_marker::Union{AbstractString, Regex, Nothing} = nothing
)
    snippet = BasicCompGeometry.read_latex_snippet(String(source); beg_marker=beg_marker, end_marker=end_marker)
    GLOBAL_MACROS[] = snippet
    empty!(LATEX_SNIPPET_CACHE)
    return GLOBAL_MACROS[]
end

"""
    BasicCompGeometry.add_latex_packages!(pkgs...)

Add one or more LaTeX package names to the global session package list.
"""
function BasicCompGeometry.add_latex_packages!(pkgs...)
    for p in pkgs
        if p isa AbstractVector || p isa Tuple
            for sub in p
                push!(GLOBAL_PACKAGES, String(sub))
            end
        else
            push!(GLOBAL_PACKAGES, String(p))
        end
    end
    empty!(LATEX_SNIPPET_CACHE)
    return collect(GLOBAL_PACKAGES)
end

"""
    BasicCompGeometry.reset_latex_preamble!()

Reset session LaTeX packages and custom macros to library defaults.
"""
function BasicCompGeometry.reset_latex_preamble!()
    empty!(GLOBAL_PACKAGES)
    for p in DEFAULT_PACKAGES
        push!(GLOBAL_PACKAGES, p)
    end
    GLOBAL_MACROS[] = ""
    empty!(LATEX_SNIPPET_CACHE)
    return
end

"""
    BasicCompGeometry.get_latex_preamble(; extra_packages=[], extra_preamble="")

Retrieve the combined LaTeX preamble string for the active session.
"""
function BasicCompGeometry.get_latex_preamble(;
    extra_packages::Vector{String} = String[],
    extra_preamble::AbstractString = ""
)
    all_pkgs = unique(vcat(collect(GLOBAL_PACKAGES), extra_packages))
    pkg_lines = join(["\\usepackage{$p}" for p in all_pkgs], "\n")
    macros = GLOBAL_MACROS[]
    combined = string(pkg_lines, "\n", macros, "\n", extra_preamble)
    return strip(combined)
end

"""
    BasicCompGeometry.cairo_draw_latex(cr, x, y, latex_str;
                                       fontsize=12.0,
                                       halign=:left,
                                       valign=:baseline,
                                       dpi=300,
                                       scale=1.0,
                                       compiler="lualatex",
                                       extra_packages=String[],
                                       preamble="")

Render a LaTeX snippet or mathematical formula `latex_str` onto the Cairo context `cr` at position `(x, y)`.
"""
function _run_silent(cmd::Cmd)
    buf = IOBuffer()
    try
        run(pipeline(cmd, stdout=buf, stderr=buf))
    catch e
        output = String(take!(buf))
        !isempty(output) && println(stderr, output)
        rethrow(e)
    end
end

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
    extra_packages::Vector{String} = String[],
    preamble::AbstractString = ""
)
    actual_cr = cr isa Cairo.CairoContext ? cr : (hasproperty(cr, :cr) ? cr.cr : cr)
    str = String(latex_str)
    full_preamble = BasicCompGeometry.get_latex_preamble(extra_packages=extra_packages, extra_preamble=preamble)
    cache_key = (str, Float64(fontsize), dpi, compiler, full_preamble)

    img = get!(LATEX_SNIPPET_CACHE, cache_key) do
        temp_dir = mktempdir()
        tex_path = joinpath(temp_dir, "snippet.tex")
        pdf_path = joinpath(temp_dir, "snippet.pdf")
        png_base = joinpath(temp_dir, "snippet")
        png_path = joinpath(temp_dir, "snippet.png")

        content = strip(str)
        is_math = startswith(content, "\$") || startswith(content, "\\[") || startswith(content, "\\begin{equation")
        body = is_math ? content : "\\textbf{\\textnormal{$content}}"

        tex_doc = """
\\documentclass[border=0pt,varwidth,12pt]{standalone}
$full_preamble
\\begin{document}
$body
\\end{document}
"""
        open(tex_path, "w") do io
            write(io, tex_doc)
        end

        try
            _run_silent(Cmd(`$compiler -interaction=nonstopmode snippet.tex`, dir=temp_dir))
            _run_silent(`pdftocairo -png -r $dpi -singlefile $pdf_path $png_base`)
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

    Cairo.save(actual_cr)
    Cairo.translate(actual_cr, x + x_offset, y + y_offset)
    Cairo.scale(actual_cr, pt_scale, pt_scale)
    Cairo.set_source_surface(actual_cr, img, 0, 0)
    Cairo.paint(actual_cr)
    Cairo.restore(actual_cr)
    return
end

"""
    BasicCompGeometry.latex_to_pdf(pages, output_path;
                                  paperwidth=800,
                                  paperheight=800,
                                  margin=50,
                                  extra_packages=String[],
                                  compiler="lualatex",
                                  passes=2,
                                  preamble="")

Compile one or more LaTeX pages into a multi-page PDF document.
Page dimensions are specified in big points (bp) to match Cairo points.
"""
function BasicCompGeometry.latex_to_pdf(
    pages::Union{AbstractVector{<:AbstractString}, AbstractString},
    output_path::String;
    paperwidth::Real = 800,
    paperheight::Real = 800,
    margin::Real = 50,
    extra_packages::Vector{String} = String[],
    compiler::String = "lualatex",
    passes::Int = 2,
    preamble::AbstractString = ""
)
    temp_dir = mktempdir()
    tex_path = joinpath(temp_dir, "document.tex")
    pdf_path = joinpath(temp_dir, "document.pdf")

    body = if pages isa AbstractVector
        join(pages, "\n\\clearpage\n")
    else
        String(pages)
    end

    full_preamble = BasicCompGeometry.get_latex_preamble(extra_packages=extra_packages, extra_preamble=preamble)

    tex_content = """
\\documentclass[12pt]{article}
\\usepackage[paperwidth=$(paperwidth)bp,paperheight=$(paperheight)bp,margin=$(margin)bp]{geometry}
$full_preamble
\\pagestyle{empty}
\\begin{document}
$body
\\end{document}
"""

    open(tex_path, "w") do io
        write(io, tex_content)
    end

    for _ in 1:max(1, passes)
        _run_silent(Cmd(`$compiler -interaction=nonstopmode -output-directory=$temp_dir $tex_path`, dir=temp_dir))
    end

    if !isfile(pdf_path)
        log_file = joinpath(temp_dir, "document.log")
        log_msg = isfile(log_file) ? "\n" * read(log_file, String) : ""
        error("LaTeX compilation failed to produce $pdf_path using $compiler$log_msg")
    end

    mkpath(dirname(output_path))
    cp(pdf_path, output_path, force=true)
    rm(temp_dir, recursive=true, force=true)
    return output_path
end

"""
    BasicCompGeometry.pdf_merge(output_pdf::String, sources...)

Merge multiple PDF files or specific pages into `output_pdf` using `qpdf`.
"""
function BasicCompGeometry.pdf_merge(output_pdf::String, sources...)
    mkpath(dirname(output_pdf))
    has_qpdf = Sys.which("qpdf") !== nothing
    has_mutool = Sys.which("mutool") !== nothing

    if has_qpdf
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
                tokens = split(String(src))
                append!(qpdf_args, tokens)
            else
                push!(qpdf_args, string(src))
            end
        end
        _run_silent(`qpdf --empty --pages $qpdf_args -- $output_pdf`)
    elseif has_mutool
        mutool_args = String["merge", "-o", output_pdf]
        for src in sources
            if src isa Pair || src isa Tuple
                file, spec = src
                push!(mutool_args, String(file))
                if spec isa AbstractVector{<:Integer}
                    push!(mutool_args, join(spec, ","))
                else
                    push!(mutool_args, string(spec))
                end
            elseif src isa AbstractString
                tokens = split(String(src))
                append!(mutool_args, tokens)
            else
                push!(mutool_args, string(src))
            end
        end
        _run_silent(`mutool $mutool_args`)
    else
        error("pdf_merge requires 'qpdf' or 'mutool' to be installed and in PATH.")
    end
    return output_pdf
end

end # module

