# -----------------------------------------------------------------------------
# XML Generation & Compilation Pipelines
# -----------------------------------------------------------------------------

const DEFAULT_FALLBACK_IPESTYLE = raw"""<ipestyle name="basic">
<symbol name="mark/disk(sx)" transformations="translations">
<path fill="sym-stroke">
0.6 0 0 0.6 0 0 e
</path>
</symbol>
<symbol name="mark/circle(sx)" transformations="translations">
<path fill="sym-fill" stroke="sym-stroke">
0.6 0 0 0.6 0 0 e
</path>
</symbol>
<symbol name="mark/box(sx)" transformations="translations">
<path fill="sym-fill" stroke="sym-stroke">
-0.6 -0.6 m
0.6 -0.6 l
0.6 0.6 l
-0.6 0.6 l
h
</path>
</symbol>
<symbol name="arrow/pointed(spx)">
<path stroke="sym-stroke" fill="sym-stroke" pen="sym-pen">
0 0 m
-1 0.333 l
-0.8 0 l
-1 -0.333 l
h
</path>
</symbol>
<symbol name="arrow/normal(spx)">
<path stroke="sym-stroke" fill="sym-stroke" pen="sym-pen">
0 0 m
-1 0.333 l
-1 -0.333 l
h
</path>
</symbol>
<color name="red" value="1 0 0"/>
<color name="green" value="0 1 0"/>
<color name="blue" value="0 0 1"/>
<color name="yellow" value="1 1 0"/>
<color name="darkred" value="0.7 0 0"/>
<color name="darkgreen" value="0 0.5 0"/>
<color name="darkblue" value="0 0 0.6"/>
<color name="lightgreen" value="0.8 1 0.8"/>
<color name="lightblue" value="0.85 0.9 1"/>
<color name="lightred" value="1 0.85 0.85"/>
<color name="lightgray" value="0.9 0.9 0.9"/>
<color name="gray1" value="0.1 0.1 0.1"/>
<color name="gray2" value="0.2 0.2 0.2"/>
<color name="gray3" value="0.3 0.3 0.3"/>
<color name="gray4" value="0.4 0.4 0.4"/>
<color name="gray5" value="0.5 0.5 0.5"/>
<color name="gray6" value="0.6 0.6 0.6"/>
<color name="gray7" value="0.88 0.88 0.88"/>
<color name="gray8" value="0.8 0.8 0.8"/>
<color name="gray9" value="0.9 0.9 0.9"/>
<pen name="extra thin" value="0.05"/>
<pen name="thin" value="0.2"/>
<pen name="heavier" value="0.8"/>
<pen name="fat" value="1.2"/>
<pen name="ultrafat" value="2"/>
<pen name="ultrafat 4.0" value="4"/>
<pen name="ultrathin 0.5" value="0.5"/>
<pen name="ultrathin 0.25" value="0.25"/>
<opacity name="10%" value="0.1"/>
<opacity name="20%" value="0.2"/>
<opacity name="30%" value="0.3"/>
<opacity name="40%" value="0.4"/>
<opacity name="50%" value="0.5"/>
<opacity name="60%" value="0.6"/>
<opacity name="70%" value="0.7"/>
<opacity name="80%" value="0.8"/>
<opacity name="90%" value="0.9"/>
<textsize name="Huge" value="\Huge"/>
<textsize name="LARGE" value="\LARGE"/>
<textsize name="Large" value="\Large"/>
<textsize name="large" value="\large"/>
<textsize name="small" value="\small"/>
<textsize name="tiny" value="\tiny"/>
<textsize name="footnote" value="\footnotesize"/>
<tiling name="falling" angle="-60" step="4" width="1"/>
<tiling name="rising" angle="30" step="4" width="1"/>
<tiling name="hatch" angle="45" step="4" width="0.5"/>
<tiling name="crosshatch" angle="45" step="4" width="0.5"/>
<tiling name="vertical" angle="90" step="4" width="0.5"/>
<tiling name="horizontal" angle="0" step="4" width="0.5"/>
<layout paper="576 504" origin="0 0" frame="576 504" crop="yes"/>
</ipestyle>"""

"""
    to_xml(canvas)

Generate the complete, self-contained Ipe 7 XML document string.
"""
function to_xml(canvas::IpeCanvas)
    template =
        if canvas.stylesheet_template !== nothing && isfile(canvas.stylesheet_template)
            read(canvas.stylesheet_template, String)
        elseif isfile(DEFAULT_STYLE_FILE)
            read(DEFAULT_STYLE_FILE, String)
        else
            "<?xml version=\"1.0\"?>\n<!DOCTYPE ipe SYSTEM \"ipe.dtd\">\n<ipe version=\"70218\">\n<info tex=\"pdflatex\"/>\n" *
            DEFAULT_FALLBACK_IPESTYLE *
            "\n<page/>\n</ipe>"
        end

    # 1. Ensure bbox="cropbox" in info tag
    if occursin(r"<info\b", template)
        if !occursin("bbox=", template)
            template =
                replace(template, r"<info\b([^>]*)/>" => s"<info\1 bbox=\"cropbox\"/>")
        end
    end

    # 2. Ensure layout has crop="yes" and canvas paper sizing
    if occursin(r"<layout\b[^>]*/>", template)
        template = replace(
            template,
            r"<layout\b[^>]*/>" => "<layout paper=\"$(canvas.paper)\" origin=\"0 0\" frame=\"$(canvas.paper)\" crop=\"yes\"/>",
        )
    elseif occursin(r"</ipestyle>\s*<page>", template)
        style_layout = "<ipestyle name=\"page_layout\">\n<layout paper=\"$(canvas.paper)\" origin=\"0 0\" frame=\"$(canvas.paper)\" crop=\"yes\"/>\n</ipestyle>"
        template = replace(
            template,
            r"</ipestyle>\s*<page>" => "</ipestyle>\n" * style_layout * "\n<page>",
        )
    end

    # 3. Enhanced preamble
    preamble_xml = "<preamble>\n$(canvas.preamble)\n</preamble>"
    if occursin(r"<preamble>.*?</preamble>"s, template)
        template = replace(template, r"<preamble>.*?</preamble>"s => preamble_xml)
    end

    # 4. Build Page Content
    buf = IOBuffer()
    println(buf, "<page>")
    for layer in canvas.layers
        println(buf, "<layer name=\"$layer\"/>")
    end
    if isempty(canvas.views)
        println(
            buf,
            "<view layers=\"$(join(canvas.layers, " "))\" active=\"$(canvas.active_layer)\"/>",
        )
    else
        for v in canvas.views
            println(buf, v)
        end
    end
    for elem in canvas.elements
        println(buf, elem)
    end
    println(buf, "</page>")
    page_xml = String(take!(buf))

    # Replace <page>...</page> block using a function to avoid backslash escape issues
    result = replace(template, r"<page>.*?</page>"s => (m -> page_xml))
    return result
end

"""
    save_ipe(canvas, filename)

Save the canvas content to an `.ipe` XML file.
"""
function save_ipe(canvas::IpeCanvas, filename::String)
    mkpath(dirname(abspath(filename)))
    xml_data = to_xml(canvas)
    write(filename, xml_data)
    return filename
end

"""
    compile_pdf(canvas_or_file, output_pdf=nothing)

Compile the canvas or an `.ipe` file to a cropped vector PDF using `ipetoipe -pdf`.
"""
function compile_pdf(canvas::IpeCanvas, output_pdf::String)
    if Sys.which("ipetoipe") === nothing
        @warn "ipetoipe command not found in PATH. Please install Ipe 7 to compile PDF figures."
        return false
    end
    mkpath(dirname(abspath(output_pdf)))
    ipe_temp = tempname() * ".ipe"
    save_ipe(canvas, ipe_temp)
    success = run(`ipetoipe -pdf $ipe_temp $output_pdf`).exitcode == 0
    rm(ipe_temp, force = true)
    return success
end

function compile_pdf(
    ipe_file::String,
    output_pdf::String = replace(ipe_file, r"\.ipe$" => ".pdf"),
)
    if Sys.which("ipetoipe") === nothing
        @warn "ipetoipe command not found in PATH. Please install Ipe 7 to compile PDF figures."
        return false
    end
    mkpath(dirname(abspath(output_pdf)))
    run(`ipetoipe -pdf $ipe_file $output_pdf`).exitcode == 0
end

"""
    save_figure_tex(filename; figure_name, caption="", label="")

Emit a standalone companion LaTeX figure fragment (`_fig.tex`) ready for `\\input`.
"""
function save_figure_tex(
    filename::String;
    figure_name::String = replace(basename(filename), r"(_fig)?\.tex$" => ""),
    caption::String = "",
    label::String = figure_name,
)
    mkpath(dirname(abspath(filename)))
    tex = """
\\begin{figure}[t]
    \\centering
    \\IncludeGraphics{\\File{figs/$(figure_name)}}
    \\caption{$(caption)}
    \\figlab{$(label)}
\\end{figure}
"""
    write(filename, tex)
    return filename
end

"""
    export_figure(canvas, base_path; caption="", label="")

Export all three companion artifacts in one call:
`<base_path>.ipe`, `<base_path>.pdf`, and `<base_path>_fig.tex`.
"""
function export_figure(
    canvas::IpeCanvas,
    base_path::String;
    caption::String = "",
    label::String = replace(basename(base_path), "_" => ":"),
    outputs = (:ipe, :pdf, :tex),
)
    requested = Set(Symbol.(outputs isa Symbol ? (outputs,) : outputs))
    all(output in (:ipe, :pdf, :tex) for output in requested) ||
        throw(ArgumentError("outputs may contain only :ipe, :pdf, and :tex"))
    ipe_path = base_path * ".ipe"
    pdf_path = base_path * ".pdf"
    tex_path = base_path * "_fig.tex"

    if :ipe in requested || :pdf in requested
        save_ipe(canvas, ipe_path)
    end
    compiled = :pdf in requested ? compile_pdf(ipe_path, pdf_path) : true
    :tex in requested && save_figure_tex(
        tex_path;
        figure_name = basename(base_path),
        caption = caption,
        label = label,
    )
    if !compiled
        @warn "PDF compilation failed; preserving Ipe source" path=ipe_path
        throw(ErrorException("PDF output was requested, but Ipe compilation failed"))
    elseif !(:ipe in requested)
        rm(ipe_path; force = true)
    end

    return (
        ipe = :ipe in requested ? ipe_path : nothing,
        pdf = :pdf in requested ? pdf_path : nothing,
        tex = :tex in requested ? tex_path : nothing,
    )
end

"""
    open_ipe(base_path; caption="", label="", preview=false, kwargs...) do canvas
        ...
    end

Convenient block syntax to construct, save, compile, and generate LaTeX fragments
in one go. With `preview=true`, open the retained `.ipe` source in Ipe.
"""
function open_ipe(
    f::Function,
    base_path::String;
    caption::String = "",
    label::String = "",
    outputs = (:ipe, :pdf, :tex),
    fit = nothing,
    margin::Real = 30.0,
    flip_y::Bool = false,
    preview::Bool = false,
    kwargs...,
)
    requested = Set(Symbol.(outputs isa Symbol ? (outputs,) : outputs))
    preview &&
        !(:ipe in requested) &&
        throw(ArgumentError("preview=true requires :ipe in outputs"))
    clean_base = replace(base_path, r"(\.ipe|\.pdf|_fig\.tex)$" => "")
    canvas = IpeCanvas(; kwargs...)
    fit === nothing || fit!(canvas, fit; margin = margin, flip_y = flip_y)
    f(canvas)
    lbl = isempty(label) ? replace(basename(clean_base), "_" => ":") : label
    artifacts =
        export_figure(canvas, clean_base; caption = caption, label = lbl, outputs = outputs)
    if preview
        edit_ipe(artifacts.ipe)
    end
    return artifacts
end

"""
    figure(f, path; keep_source=true, source=:ipe, preview=false, kwargs...)

Generate a figure from geometry. A `.pdf` target is rendered through Ipe so that
geometry and LaTeX labels remain vector content. The editable `.ipe` source is
retained by default.
Set `preview=true` to retain and open the editable source after generation.
"""
function figure(
    f::Function,
    path::String;
    keep_source::Bool = true,
    source::Symbol = :ipe,
    preview::Bool = false,
    kwargs...,
)
    source == :ipe || throw(ArgumentError("the only supported figure source is :ipe"))
    base, extension = splitext(path)
    if isempty(extension)
        base = path
        outputs = keep_source || preview ? (:ipe, :pdf) : (:pdf,)
    elseif extension == ".pdf"
        outputs = keep_source || preview ? (:ipe, :pdf) : (:pdf,)
    elseif extension == ".ipe"
        outputs = (:ipe,)
    else
        throw(ArgumentError("figure target must have extension .pdf or .ipe"))
    end
    return open_ipe(f, base; outputs = outputs, preview = preview, kwargs...)
end

"""
    edit_ipe(filename::String)

Launch the Ipe GUI editor on `filename` in the background.
"""
function edit_ipe(filename::String)
    if Sys.which("ipe") === nothing
        @warn "ipe executable not found in PATH. Please install Ipe 7."
        return nothing
    end
    run(`ipe $filename`, wait = false)
end
