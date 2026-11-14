#!/usr/bin/env julia

include(joinpath(@__DIR__, "support/runtime.jl"))

using BasicCompGeometry
using BasicCompGeometry.IpeDraw

#=
Independent cropping of several views on one Ipe page.

The page holds a triangle, its incircle, and its circumcircle on separate
layers, and three views that reveal them one after the other. `ipetoipe -pdf`
turns each view into one PDF page. Ipe normally crops every view to the bounding
box of the whole page, so the first view (the triangle alone) would carry the
blank space of the circumcircle. A view that contains Ipe's special `VIEWBBOX`
layer is instead cropped to its own visible objects; `add_view!(...;
autocrop=true)` adds that layer. The script writes both versions for comparison
and prints the page sizes of each PDF.
=#

output_dir = example_output()
mkpath(output_dir)

a, b, c = point(0.0, 0.0), point(4.0, 0.0), point(1.0, 2.0)

function circumcircle(a, b, c)
    d = 2 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y))
    sq(p) = p.x^2 + p.y^2
    ux = (sq(a) * (b.y - c.y) + sq(b) * (c.y - a.y) + sq(c) * (a.y - b.y)) / d
    uy = (sq(a) * (c.x - b.x) + sq(b) * (a.x - c.x) + sq(c) * (b.x - a.x)) / d
    center = point(ux, uy)
    return Circle(center, dist(center, a))
end

function incircle(a, b, c)
    la, lb, lc = dist(b, c), dist(c, a), dist(a, b)
    perimeter = la + lb + lc
    center = point(
        (la * a.x + lb * b.x + lc * c.x) / perimeter,
        (la * a.y + lb * b.y + lc * c.y) / perimeter,
    )
    s = perimeter / 2
    area = sqrt(s * (s - la) * (s - lb) * (s - lc))
    return Circle(center, area / s)
end

outer = circumcircle(a, b, c)
inner = incircle(a, b, c)

function draw_views(path::String; autocrop::Bool)
    return figure(path; fit = [outer], margin = 12) do fig
        add_layer!(fig, "triangle")
        add_layer!(fig, "incircle")
        add_layer!(fig, "circumcircle")
        layer(fig, "triangle") do f
            draw!(f, [a, b, c]; stroke = :black, pen = :heavier)
            mark!(f, [a, b, c])
        end
        layer(fig, "incircle") do f
            draw!(f, inner; stroke = :blue, fill = :lightblue, fill_opacity = 0.3)
        end
        layer(fig, "circumcircle") do f
            draw!(f, outer; stroke = :darkred, pen = :heavier)
        end
        # Equivalent to: add_layer!(fig, "VIEWBBOX");
        #                add_view!(fig, ["triangle", "VIEWBBOX"]; active = "triangle")
        add_view!(fig, ["triangle"]; autocrop = autocrop)
        add_view!(fig, ["triangle", "incircle"]; active = "incircle", autocrop = autocrop)
        add_view!(
            fig,
            ["triangle", "incircle", "circumcircle"];
            active = "circumcircle",
            autocrop = autocrop,
        )
    end
end

function page_sizes(pdf::String)
    Sys.which("pdfinfo") === nothing && return "(pdfinfo not found)"
    info = read(`pdfinfo -f 1 -l 99 $pdf`, String)
    sizes = [m.captures[1] for m in eachmatch(r"Page\s+\d+ size:\s+(.*?) pts", info)]
    return join(sizes, "; ")
end

for (name, crop) in (("ipe_view_autocrop", true), ("ipe_view_nocrop", false))
    pdf = joinpath(output_dir, name * ".pdf")
    draw_views(pdf; autocrop = crop)
    isfile(pdf) && println(rpad(name, 20), page_sizes(pdf))
end
