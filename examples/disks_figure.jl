#!/usr/bin/env julia

using BasicCompGeometry
using BasicCompGeometry.IpeDraw

output_dir = joinpath(@__DIR__, "..", "output")
mkpath(output_dir)

left = Circle(point(0.0, 0.0), 1.0)
right = Circle(point(1.0, 0.0), 1.0)
upper, lower = sort(intersections(left, right); by=p -> -p.y)
theme = Theme(
    left=Style(fill=:lightblue, stroke=:blue, fill_opacity=0.2, pen=:heavier),
    right=Style(fill=:lightgreen, stroke=:darkgreen, fill_opacity=0.2, pen=:heavier),
)

figure(joinpath(output_dir, "disks_figure.pdf"); fit=(left, right), margin=18,
       preview="--preview" in ARGS) do fig
    draw!(fig, left; style=theme.left)
    draw!(fig, right; style=theme.right)
    mark!(fig, (upper, lower))
    label!(fig, upper, raw"p"; offset=(6, 5), anchor=:southwest)
    label!(fig, lower, raw"q"; offset=(6, -5), anchor=:northwest)
    scale_bar!(fig, 0.5; label=raw"1/2")
    legend!(fig, ["left disk" => theme.left, "right disk" => theme.right];
            position=:northwest, width=115)

    detail = BBox(point(upper.x - 0.25, upper.y - 0.2),
                  point(upper.x + 0.25, upper.y + 0.2))
    inset(fig, PageBox(410, 330, 145, 145); fit=detail) do zoom
        draw!(zoom, left; style=theme.left)
        draw!(zoom, right; style=theme.right)
        mark!(zoom, upper)
    end
end
