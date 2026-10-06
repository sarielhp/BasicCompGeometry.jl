#!/usr/bin/env julia

using BasicCompGeometry
using BasicCompGeometry.IpeDraw

output_dir = joinpath(@__DIR__, "..", "output")
mkpath(output_dir)

left = Circle(point(0.0, 0.0), 1.0)
right = Circle(point(1.0, 0.0), 1.0)
upper, lower = sort(intersections(left, right); by=p -> -p.y)

figure(joinpath(output_dir, "disks_figure.pdf"); fit=(left, right), margin=18) do fig
    disk = Style(fill_opacity=0.2, pen=:heavier)
    draw!(fig, left; style=disk, stroke=:blue, fill=:lightblue)
    draw!(fig, right; style=disk, stroke=:darkgreen, fill=:lightgreen)
    mark!(fig, (upper, lower))
    label!(fig, upper, raw"p"; offset=(6, 5), anchor=:southwest)
    label!(fig, lower, raw"q"; offset=(6, -5), anchor=:northwest)
end
