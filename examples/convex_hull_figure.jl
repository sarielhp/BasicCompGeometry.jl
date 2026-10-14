#!/usr/bin/env julia

# Using QuickEnv for zero-ceremony standalone execution without managing environments.
include(joinpath(@__DIR__, "support/runtime.jl"))

using QuickEnv
using BasicCompGeometry
using BasicCompGeometry.IpeDraw
using Random

# 1. Deterministic point generation
Random.seed!(42)
n = 35
pts = [point(rand() * 10.0, rand() * 10.0) for _ in 1:n]

# 2. Geometric computations
ch = convex_hull(pts)

# Compute the diameter pair among hull vertices
max_d = 0.0
p_diam, q_diam = pts[1], pts[2]
for i in 1:(length(ch) - 1), j in (i + 1):length(ch)
    d = dist(ch[i], ch[j])
    if d > max_d
        global max_d = d
        global p_diam, q_diam = ch[i], ch[j]
    end
end
diam_seg = Segment(p_diam, q_diam)
diam_mid = point((p_diam.x + q_diam.x) / 2, (p_diam.y + q_diam.y) / 2)

# 3. Output directory setup
output_dir = example_output()
mkpath(output_dir)
pdf_path = joinpath(output_dir, "convex_hull_figure.pdf")

# 4. Styling theme
theme = Theme(
    hull = Style(fill=:lightblue, stroke=:navy, fill_opacity=0.2, pen=:heavier),
    diam = Style(stroke=:red, pen=:heavy, dashstyle=:dashed),
)

# 5. Render publication figure
figure(pdf_path; fit=pts, margin=20, preview="--preview" in ARGS) do fig
    # Draw convex hull and interior points
    draw!(fig, ch; style=theme.hull)
    mark!(fig, pts; stroke=:gray, fill=:gray, shape=:disk, size=:tiny)
    mark!(fig, ch; stroke=:navy, fill=:navy, shape=:disk, size=:small)

    # Draw diameter segment and endpoint marks
    draw!(fig, diam_seg; style=theme.diam)
    mark!(fig, (p_diam, q_diam); stroke=:red, fill=:red, shape=:disk, size=:normal)

    # LaTeX mathematical labels
    label!(fig, p_diam, raw"p"; offset=(-8, -8), anchor=:northeast)
    label!(fig, q_diam, raw"q"; offset=(8, 8), anchor=:southwest)
    label!(fig, diam_mid, raw"D = \mathrm{diam}(P)"; offset=(8, -8), anchor=:northwest)

    # Legend and scale bar in page coordinates
    legend!(fig, [raw"Convex Hull $\mathcal{CH}(P)$" => theme.hull, raw"Diameter $D$" => theme.diam];
            position=:northwest, width=155)
    scale_bar!(fig, 2.0; label=raw"2.0")

    # Corner inset zooming into diameter endpoint p
    zoom_box = BBox(point(p_diam.x - 1.2, p_diam.y - 1.2),
                    point(p_diam.x + 1.2, p_diam.y + 1.2))
    inset(fig, PageBox(400, 320, 150, 150); fit=zoom_box) do zoom
        draw!(zoom, ch; style=theme.hull)
        mark!(zoom, pts; stroke=:gray, fill=:gray, shape=:disk, size=:tiny)
        draw!(zoom, diam_seg; style=theme.diam)
        mark!(zoom, p_diam; stroke=:red, fill=:red, shape=:disk, size=:normal)
        label!(zoom, p_diam, raw"p"; offset=(-6, 6), anchor=:southeast)
    end
end

println("Generated: ", pdf_path)
