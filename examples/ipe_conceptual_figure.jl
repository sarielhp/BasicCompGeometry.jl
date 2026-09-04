#!/usr/bin/env julia

using QuickEnv
using BasicCompGeometry
using BasicCompGeometry.IpeDraw
using LaTeXStrings

"""
    Example: Conceptual Ipe Figure Generation (Inverse Estimation via Sampling)

Demonstrates how to use `BasicCompGeometry.IpeDraw` to programmatically generate
publication-quality, vector conceptual figures with native LaTeX math annotations.

Outputs:
  - `output/inverse_estimation_demo.ipe`   (Fully editable Ipe 7 XML document)
  - `output/inverse_estimation_demo.pdf`   (Tight, cropped publication PDF)
  - `output/inverse_estimation_demo_fig.tex` (Companion LaTeX figure environment)

Usage:
  ./examples/ipe_conceptual_figure.jl
"""

function generate_inverse_estimation_figure(output_base::String)
    open_ipe(output_base;
        caption = "Inverse estimation via sampling (Lemma 1.3.1). A good sample with error \$\\le \\eps \\mu\$ in the active range yields an estimate \$\\in [x - \\eps, x + \\eps]\$.",
        label = "inverse:estimation:demo",
        width = 600.0,
        height = 360.0
    ) do cv
        # 1. Register custom LaTeX macros for this figure
        add_preamble!(cv, raw"""
        \newcommand{\eps}{\varepsilon}
        \newcommand{\E}{\mathcal{E}}
        """)

        # 2. Layers
        add_layer!(cv, "background")
        add_layer!(cv, "intervals")
        add_layer!(cv, "markers")
        add_layer!(cv, "labels")

        # -------------------------------------------------------------
        # Level 1 (Top): Transformed Space S_i
        # -------------------------------------------------------------
        set_layer!(cv, "background")
        y1 = 230.0
        x_min, x_max = 110.0, 470.0
        draw_bar!(cv, x_min, x_max, y1;
            stroke = :black,
            fill = :gray7,
            pen = :heavier,
            label_left = L"0",
            label_right = L"M = \frac{1}{\eps}"
        )

        set_layer!(cv, "intervals")
        # Danger / Error zone
        draw_span!(cv, 180.0, 400.0, y1; fill=:lightred, stroke=:darkred, pen=:heavier)
        # Success zone
        draw_span!(cv, 230.0, 350.0, y1; fill=:lightgreen, stroke=:darkgreen, pen=:heavier)

        # Dimension double-arrow
        draw_dimension!(cv, 230.0, 350.0, y1 + 30.0;
            label = L"\Delta \le 2\eps \mu",
            arrow = :both,
            stroke = :darkgreen,
            label_offset = 12.0
        )

        set_layer!(cv, "markers")
        # Center target point
        draw_point!(cv, 290.0, y1; stroke=:darkblue, fill=:darkblue, size=:large, shape=:disk)

        set_layer!(cv, "labels")
        draw_label!(cv, 80.0, y1 - 4.0, L"\mathcal{E}_{\le i}"; size=:LARGE, halign=:center)
        draw_label!(cv, 180.0, y1 - 22.0, L"s_{i-1}"; stroke=:darkred, size=:normal)
        draw_label!(cv, 400.0, y1 - 22.0, L"s_{i+1}"; stroke=:darkred, size=:normal)
        draw_label!(cv, 290.0, y1 + 14.0, L"\mu"; stroke=:darkblue, size=:normal)

        # -------------------------------------------------------------
        # Transition Arrows (Lifting / Conditioning)
        # -------------------------------------------------------------
        set_layer!(cv, "markers")
        draw_arrow!(cv, Point(230.0, 210.0), Point(200.0, 105.0); stroke=:darkgreen, pen=:fat, dash=:dashed)
        draw_arrow!(cv, Point(350.0, 210.0), Point(380.0, 105.0); stroke=:darkgreen, pen=:fat, dash=:dashed)

        # -------------------------------------------------------------
        # Level 2 (Bottom): Original Search Space
        # -------------------------------------------------------------
        set_layer!(cv, "background")
        y2 = 90.0
        draw_bar!(cv, 150.0, 430.0, y2;
            stroke = :black,
            fill = :gray7,
            pen = :heavier,
            label_left = L"0",
            label_right = L"1"
        )

        set_layer!(cv, "intervals")
        draw_span!(cv, 200.0, 380.0, y2; fill=:lightgreen, stroke=:darkgreen, pen=:heavier)

        # Bottom dimension line
        draw_dimension!(cv, 200.0, 380.0, y2 - 25.0;
            label = L"2\eps",
            arrow = :both,
            stroke = :darkgreen,
            label_offset = -12.0
        )

        set_layer!(cv, "markers")
        draw_point!(cv, 290.0, y2; stroke=:darkblue, fill=:darkblue, size=:large, shape=:disk)

        set_layer!(cv, "labels")
        draw_label!(cv, 105.0, y2 - 4.0, L"\mathrm{Output}"; size=:LARGE, halign=:center)
        draw_label!(cv, 200.0, y2 + 14.0, L"x - \eps"; stroke=:darkgreen, size=:normal)
        draw_label!(cv, 380.0, y2 + 14.0, L"x + \eps"; stroke=:darkgreen, size=:normal)
        draw_label!(cv, 290.0, y2 - 20.0, L"x"; stroke=:darkblue, size=:normal)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    out_dir = normpath(joinpath(@__DIR__, "..", "output"))
    mkpath(out_dir)
    target = joinpath(out_dir, "inverse_estimation_demo")
    artifacts = generate_inverse_estimation_figure(target)
    println("Generated artifacts:")
    println("  - Ipe XML : ", artifacts.ipe)
    println("  - PDF     : ", artifacts.pdf)
    println("  - LaTeX   : ", artifacts.tex)
end
