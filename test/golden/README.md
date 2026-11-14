# Golden fixtures

Byte-for-byte reference outputs for the package's file writers and readers,
checked by [`test/test_golden.jl`](../test_golden.jl). They were recorded from
known-good output. Any change to an output format (coordinate rounding, line
layout, attribute order, binary layout, and so on) makes the test fail and
report the first line that differs.

| Fixture | Produced by |
|---|---|
| `ipe_primitives.ipe` | IpeDraw: every primitive (marks, segments, boxes, circles, arcs, ellipses, elliptic arcs, Bézier, Catmull-Rom/natural/explicit splines, open/closed B-splines, polygons, polygons with holes) at coordinates with rounding ties, signed zeros, tiny negatives, and large values |
| `ipe_viewport.ipe`, `ipe_viewport_flip.ipe` | IpeDraw: the same objects drawn through `fit!`, without and with `flip_y` |
| `ipe_labels_styles.ipe` | IpeDraw: labels (math wrapping, XML escaping, anchors, offsets), `Style`, `Theme`, `publication_theme`, `with_style`, `fill_opacity`, `draw!`/`mark!`/`label!`, layers, views, preamble |
| `ipe_diagrams.ipe` | IpeDraw: groups, `clip_to`, `page_space`, `inset`, `legend!`, `scale_bar!`, bars, spans, dimensions, arrows, curved arrows |
| `figure_fig.tex`, `export_fig.tex` | `save_figure_tex`, `export_figure(...; outputs = (:tex,))` |
| `matrix.fbin`, `matrix.ibin`, `empty.fbin` | `write_fbin`, `write_ibin` |
| `readers_dump.txt` | Exact dump of `read_plt`, `read_geolife_plt`, and `read_fvecs` applied to the inputs below |
| `write_plt.plt`, `write_plt_3d.plt` | `write_plt` |
| `halfplanes.txt` | `write_halfplanes` |

The files `input_points.plt`, `input_geolife.plt`, and `input.fvecs` are
hand-written reader inputs, not generated outputs. `input.fvecs` holds three
2-dimensional vectors: `[0.0025, -1.5]`, `[1f-30, Inf]`, `[-0.0, 3.1415927]`.

## Regenerating

Only regenerate when an output change is intended:

```bash
BCG_UPDATE_GOLDEN=1 julia --project=. -e 'using Test; include("test/test_golden.jl")'
```

In this mode the test writes the fixtures instead of comparing them and
records each one as broken (nothing was verified). Review the change with
`git diff test/golden/` before committing. Every changed line is an output
change that users will see. Do not regenerate merely to make the test pass.
