# Vector figures

`IpeDraw` writes editable Ipe 7 figures with geometry and LaTeX labels. Its
high-level API fits world coordinates to a page while keeping label offsets,
marks, and annotations in page units.

## A first figure

Target a `.ipe` path to write editable source without external compilation:

```julia
using BasicCompGeometry
using BasicCompGeometry.IpeDraw

c1 = Circle(point(0.0, 0.0), 1.0)
c2 = Circle(point(1.0, 0.0), 1.0)
p, q = sort(intersections(c1, c2); by=p -> -p.y)

figure("output/disks.ipe"; fit=[c1, c2], margin=12) do fig
    disk = Style(fill_opacity=0.2, pen=:heavier)
    draw!(fig, c1; style=disk, stroke=:blue, fill=:lightblue)
    draw!(fig, c2; style=disk, stroke=:darkgreen, fill=:lightgreen)
    mark!(fig, [p, q])
    label!(fig, p, raw"p"; offset=(6, 4), anchor=:southwest)
    label!(fig, q, raw"\sqrt{q}"; offset=(6, -6), anchor=:northwest)
end
```

Change the target to `output/disks.pdf` to compile a cropped vector PDF and keep
`output/disks.ipe`. PDF targets require Ipe's `ipetoipe` command on `PATH` and
the LaTeX tools needed by the labels. Pass `keep_source=false` when only the PDF
is wanted, or `preview=true` to open the source in the Ipe editor.
Plain LaTeX strings work without `LaTeXStrings`; load that package if using its
`L"..."` literals.

## Styles and page composition

`Style` collects reusable drawing options. `Theme` names a collection of styles;
`publication_theme()` supplies restrained defaults. In an open `figure` block,
insets and clipping preserve the surrounding viewport:

```julia
theme = Theme(region=Style(fill=:lightblue, stroke=:blue, fill_opacity=0.2))
detail_bbox = BBox(point(-0.2, -0.2), point(1.2, 1.2))

figure("output/inset.ipe"; fit=[c1, c2]) do fig
    draw!(fig, [c1, c2]; style=theme.region)
    inset(fig, PageBox(400, 320, 150, 140); fit=detail_bbox) do zoom
        clip_to(zoom, detail_bbox) do clipped
            draw!(clipped, [c1, c2]; style=theme.region)
        end
    end
    legend!(fig, ["feasible" => theme.region]; position=:northwest)
    scale_bar!(fig, 0.5; label=raw"1/2")
end
```

Drawing methods cover points, segments, boxes, circles and arcs, ellipses,
Bézier curves, and splines. Other helpers draw polygons with holes, B-splines,
scoped affine groups, intervals, dimensions, arrows, hatch patterns, layers,
and progressive presentation views.

See the [canvas and composition reference](reference/ipe.md) and the
[drawing and export reference](reference/ipe-drawing.md).

## Low-level page coordinates

`open_ipe` exposes the page-coordinate drawing interface:

```julia
using BasicCompGeometry
using BasicCompGeometry.IpeDraw
using LaTeXStrings

open_ipe("output/interval_demo"; caption="Interval estimation", label="fig:demo") do cv
    add_layer!(cv, "intervals")
    add_layer!(cv, "labels")
    set_layer!(cv, "intervals")
    draw_bar!(cv, 50.0, 350.0, 100.0; label_left=L"0", label_right=L"1")
    draw_span!(cv, 100.0, 250.0, 100.0; fill=:lightgreen, stroke=:darkgreen)
    draw_dimension!(cv, 100.0, 250.0, 130.0; label=L"\Delta \le \varepsilon", arrow=:both)
    set_layer!(cv, "labels")
    draw_point!(cv, 175.0, 100.0; stroke=:darkblue, fill=:darkblue, shape=:disk)
    draw_label!(cv, 175.0, 85.0, L"\mu"; halign=:center)
end
```

## Layers, views, and per-view cropping

An Ipe page can hold several views, each showing a subset of the page's layers;
`ipetoipe -pdf` writes every view as a separate PDF page, which suits
step-by-step slides and figure sequences. By default Ipe crops all views of a
page to the bounding box of the whole page, so an early view showing a small
part of the figure carries the blank space of the later ones.

Ipe crops a view that contains the special `VIEWBBOX` layer to its own visible
objects instead. The layer may stay empty. Add it to each view that should be
cropped independently, either explicitly:

```julia
add_layer!(fig, "VIEWBBOX")
add_view!(fig, ["alpha", "VIEWBBOX"]; active="alpha")
```

or with the `autocrop` flag, which registers the layer and appends it to the
view:

```julia
figure("output/views.pdf"; fit=[outer]) do fig
    layer(fig, "triangle") do f
        draw!(f, [a, b, c])
    end
    layer(fig, "circumcircle") do f
        draw!(f, outer)
    end
    add_view!(fig, ["triangle"]; autocrop=true)                  # tight around the triangle
    add_view!(fig, ["triangle", "circumcircle"]; autocrop=true)  # tight around the circle
end
```

When a page has no explicit views, IpeDraw emits a single view of all layers.
The [view autocrop example](https://github.com/sarielhp/BasicCompGeometry.jl/blob/main/examples/ipe_view_autocrop.jl)
writes the same three-view page with and without cropping and prints the PDF
page sizes of both.

## Complete examples

The [convex hull figure](https://github.com/sarielhp/BasicCompGeometry.jl/blob/main/examples/convex_hull_figure.jl)
draws the hull and diameter of 35 deterministic points, with math labels, a
legend, a scale bar, and an inset. Run `./examples/convex_hull_figure.jl` from
the repository to produce `output/convex_hull_figure.ipe` and
`output/convex_hull_figure.pdf`. The script was generated with Gemini Flash and
uses [QuickEnv.jl](https://github.com/sarielhp/QuickEnv.jl) for dependencies.

See also the
[conceptual figure example](https://github.com/sarielhp/BasicCompGeometry.jl/blob/main/examples/ipe_conceptual_figure.jl)
and the [full example catalog](https://github.com/sarielhp/BasicCompGeometry.jl/blob/main/examples/README.md).
