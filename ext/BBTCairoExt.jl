module BBTCairoExt

using BasicCompGeometry
import BasicCompGeometry: cairo_draw_setup, cairo_set_line_width, description, open_canvas
import BasicCompGeometry: width, height
using BasicCompGeometry.BBT
import Cairo
import Colors

const _BBT_FILL_COLOR = Colors.RGBA{Float64}(1.0, 1.0, 0.0, 0.1)
const _BBT_STROKE_COLOR = Colors.RGB{Float64}(0.0, 0.0, 1.0)
const _BBT_BG_COLOR = Colors.RGB{Float64}(173 / 255, 216 / 255, 230 / 255)

function _draw_bbt_bg(canvas, bb, color)
    Cairo.set_source(canvas.cr, color)
    bl = bottom_left(bb)
    Cairo.rectangle(canvas, bl.x, bl.y, width(bb, 1), height(bb, 2))
    Cairo.fill(canvas)
end

function _draw_bbt_node(canvas, node, fill_color, stroke_color)
    Cairo.set_source(canvas.cr, fill_color)
    bb = node.bb + (diam(node.bb) * 0.01)
    bl = bottom_left(bb)
    Cairo.rectangle(canvas, bl.x, bl.y, width(bb, 1), height(bb, 2))
    Cairo.fill_preserve(canvas)
    Cairo.set_source(canvas.cr, stroke_color)
    cairo_set_line_width(canvas, 1.5)
    Cairo.stroke(canvas)
end

function _collect_bbt_levels!(levels::Vector{Vector{N}}, node::N, level::Int) where {N}
    node.f_leaf && return
    push!(levels[level+1], node)
    !isnothing(node.left) && _collect_bbt_levels!(levels, node.left, level + 1)
    !isnothing(node.right) && _collect_bbt_levels!(levels, node.right, level + 1)
end

function _draw_bbt_levels!(canvas, levels, level_range)
    for lvl in level_range
        for node in levels[lvl+1]
            _draw_bbt_node(canvas, node, _BBT_FILL_COLOR, _BBT_STROKE_COLOR)
        end
    end
end

"""
    BBT.Tree_draw(tree::Tree{2,T,S,V}, filename::String) where {T,S,V}

Implementation of `Tree_draw` using Cairo and Colors.
Only available when `Cairo` and `Colors` are loaded.
"""
function BBT.Tree_draw(tree::Tree{2,T,S,V}, filename::String) where {T,S,V}
    isnothing(tree.root) && throw(ArgumentError("cannot draw an empty BBT tree"))

    # Expand root BB for margins
    bb_root = tree.root.bb + (diam(tree.root.bb) * 0.2)
    cw, ch = 800, 800

    d = BBT.depth(tree.root)
    NodeT = BBT.Node{2,T,S,V}
    levels = [NodeT[] for _ = 0:d]
    _collect_bbt_levels!(levels, tree.root, 0)

    open_canvas(filename, cw, ch; title = "Bounding Box Tree Visualization") do canvas
        cairo_draw_setup(canvas, bb_root, cw, ch, 40)

        # Page 1: Overview
        _draw_bbt_bg(canvas, bb_root, _BBT_BG_COLOR)
        _draw_bbt_levels!(canvas, levels, 0:d)
        description(canvas, "Bounding Box Tree Overview (all levels)")
        Cairo.show_page(canvas)

        # Subsequent pages: Level-by-level visualization
        for i = (d-1):-1:0
            _draw_bbt_bg(canvas, bb_root, _BBT_BG_COLOR)
            _draw_bbt_levels!(canvas, levels, i:(i+1))
            description(canvas, "Bounding Box Tree Levels $i to $(i+1)")
            Cairo.show_page(canvas)
        end
    end
    println("Created tree visualization: $filename")
end

end # module
