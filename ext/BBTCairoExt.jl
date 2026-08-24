module BBTCairoExt

using BasicCompGeometry
using BasicCompGeometry.BBT
import Cairo
import Colors

"""
    BBT.Tree_draw(tree::Tree{2,T,S,V}, filename::String) where {T,S,V}

Implementation of `Tree_draw` using Cairo and Colors.
Only available when `Cairo` and `Colors` are loaded.
"""
function BBT.Tree_draw(tree::Tree{2,T,S,V}, filename::String) where {T,S,V}
    C = Cairo
    Clr = Colors

    # Expand root BB for margins
    bb_root = tree.root.bb + (diam(tree.root.bb) * 0.2)
    cw, ch = 800, 800

    # Utility to draw a single bounding box
    function bbox_draw(canvas, bb, color)
        C.set_source(canvas, color)
        bl = bottom_left(bb)
        w = width(bb, 1)
        h = height(bb, 2)
        C.rectangle(canvas, bl[1], bl[2], w, h)
        C.fill(canvas)
    end

    # Recursive utility to draw nodes at specific levels
    function node_draw(canvas, node, level, range)
        node.f_leaf && return
        if level ∈ range
            yellow_transparent = Clr.coloralpha(Clr.parse(Clr.Colorant, "yellow"), 0.1)
            C.set_source(canvas, yellow_transparent)

            # Expand slightly for better visibility
            bb = node.bb + (diam(node.bb) * 0.01)
            bl = bottom_left(bb)
            w = width(bb, 1)
            h = height(bb, 2)

            C.rectangle(canvas, bl[1], bl[2], w, h)
            C.fill_preserve(canvas)

            C.set_source(canvas, Clr.parse(Clr.Colorant, "blue"))
            cairo_set_line_width(canvas, 1.5)
            C.stroke(canvas)
        end
        if !isnothing(node.left)
            node_draw(canvas, node.left, level + 1, range)
        end
        if !isnothing(node.right)
            node_draw(canvas, node.right, level + 1, range)
        end
    end

    d = BBT.depth(tree.root)
    open_canvas(filename, cw, ch; title="Bounding Box Tree Visualization") do canvas
        cairo_draw_setup(canvas, bb_root, cw, ch, 40)
        
        # Page 1: Overview
        bbox_draw(canvas, bb_root, Clr.parse(Clr.Colorant, "lightblue"))
        node_draw(canvas, tree.root, 0, 0:d)
        description(canvas, "Bounding Box Tree Overview (all levels)")
        C.show_page(canvas)

        # Subsequent pages: Level-by-level visualization
        for i = (d-1):-1:0
            bbox_draw(canvas, bb_root, Clr.parse(Clr.Colorant, "lightblue"))
            node_draw(canvas, tree.root, 0, i:(i+1))
            description(canvas, "Bounding Box Tree Levels $i to $(i+1)")
            C.show_page(canvas)
        end
    end
    println("Created tree visualization: $filename")
end

end # module
