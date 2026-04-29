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

    # Utility to draw a single bounding box
    function bbox_draw(context, bb, color)
        C.set_source(context, color)
        bl = bottom_left(bb)
        w = width(bb, 1)
        h = height(bb, 2)
        C.rectangle(context, bl[1], bl[2], w, h)
        C.fill(context)
    end

    # Recursive utility to draw nodes at specific levels
    function node_draw(context, node, level, range)
        node.f_leaf && return
        if level ∈ range
            yellow_transparent = Clr.coloralpha(Clr.parse(Clr.Colorant, "yellow"), 0.1)
            C.set_source(context, yellow_transparent)

            # Expand slightly for better visibility
            bb = node.bb + (diam(node.bb) * 0.01)
            bl = bottom_left(bb)
            w = width(bb, 1)
            h = height(bb, 2)

            C.rectangle(context, bl[1], bl[2], w, h)
            C.fill_preserve(context)

            C.set_source(context, Clr.parse(Clr.Colorant, "blue"))
            C.set_line_width(context, 1.0)
            C.stroke(context)
        end
        if !isnothing(node.left)
            node_draw(context, node.left, level + 1, range)
        end
        if !isnothing(node.right)
            node_draw(context, node.right, level + 1, range)
        end
    end

    # Main drawing logic
    # Expand root BB for margins
    bb_root = tree.root.bb + (diam(tree.root.bb) * 0.2)

    # We use the max coordinates for surface size
    tr = top_right(bb_root)
    surface = C.CairoPDFSurface(filename, tr[1], tr[2])
    context = C.CairoContext(surface)

    # Page 1: Overview
    bbox_draw(context, bb_root, Clr.parse(Clr.Colorant, "lightblue"))
    d = BBT.depth(tree.root)
    node_draw(context, tree.root, 0, 0:d)
    C.show_page(context)

    # Subsequent pages: Level-by-level visualization
    for i = (d-1):-1:0
        bbox_draw(context, bb_root, Clr.parse(Clr.Colorant, "lightblue"))
        node_draw(context, tree.root, 0, i:(i+1))
        C.show_page(context)
    end

    C.finish(surface)
    println("Created tree visualization: $filename")
end

end # module
