# Systems Code Review Report #010 (Lens: performance)

- **Date**: 2026-11-14
- **Auditor**: Claude / Codex (Tier 1 Standard) via `review-audit`
- **Focus Lens**: `performance`
- **Model Tier**: `STANDARD`
- **Backend**: `auto`
- **Scope**: `src/IpeDraw.jl src/IpeDraw/diagrams.jl src/IpeDraw/interface.jl src/IpeDraw/output.jl src/IpeDraw/primitives.jl`
- **Status**: Action Required

---

Auditing via Codex CLI (gpt-5.6-sol, tier: standard) [Profile: performance]...
[SEVERITY]: Moderate  
[LOCATION]: `src/IpeDraw.jl:154-159`, `_apply_tf`  
[ROOT CAUSE]: `_apply_tf` returns the original coordinate type when no viewport is installed, but returns `Point{2,Float64}` when a viewport exists. For non-`Float64` inputs, its inferred return is therefore a value-dependent union such as `Union{Point{2,Int},Point{2,Float64}}`. This instability propagates into polygon comprehensions and formatting loops.  
[FAILURE TRACE]: Draw `Point(1, 2)` or a polygon of integer points. With `canvas.viewport === nothing`, transformed points remain integer-valued; after `fit!`, the same call returns floating-point coordinates. Large integer or rational point collections consequently execute union-split formatting paths and can produce union-typed temporary arrays.  
[REMEDIATION]:
```julia
function _apply_tf(canvas::IpeCanvas, p::Point{2})
    x = Float64(p.x)
    y = Float64(p.y)
    vp = canvas.viewport
    vp === nothing && return Point(x, y)

    yscale = vp.flip_y ? -vp.scale : vp.scale
    return Point(vp.tx + vp.scale * x, vp.ty + yscale * y)
end
```

[SEVERITY]: Moderate  
[LOCATION]: `src/IpeDraw.jl:50-52`, `src/IpeDraw.jl:107`, and `src/IpeDraw/interface.jl:14-16`  
[ROOT CAUSE]: `Style.values::NamedTuple` is abstractly typed: `NamedTuple` omits both its key tuple and value tuple parameters. Every drawing call reads this field and then performs `merge`, iteration, filtering, and keyword splatting on a value whose concrete layout inference cannot recover from the field declaration. `canvas.active_style.values` makes this cost apply even when no explicit style is supplied.  
[FAILURE TRACE]: Draw thousands of points or segments through `draw!`. Each element calls `_style_options`; retrieving either `style.values` or `canvas.active_style.values` loses the concrete named-tuple type, inhibiting specialization of the merge/filter/keyword path and introducing dynamic dispatch and boxing around a per-element operation.  
[REMEDIATION]: Parameterize reusable styles, and keep the active style typed at the scope executing the drawing callback rather than storing differently shaped named tuples behind an abstract field:
```julia
struct Style{NT<:NamedTuple}
    values::NT
end

Style(values::NT) where {NT<:NamedTuple} = Style{NT}(values)
Style(; kwargs...) = Style((; kwargs...))

function _style_options(
    active::Style{A},
    style::Union{Nothing,Style{S}},
    overrides::O,
    allowed,
) where {A,S,O<:NamedTuple}
    local_values = style === nothing ? NamedTuple() : style.values
    merged = merge(active.values, local_values, overrides)
    return (; (k => v for (k, v) in pairs(merged) if k in allowed)...)
end
```
The canvas/style-scope interface should pass the concrete active `Style{A}` through a function barrier; merely changing the field to `active_style::Style` would reintroduce the same abstract-field problem.

[SEVERITY]: Moderate  
[LOCATION]: `src/IpeDraw/primitives.jl:125-149`, `src/IpeDraw/primitives.jl:364-399`, and `src/IpeDraw/primitives.jl:420-444`  
[ROOT CAUSE]: Path serialization materializes complete intermediate collections before creating the final XML. `draw_polygon!` allocates a transformed-point vector. `draw_spline!` allocates a growing `Vector{String}` and one formatted string per segment. `draw_bspline!` allocates transformed points, a `tf_pts[2:end]` slice, formatted-string arrays, joined strings, and the final joined path. Peak memory is multiple times the unavoidable XML output size.  
[FAILURE TRACE]: Render a polygon, spline, or B-spline containing hundreds of thousands of vertices/control points, or render many medium-sized paths in animation frames. Each call retains all transformed points and per-point strings until the final `join`, producing substantial garbage-collector pressure and memory spikes despite serialization requiring only a single pass.  
[REMEDIATION]: Stream transformations and formatting directly into one `IOBuffer`; do not construct point slices or string arrays:
```julia
function draw_bspline!(canvas::IpeCanvas, pts; closed::Bool=false, kwargs...)
    n = length(pts)
    n >= 3 || throw(ArgumentError("B-spline requires at least 3 control points"))

    buf = IOBuffer()
    println(buf, "<path layer=\"$(canvas.active_layer)\">")

    first_index = firstindex(pts)
    if !closed
        p = _apply_tf(canvas, pts[first_index])
        @printf(buf, "%.3f %.3f m\n", p.x, p.y)
    end

    start = closed ? first_index : nextind(pts, first_index)
    for i in start:lastindex(pts)
        p = _apply_tf(canvas, pts[i])
        @printf(buf, "%.3f %.3f\n", p.x, p.y)
    end
    println(buf, closed ? "u" : "s")
    print(buf, "</path>")

    push!(canvas.elements, String(take!(buf)))
    return canvas
end
```
The same direct-to-buffer pattern applies to polygons and cubic splines.

[SEVERITY]: Moderate  
[LOCATION]: `src/IpeDraw/diagrams.jl:5-10` and `src/IpeDraw/diagrams.jl:32-45`  
[ROOT CAUSE]: Every anchored label constructs and populates a new nine-entry `Dict`. XML escaping then performs three successive `replace` calls, each potentially allocating and copying the entire label. These costs occur per label, even though the anchor mapping is immutable.  
[FAILURE TRACE]: Annotate a dense arrangement with 10,000 labels using `anchor=:north` and text containing XML-sensitive characters. This creates 10,000 dictionaries and up to three full-size escaped-string intermediates before constructing each final XML element.  
[REMEDIATION]:
```julia
const _LABEL_ALIGNMENTS = Dict{Symbol,Tuple{Symbol,Symbol}}(
    :center => (:center, :center),
    :north => (:center, :top),
    :south => (:center, :bottom),
    :east => (:right, :center),
    :west => (:left, :center),
    :northeast => (:right, :top),
    :northwest => (:left, :top),
    :southeast => (:right, :bottom),
    :southwest => (:left, :bottom),
)

function _anchor_alignment(anchor::Symbol)
    alignment = get(_LABEL_ALIGNMENTS, anchor, nothing)
    alignment === nothing &&
        throw(ArgumentError("unknown label anchor: $anchor"))
    return alignment
end

# In draw_label!:
if anchor !== nothing
    halign, valign = _anchor_alignment(anchor)
end
```
XML escaping should likewise write escaped characters directly into the final label `IOBuffer`, avoiding the three complete intermediate strings.
