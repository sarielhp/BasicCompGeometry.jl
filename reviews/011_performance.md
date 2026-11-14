# Systems Code Review Report #011 (Lens: performance)

- **Date**: 2026-11-14
- **Auditor**: Gemini 3.8 Flash (Tier 0 Workhorse) via `review-audit`
- **Focus Lens**: `performance`
- **Model Tier**: `TIER0`
- **Backend**: `gemini`
- **Scope**: `src/CairoDraw.jl src/ReadWrite.jl ext/CairoExt.jl ext/LaTeXCairoExt.jl ext/BBTCairoExt.jl`
- **Status**: Action Required

---

Auditing via Gemini Flash (bws run) [Profile: performance]...
[SEVERITY]: Moderate
[LOCATION]: [src/ReadWrite.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/ReadWrite.jl#L203-L225) (`_parse_numbers`, `read_plt`)
[ROOT CAUSE]: `read_plt` takes `dim::Int = 2` as a runtime value and passes it directly as a type parameter to `PntSeq{dim,Float64}()` and `Point{dim,Float64}(ntuple(i -> pieces[i], dim))`. Inspecting `BasicCompGeometry.ReadWrite.read_plt("x.plt")` and its keyword body with `@code_warntype` confirms `P::PntSeq{_A, Float64} where _A` (non-concrete `UnionAll`), `ntuple(...)` inferred as `Tuple{Vararg{Float64}}`, `Point{dim,Float64}(...)` inferred as `Any`, and dynamic dispatch on `push_distinct!(P, ...)` on every point. Additionally, `readlines(filename)` eagerly allocates a `Vector{String}` for the whole file, and `_parse_numbers` allocates both a `Vector{SubString{String}}` (from `split`) and a `Vector{Float64}` (parsing all tokens on the line even though only the first `dim` are used) per line.
[FAILURE TRACE]: Calling `read_plt("curve.plt")` on a file with $N$ points materializes all $N$ lines in a `Vector{String}`, allocates a `Vector{SubString{String}}` and `Vector{Float64}` on every non-comment line, boxes each `Point{dim,Float64}`, and performs $N$ dynamic dispatches to `push_distinct!`.
[REMEDIATION]:
```julia
read_plt(filename::AbstractString; dim::Int = 2) = read_plt(filename, Val(dim))

function read_plt(filename::AbstractString, ::Val{D}) where {D}
    P = PntSeq{D,Float64}()
    for line in eachline(filename)
        text = strip(line)
        (isempty(text) || startswith(text, '#') || text == "x y k tid") && continue
        coords = _parse_first_n(text, Val(D), line)
        push_distinct!(P, Point{D,Float64}(coords))
    end
    return P
end

function _parse_first_n(text::AbstractString, ::Val{D}, raw_line::AbstractString) where {D}
    idx = firstindex(text)
    last_idx = lastindex(text)
    coords = ntuple(D) do _
        while idx <= last_idx && (text[idx] == ',' || isspace(text[idx]))
            idx = nextind(text, idx)
        end
        idx > last_idx && error("expected at least $D numbers in line \"$raw_line\"")
        start_idx = idx
        while idx <= last_idx && text[idx] != ',' && !isspace(text[idx])
            idx = nextind(text, idx)
        end
        parse(Float64, SubString(text, start_idx, prevind(text, idx)))
    end
    return coords
end
```

[SEVERITY]: Moderate
[LOCATION]: [src/ReadWrite.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/ReadWrite.jl#L234-L242) (`read_geolife_plt`)
[ROOT CAUSE]: `Iterators.drop(readlines(filename), 6)` eagerly reads the entire `.plt` file into a `Vector{String}` before skipping the first 6 lines, and `split(line, ",")` allocates a `Vector{SubString{String}}` across all comma-separated fields on every line even though only the first two coordinates (`parts[1]` and `parts[2]`) are needed.
[FAILURE TRACE]: Loading large GeoLife `.plt` trajectory files with $N$ records allocates an $N$-element `Vector{String}` upfront plus a 7-element `Vector{SubString{String}}` per data line.
[REMEDIATION]:
```julia
function read_geolife_plt(filename::AbstractString)
    P = PntSeq2F()
    for line in Iterators.drop(eachline(filename), 6)
        isempty(strip(line)) && continue
        c1 = findfirst(',', line)
        c2 = c1 === nothing ? nothing : findnext(',', line, nextind(line, c1))
        (c1 === nothing || c2 === nothing) && error("malformed GeoLife record: \"$line\"")
        x = parse(Float64, SubString(line, firstindex(line), prevind(line, c1)))
        y = parse(Float64, SubString(line, nextind(line, c1), prevind(line, c2)))
        push_distinct!(P, point(x, y))
    end
    return P
end
```

[SEVERITY]: Moderate
[LOCATION]: [src/ReadWrite.jl](file:///home/sariel/prog/26/BasicCompGeometry/src/ReadWrite.jl#L114-L140) (`write_fbin`, `write_ibin`)
[ROOT CAUSE]: `write(f, vecs')` passes a lazy `Adjoint` wrapper (`Adjoint{Float32, Matrix{Float32}}` / `Adjoint{Int32, Matrix{Int32}}`) to `Base.write`. As verified via `@which write(devnull, zeros(Float32, 2, 2)')`, this dispatches to `write(s::IO, A::AbstractArray)` (`Base io.jl:844`), which iterates element-by-element instead of performing a bulk `unsafe_write`, issuing `nvecs * dim` individual scalar `write` calls while reading `vecs` with non-contiguous stride `nvecs`.
[FAILURE TRACE]: Serializing a large matrix (e.g., $10^6 \times 128$ `Float32`) with `write_fbin` executes $1.28 \times 10^8$ individual `write` calls with column-strided memory reads rather than a single contiguous buffer write.
[REMEDIATION]:
```julia
function _write_bin(filename::String, vecs::Matrix{T}) where {T<:Union{Float32,Int32}}
    open(filename, "w") do f
        nvecs, dim = size(vecs)
        write(f, Int32(nvecs))
        write(f, Int32(dim))
        write(f, permutedims(vecs))
    end
end

write_fbin(filename::String, vecs::Matrix{Float32}) = _write_bin(filename, vecs)
write_ibin(filename::String, vecs::Matrix{Int32}) = _write_bin(filename, vecs)
```

[SEVERITY]: Moderate
[LOCATION]: [ext/BBTCairoExt.jl](file:///home/sariel/prog/26/BasicCompGeometry/ext/BBTCairoExt.jl#L35-L78) (`BBT.Tree_draw`, `node_draw`)
[ROOT CAUSE]:
1. **Unpruned `O(d * N)` / `O(N^2)` traversal**: In `for i = (d-1):-1:0`, `node_draw(canvas, tree.root, 0, i:(i+1))` walks the entire `N`-node tree down to every leaf on each page because it continues recursing into `node.left` and `node.right` even when `level >= last(range)`. For a tree of depth `d`, this visits $\Theta(d \cdot N)$ nodes ($\Theta(N^2)$ when unbalanced with $d = \Theta(N)$) instead of $O(N)$ total across all pages.
2. **Type-unstable per-node color parsing**: Inside `if level ∈ range`, `Clr.parse(Clr.Colorant, "yellow")` and `Clr.parse(Clr.Colorant, "blue")` are called for every drawn internal node. Passing the abstract type `Colors.Colorant` and a runtime `String` is type-unstable and forces repeated string dictionary lookups, color allocation, and dynamic dispatch on `Clr.coloralpha` and `C.set_source` on every node.
[FAILURE TRACE]: Calling `BBT.Tree_draw(tree, "tree.pdf")` on a `Tree{2,...}` with $N$ nodes and depth $d$ traverses all $N$ nodes on every one of the $d + 1$ pages and dynamically parses `"yellow"` and `"blue"` as `Colors.Colorant` on every internal node rendered.
[REMEDIATION]:
```julia
    fill_color = Clr.RGBA(1.0, 1.0, 0.0, 0.1)
    stroke_color = Clr.RGB(0.0, 0.0, 1.0)
    bg_color = Clr.RGB(0.6784313725490196, 0.8470588235294118, 0.9019607843137255)

    d = BBT.depth(tree.root)
    NodeT = typeof(tree.root)
    levels = [NodeT[] for _ = 0:d]
    function collect_levels!(node, level)
        node.f_leaf && return
        push!(levels[level + 1], node)
        !isnothing(node.left) && collect_levels!(node.left, level + 1)
        !isnothing(node.right) && collect_levels!(node.right, level + 1)
    end
    collect_levels!(tree.root, 0)

    function draw_node(canvas, node)
        C.set_source(canvas.cr, fill_color)
        bb = node.bb + (diam(node.bb) * 0.01)
        bl = bottom_left(bb)
        C.rectangle(canvas, bl[1], bl[2], width(bb, 1), height(bb, 2))
        C.fill_preserve(canvas)
        C.set_source(canvas.cr, stroke_color)
        cairo_set_line_width(canvas, 1.5)
        C.stroke(canvas)
    end
    # Page 1: iterate levels 0:d; level pages: iterate i:(i+1) using levels[lvl + 1]
```
