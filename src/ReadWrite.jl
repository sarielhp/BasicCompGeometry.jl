"""
    ReadWrite

A module for reading and writing various binary vector formats common in 
computational geometry and similarity search.
Supports .fbin, .ibin, and .fvecs formats, and text files of planar curves
(plain coordinate lists and GeoLife `.plt` trajectories).
"""
module ReadWrite

using LinearAlgebra
import ..BasicCompGeometry: MVector, Point, PntSeq, PntSeq2F, point, push_distinct!

export read_fbin, read_fbin_n, read_ibin
export write_fbin, write_ibin
export read_fvecs
export read_plt, read_geolife_plt

"""
    read_fbin(filename::String, start_idx::Int = 0, chunk_size::Union{Int,Nothing} = nothing)

Read an .fbin file containing Float32 vectors.
The format starts with two Int32 values: number of vectors and their dimension.

# Arguments
- `filename`: path to .fbin file.
- `start_idx`: start reading vectors from this index (0-based).
- `chunk_size`: number of vectors to read. If nothing, reads until the end.
"""
function read_fbin(
    filename::String,
    start_idx::Int = 0,
    chunk_size::Union{Int,Nothing} = nothing,
)
    f = open(filename, "r")
    nvecs = read(f, Int32)
    dim = read(f, Int32)

    seek(f, 8 + start_idx * sizeof(Float32) * dim)

    nvecs_to_read = isnothing(chunk_size) ? nvecs - start_idx : chunk_size

    arr = Vector{Float32}(undef, nvecs_to_read * dim)
    readbytes!(f, reinterpret(UInt8, arr))

    close(f)

    m = reshape(arr, (dim, nvecs_to_read))

    println("Rows : ", size(m, 1))
    println("Cols : ", size(m, 2))

    return m
end

"""
    read_fbin_n(filename::String, n::Int = 0)

Read at most `n` Float32 vectors from an .fbin file. If `n` is 0, reads all.
"""
function read_fbin_n(filename::String, n::Int = 0)
    f = open(filename, "r")
    nvecs = read(f, Int32)
    dim = read(f, Int32)

    if n != 0
        nvecs = min(nvecs, n)
    end
    seek(f, 8)

    arr = Vector{Float32}(undef, nvecs * dim)
    readbytes!(f, reinterpret(UInt8, arr))

    close(f)
    m = reshape(arr, (dim, nvecs))

    return m
end

"""
    read_ibin(filename::String, start_idx::Int = 0, chunk_size::Union{Int,Nothing} = nothing)

Read an .ibin file containing Int32 vectors. 
The format starts with two Int32 values: number of vectors and their dimension.

Returns a matrix where each row is a vector (consistent with original implementation).
"""
function read_ibin(
    filename::String,
    start_idx::Int = 0,
    chunk_size::Union{Int,Nothing} = nothing,
)
    f = open(filename, "r")
    nvecs = read(f, Int32)
    dim = read(f, Int32)

    seek(f, 8 + start_idx * sizeof(Int32) * dim)

    nvecs_to_read = isnothing(chunk_size) ? nvecs - start_idx : chunk_size

    arr = Vector{Int32}(undef, nvecs_to_read * dim)
    readbytes!(f, reinterpret(UInt8, arr))

    close(f)
    return reshape(arr, (dim, nvecs_to_read))'
end

function _write_bin(filename::String, vecs::Matrix{T}) where {T<:Union{Float32,Int32}}
    open(filename, "w") do f
        nvecs, dim = size(vecs)
        write(f, Int32(nvecs))
        write(f, Int32(dim))
        (nvecs == 0 || dim == 0) && return
        # Transpose into a <=1 MiB chunk buffer so bulk `unsafe_write` avoids scalar IO and full-matrix duplication.
        chunk_rows = clamp(div(262144, dim), 1, nvecs)
        buf = Matrix{T}(undef, dim, chunk_rows)
        r = 1
        while r <= nvecs
            count = min(chunk_rows, nvecs - r + 1)
            @inbounds for j in 1:dim, i in 1:count
                buf[j, i] = vecs[r+i-1, j]
            end
            write(f, @view(buf[:, 1:count]))
            r += count
        end
        return
    end
end

"""
    write_fbin(filename::String, vecs::Matrix{Float32})

Write a matrix of Float32 vectors to an .fbin file. 
`vecs` is expected to have vectors as rows (to match `read_fbin`'s metadata logic).
"""
write_fbin(filename::String, vecs::Matrix{Float32}) = _write_bin(filename, vecs)

"""
    write_ibin(filename::String, vecs::Matrix{Int32})

Write a matrix of Int32 vectors to an .ibin file.
`vecs` is expected to have vectors as rows.
"""
write_ibin(filename::String, vecs::Matrix{Int32}) = _write_bin(filename, vecs)

"""
    read_fvecs(filename::String)

Reads a binary .fvec file and returns a matrix where each column is a vector.
The .fvec format consists of a 4-byte integer for the number of dimensions,
followed by the single-precision floating-point vector data. This pattern
repeats for each vector in the file.
"""
function read_fvecs(filename::String)
    try
        io = open(filename, "r")
        # Read the number of dimensions from the first 4 bytes (as an Int32)
        if eof(io)
            close(io)
            throw(ErrorException("File is empty or malformed."))
        end
        d = read(io, Int32)

        if d <= 0
            close(io)
            throw(ErrorException("Invalid vector dimension read from file: d = $d."))
        end

        vec_record_size_bytes = sizeof(Int32) + d * sizeof(Float32)

        seekend(io)
        file_size_bytes = position(io)

        if file_size_bytes % vec_record_size_bytes != 0
            close(io)
            throw(
                ErrorException(
                    "File size is not a multiple of the vector record size, indicating a malformed file.",
                ),
            )
        end

        seekstart(io)
        num_vectors = div(file_size_bytes, vec_record_size_bytes)
        data = zeros(Float32, d, num_vectors)

        for i in 1:num_vectors
            current_d = read(io, Int32)
            if current_d != d
                close(io)
                throw(
                    ErrorException(
                        "Inconsistent vector dimension detected at vector $i: expected $d, but got $current_d.",
                    ),
                )
            end
            read!(io, view(data, :, i))
        end

        close(io)
        return data
    catch e
        rethrow(ErrorException("Error reading .fvec file: $(e)"))
    end
end

@inline function _parse_next_float(
    text::AbstractString,
    idx::Int,
    last_idx::Int,
    dim::Int,
    raw_line::AbstractString,
)
    while idx <= last_idx && (text[idx] == ',' || isspace(text[idx]))
        idx = nextind(text, idx)
    end
    idx > last_idx && error("expected at least $dim numbers in line \"$raw_line\"")
    start_idx = idx
    while idx <= last_idx && text[idx] != ',' && !isspace(text[idx])
        idx = nextind(text, idx)
    end
    val = parse(Float64, SubString(text, start_idx, prevind(text, idx)))
    return val, idx
end

function _parse_point(text::AbstractString, ::Val{D}, raw_line::AbstractString) where {D}
    buf = MVector{D,Float64}(undef)
    idx = firstindex(text)
    last_idx = lastindex(text)
    @inbounds for i in 1:D
        buf[i], idx = _parse_next_float(text, idx, last_idx, D, raw_line)
    end
    return Point{D,Float64}(buf)
end

function _read_plt(filename::AbstractString, ::Val{D}) where {D}
    P = PntSeq{D,Float64}()
    # open/do closes the file even when a malformed line throws; a bare
    # eachline(filename) would keep it open until GC (Windows cannot delete it).
    open(filename) do io
        for line in eachline(io)
            text = strip(line)
            (isempty(text) || startswith(text, '#') || text == "x y k tid") && continue
            push_distinct!(P, _parse_point(text, Val(D), line))
        end
    end
    return P
end

"""
    read_plt(filename; dim = 2)

Read a point sequence from a text file with one point per line, as written by
`write_plt`: coordinates separated by commas and/or whitespace. The first `dim`
numbers of each line are used. Blank lines, lines starting with `#`, and the
header `x y k tid` (used by some datasets) are skipped, and consecutive repeated
points are dropped.
"""
Base.@constprop :aggressive function read_plt(filename::AbstractString; dim::Int = 2)
    dim < 1 && throw(ArgumentError("dim must be >= 1, got $dim"))
    return _read_plt(filename, Val(dim))
end

"""
    read_geolife_plt(filename)

Read a GeoLife `.plt` trajectory: six header lines, then comma separated records
whose first two fields are the coordinates. Blank lines are skipped and
consecutive repeated points are dropped.
"""
function read_geolife_plt(filename::AbstractString)
    P = PntSeq2F()
    open(filename) do io
        _read_geolife_records!(P, io)
    end
    return P
end

function _read_geolife_records!(P::PntSeq2F, io::IO)
    for line in Iterators.drop(eachline(io), 6)
        isempty(strip(line)) && continue
        c1 = findfirst(',', line)
        c1 === nothing && error("malformed GeoLife record: \"$line\"")
        y_start = nextind(line, c1)
        c2 = findnext(',', line, y_start)
        y_end = c2 === nothing ? lastindex(line) : prevind(line, c2)
        x = parse(Float64, SubString(line, firstindex(line), prevind(line, c1)))
        y = parse(Float64, SubString(line, y_start, y_end))
        push_distinct!(P, point(x, y))
    end
    return P
end

end # module ReadWrite
