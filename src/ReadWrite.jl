"""
    ReadWrite

A module for reading and writing various binary vector formats common in 
computational geometry and similarity search.
Supports .fbin, .ibin, and .fvecs formats.
"""
module ReadWrite

using LinearAlgebra

export read_fbin, read_fbin_n, read_ibin
export write_fbin, write_ibin
export read_fvecs

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

"""
    write_fbin(filename::String, vecs::Matrix{Float32})

Write a matrix of Float32 vectors to an .fbin file. 
`vecs` is expected to have vectors as rows (to match `read_fbin`'s metadata logic).
"""
function write_fbin(filename::String, vecs::Matrix{Float32})
    f = open(filename, "w")
    nvecs, dim = size(vecs)

    write(f, Int32(nvecs))
    write(f, Int32(dim))
    write(f, vecs')

    close(f)
end

"""
    write_ibin(filename::String, vecs::Matrix{Int32})

Write a matrix of Int32 vectors to an .ibin file.
`vecs` is expected to have vectors as rows.
"""
function write_ibin(filename::String, vecs::Matrix{Int32})
    f = open(filename, "w")
    nvecs, dim = size(vecs)

    write(f, Int32(nvecs))
    write(f, Int32(dim))
    write(f, vecs')

    close(f)
end

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

        for i = 1:num_vectors
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

end # module ReadWrite
