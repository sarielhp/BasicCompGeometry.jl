# Generic immutable range-minimum queries (rightmost minimum on ties).

"""
    RMQIndex

Abstract supertype of the static range-minimum indexes built by [`build_rmq`](@ref).
"""
abstract type RMQIndex end

# Blocked linear RMQ: one monotone-stack bitmask per key plus a sparse table
# over block minima (O(n / log n) blocks, hence O(n) words).
struct BlockRMQ{V,F} <: RMQIndex
    keys::V
    less::F
    block_size::Int
    masks::Vector{UInt}
    levels::Vector{Vector{Int}}
end

# Reference RMQ storing the minimum of every power-of-two interval.
struct SparseRMQ{V,F} <: RMQIndex
    keys::V
    less::F
    levels::Vector{Vector{Int}}
end

const _RMQ_WORD_BITS = 8sizeof(UInt)
_floor_log2(n::Int) = 8sizeof(Int) - leading_zeros(n) - 1

# One comparison suffices: on equivalence prefer the larger original index.
@inline function _better_index(keys, less, a::Int, b::Int)
    lo, hi = minmax(a, b)
    return less(keys[lo], keys[hi]) ? lo : hi
end

function _sparse_levels(keys, less, minima::Vector{Int})
    isempty(minima) && return Vector{Int}[]
    levels = [minima]
    half = 1
    while half <= length(minima) ÷ 2
        previous = last(levels)
        row = [
            _better_index(keys, less, previous[i], previous[i+half]) for
            i = 1:length(minima)-2half+1
        ]
        push!(levels, row)
        half *= 2
    end
    return levels
end

function _block_masks(keys, less, block_size::Int)
    masks = Vector{UInt}(undef, length(keys))
    minima = Int[]
    for first = 1:block_size:length(keys)
        stop = first + min(block_size, length(keys) - first + 1) - 1
        mask = zero(UInt)
        for j = first:stop
            while !iszero(mask)
                top = _RMQ_WORD_BITS - leading_zeros(mask) - 1
                less(keys[first+top], keys[j]) && break
                mask &= ~(one(UInt) << top)
            end
            mask |= one(UInt) << (j - first)
            masks[j] = mask
        end
        push!(minima, first + trailing_zeros(mask))
    end
    return masks, minima
end

"""
    build_rmq(keys; less=isless, backend=:blocked)

Index the one-based vector `keys` for range-minimum queries under the strict
weak order `less`. Query with [`range_argmin`](@ref). Empty input may be indexed
but has no valid query.

The default `:blocked` backend takes O(n) preprocessing time and machine words
and answers queries in O(1) comparisons and word operations. The model is a
word RAM with constant-time `leading_zeros`, `trailing_zeros` and shifts on
`UInt` words (block size is `floor(log2 n) <= 63 <= word size`); key
comparisons are charged to `less`. No hashing is used. `backend=:sparse` is the
O(n log n) preprocessing/space reference with O(1) queries.

The keys are borrowed, not copied: the caller must keep them (and the order
`less`) unchanged for the lifetime of the index. The index is never mutated by
queries.
"""
function build_rmq(keys::AbstractVector; less = isless, backend::Symbol = :blocked)
    Base.require_one_based_indexing(keys)
    backend in (:blocked, :sparse) ||
        throw(ArgumentError("unknown RMQ backend: $backend"))
    if backend == :sparse
        return SparseRMQ(keys, less, _sparse_levels(keys, less, collect(eachindex(keys))))
    end
    block_size = max(1, _floor_log2(max(1, length(keys))))
    masks, minima = _block_masks(keys, less, block_size)
    return BlockRMQ(keys, less, block_size, masks, _sparse_levels(keys, less, minima))
end

Base.length(index::RMQIndex) = length(index.keys)

function _check_interval(index::RMQIndex, lo::Integer, hi::Integer)
    1 <= lo <= hi <= length(index.keys) || throw(BoundsError(index.keys, lo:hi))
    return nothing
end

function _sparse_minimum(index::RMQIndex, lo::Int, hi::Int)
    power = _floor_log2(hi - lo + 1)
    row = index.levels[power+1]
    return _better_index(index.keys, index.less, row[lo], row[hi-(1<<power)+1])
end

# Precondition: nonempty valid interval wholly inside one block.
function _small_minimum(index::BlockRMQ, lo::Int, hi::Int)
    offset = (lo - 1) % index.block_size
    mask = index.masks[hi] & (typemax(UInt) << offset)
    return ((hi - 1) ÷ index.block_size) * index.block_size + trailing_zeros(mask) + 1
end

"""
    range_argmin(index, lo, hi)

Original index of the minimum key in the inclusive interval `lo:hi`, with
`1 <= lo <= hi <= n`; among equivalent minima (neither key less than the other)
the greatest index is returned. Throws `BoundsError` for an invalid interval.
O(1) time.
"""
function range_argmin(index::BlockRMQ, lo::Integer, hi::Integer)
    _check_interval(index, lo, hi)
    lo, hi = Int(lo), Int(hi)
    b = index.block_size
    left = (lo - 1) ÷ b + 1
    right = (hi - 1) ÷ b + 1
    left == right && return _small_minimum(index, lo, hi)
    a = _small_minimum(index, lo, left * b)
    c = _small_minimum(index, (right - 1) * b + 1, hi)
    best = _better_index(index.keys, index.less, a, c)
    if left + 1 <= right - 1
        middle = _sparse_minimum(index, left + 1, right - 1)
        best = _better_index(index.keys, index.less, best, middle)
    end
    return best
end

function range_argmin(index::SparseRMQ, lo::Integer, hi::Integer)
    _check_interval(index, lo, hi)
    return _sparse_minimum(index, Int(lo), Int(hi))
end
