"""
    VirtArray

Permutation views over a vector (`VArray`): algorithms reorder the view by swapping
entries of a permutation instead of moving the data, and `orig_index` maps a view
position back to the index in the underlying vector.
"""
module VirtArray

"""
    VArray{T, V}

A virtual array that stores a reference to an underlying vector `vec` of type `V` and 
a permutation vector `loc`. Accessing `va[i]` returns `vec[loc[i]]`.
This is useful for algorithms that need to reorder elements without 
modifying the original data or making expensive copies.
"""
struct VArray{T,V<:AbstractVector{T}} <: AbstractVector{T}
    vec::V
    loc::Vector{Int}
end

"""
    VArray(vec)

Construct a `VArray` wrapping the given vector with an identity permutation.
"""
function VArray(vec::V) where {T,V<:AbstractVector{T}}
    return VArray{T,V}(vec, collect(1:length(vec)))
end

@inline Base.size(va::VArray) = size(va.loc)
@inline Base.length(va::VArray) = length(va.loc)
@inline Base.eachindex(va::VArray) = eachindex(va.loc)

@inline function Base.first(va::VArray)
    @inbounds return va.vec[va.loc[1]]
end

@inline function Base.last(va::VArray)
    @inbounds return va.vec[va.loc[end]]
end

@inline Base.@propagate_inbounds function Base.getindex(va::VArray, i::Integer)
    @boundscheck checkbounds(va, i)
    @inbounds return va.vec[va.loc[i]]
end

"""
    orig_index(va::VArray, i) -> Int

Index in the underlying vector of the element at position `i` of the view.
"""
@inline Base.@propagate_inbounds function orig_index(va::VArray, i::Integer)
    @boundscheck checkbounds(va, i)
    @inbounds return va.loc[i]
end

"""
    swap!(va::VArray, i, j) -> va

Exchange positions `i` and `j` of the view by swapping permutation entries; the
underlying vector is not modified.
"""
@inline Base.@propagate_inbounds function swap!(va::VArray, i::Integer, j::Integer)
    @boundscheck checkbounds(va, i)
    @boundscheck checkbounds(va, j)
    @inbounds va.loc[i], va.loc[j] = va.loc[j], va.loc[i]
    return va
end

@inline function Base.iterate(va::VArray, state::Int = 1)
    if state > length(va.loc)
        return nothing
    else
        @inbounds return va.vec[va.loc[state]], state + 1
    end
end

export VArray, swap!, orig_index

end # module VirtArray
