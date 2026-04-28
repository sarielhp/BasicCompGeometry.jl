# Contains code to do hausdorff simplification of a curve

"""
    hausdorff_dist_subseg(P, rng)

Calculate the Hausdorff distance between the sub-point sequence `P[rng]` and the single segment
connecting its first and last vertices. 
This is used as an error measure for curve simplification.
"""
function hausdorff_dist_subseg(P::AbsPntSeq{D,T}, rng::UnitRange{Int} = 0:0)::T where {D,T}
    if rng == 0:0
        rng = 1:length(P)
    end
    n = length(rng)
    if n <= 2 || length(P) <= 2
        return zero(T)
    end

    s, t = P[first(rng)], P[last(rng)]
    leash = zero(T)
    for i = (first(rng)+1):(last(rng)-1)
        leash = max(leash, dist_point_segment(P[i], s, t))
    end

    return leash
end

"""
    exp_search_prefix(P, start, w)

Find a prefix of pnt_seq `P` starting at `start` that (approximately) exceeds 
the error threshold `w`. Uses exponential search for efficiency on large curves.
Returns the end index of the prefix.
"""
function exp_search_prefix(P::AbsPntSeq{D,T}, start::Int, w::Real) where {D,T}
    n = length(P)
    hi = min(start + 2, n)
    hi >= n && return hi

    while hi < n
        r = hausdorff_dist_subseg(P, start:hi)
        r > w && return min(hi + 5, n)
        hi = start + 2 * (hi - start)
    end
    return n
end

"""
    h_bin_search_inner(P, start, i, j, w)

Perform a binary search between indices `i` and `j` to find the longest prefix 
starting at `start` that has a Hausdorff distance at most `w`.
"""
function h_bin_search_inner(P::AbsPntSeq{D,T}, start::Int, i::Int, j::Int, w::T) where {D,T}
    if i >= j || (start + 1) == j
        return j
    end
    r = hausdorff_dist_subseg(P, start:j)
    r <= w && return j
    if (i + 1) == j
        return i
    end

    mid = (i + j - 1) >> 1
    r_m = hausdorff_dist_subseg(P, start:mid)
    if r_m > w
        return h_bin_search_inner(P, start, i, mid - 1, w)
    end
    return h_bin_search_inner(P, start, mid, j - 1, w)
end

"""
    find_prefix(P, i, j, w)

High-level wrapper to find the optimal split point for simplification between `i` and `j`.
"""
function find_prefix(P::AbsPntSeq{D,T}, i::Int, j::Int, w::T) where {D,T}
    r = hausdorff_dist_subseg(P, i:j)
    if r <= w
        return j
    end
    return h_bin_search_inner(P, i, i, j, w)
end

"""
    hausdorff_simplify(P, w)

Simplify the pnt_seq `P` using a greedy approach based on the Hausdorff distance.
The resulting pnt_seq vertices are a subset of the original vertices. 
For every segment `[v_i, v_{i+1}]` in the simplification, the Hausdorff distance 
to the original subcurve it replaces is at most `w`.
Returns `(simplified_pnt_seq, vertex_indices)`.
"""
function hausdorff_simplify(P::AbsPntSeq{D,T}, w::T) where {D,T}
    pout = PntSeq{D,T}()
    pindices = Int[]
    n = length(P)
    n == 0 && return pout, pindices

    push!(pout, P[1])
    push!(pindices, 1)

    curr_ind = 1
    while true
        hi = exp_search_prefix(P, curr_ind, w)
        next_ind = find_prefix(P, curr_ind, hi, w)
        @assert next_ind > curr_ind
        push!(pindices, next_ind)
        push!(pout, P[next_ind])
        if next_ind == n
            return pout, pindices
        end
        curr_ind = next_ind
    end
end

export hausdorff_simplify, hausdorff_dist_subseg
