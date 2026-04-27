raw"""
    exact_diameter(P)

Calculate the exact Euclidean diameter of a set of points `P` (Polygon or vector of points).
This uses a brute-force $O(n^2)$ approach.
"""
function exact_diameter(P::Polygon{D, T}) where {D, T}
    n = length(P)
    n <= 1 && return 0.0
    curr = 0.0
    for i in 1:n-1
        for j in i+1:n
            curr = max(curr, dist(P[i], P[j]))
        end
    end
    return curr
end

raw"""
    approx_diameter(P, ε)

Calculate a $(1+ε)$-approximation of the Euclidean diameter of a set of points `P`.
Uses Well-Separated Pairs Decomposition (WSPD) to achieve $O(n \log n + n/\epsilon^d)$ time complexity.
"""
function approx_diameter(P::Polygon{D, T}, ε::Real) where {D, T}
    n = length(P)
    n <= 1 && return 0.0
    
    c = 1.0 + ε
    W = WSPD.init(P, ε / 2.0)
    
    # Initialize with distance between first two distinct points
    curr = dist(P[1], P[end])
    
    while !isempty(W.heap)
        if WSPD.top_diam_ub(W) <= c * curr
            break
        end
        
        top = WSPD.get_top(W)
        p, q = WSPD.get_reps(W, top)
        curr = max(curr, dist(p, q))
        WSPD.top_refine!(W)
    end
    
    return curr
end

export exact_diameter, approx_diameter
