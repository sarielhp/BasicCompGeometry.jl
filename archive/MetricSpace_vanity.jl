# Historical fragment from src/MetricSpace.jl at
# f7063bef584d6dfa46ca41cbf58ba55b021814e0; removed before v0.5.0.
# Originally defined inside MetricSpace, with dist and update_distances in scope.
# Preserved unchanged for reference; this file is not included by the package.

function get_points_assigned_to_pos(M, I, D, N, pos, sz)
    J = Int[]
    x = I[pos]
    for i ∈ 1:sz
        y = I[i]
        d = dist(M, x, y)
        if d < D[i] || D[i] < 0.0
            push!(J, i)
        end
    end
    return J
end

"""
    greedy_permutation_vanity(M::AbsFMS, vanity::Vector{Float64}, n::Int)

A variant of the greedy permutation where ties (or near-ties) in distance
are broken by minimizing a "vanity" score.

The algorithm first identifies the point furthest from current centers,
finds all points currently "assigned" to that furthest point (i.e., would
be moved to the center set if it were chosen), and then picks the point
within that set that has the minimum vanity score.
"""
function greedy_permutation_vanity(M::AbsFMS, vanity::Vector{Float64}, n::Int)
    sz = size(M)
    I = collect(1:sz)
    N = fill(1, sz)
    D = fill(-1.0, sz)

    update_distances(M, I, D, N, 1, sz)

    for i ∈ 2:n
        tpos = argmax(view(D, i:sz)) + i - 1
        J = get_points_assigned_to_pos(M, I, D, N, tpos, sz)
        j_min = argmin(vanity[j] for j ∈ J)
        pos = J[j_min]

        D[i-1] = D[pos]
        I[i], I[pos] = I[pos], I[i]
        D[i], D[pos] = D[pos], D[i]
        vanity[i], vanity[pos] = vanity[pos], vanity[i]
    end

    if n == sz
        D[sz] = 0.0
    end

    return I, D, N
end
