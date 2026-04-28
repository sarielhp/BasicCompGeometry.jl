"""
    LongestConvexSubset

This module implements an algorithm to find the largest convex subset of a given set of points.
The algorithm uses dynamic programming (DP) to compute the longest convex chain 
between pairs of points and then combines upper and lower chains to form the longest convex polygon.

The running time of the algorithm is O(n^3 log n).
"""
module LongestConvexSubset

using ..BasicCompGeometry
using Printf

import ..BasicCompGeometry: dist, is_left_turn, is_right_turn, Points

export longest_convex_chain, longest_concave_chain, compute_largest_convex_subset

################################################################################
# Data Structures
################################################################################

"""
    PointTurn

Stores information about possible next points in a convex chain from a given point.
Used to represent the state in the dynamic programming approach.
"""
mutable struct PointTurn
    tail::Vector{Int}
    tvals::Vector{Int}
    jump::Vector{Bool}
    f_terminal::Bool
end

function PointTurn()
    return PointTurn(Int[], Int[], Bool[], false)
end

"""
    LCSContext{T}

Main context for the Longest Convex Subset algorithm.
"""
mutable struct LCSContext{T}
    PS::AbsPolygon{2,T}
    turns::Vector{PointTurn}
    Lens::Vector{Int}
end

struct RotateOrder{T} <: Base.Order.Ordering
    PS::AbsPolygon{2,T}
    base::Point{2,T}
end

function Base.Order.lt(o::RotateOrder, a::Int, b::Int)
    return is_left_turn(o.PS[a], o.base, o.PS[b])
end

################################################################################
# Initialization
################################################################################

function pt_init(PS::AbsPolygon{2,T}, loc::Int) where {T}
    pt = PointTurn()

    for j = (loc+1):length(PS)
        push!(pt.tail, j)
        push!(pt.tvals, -1)
        push!(pt.jump, false)
    end

    ro = RotateOrder(PS, PS[loc])
    sort!(pt.tail, order = ro)

    return pt
end

function sort_points!(ps::Polygon{2,T}) where {T}
    sort!(ps.pnts, lt = (p, q) -> p[1] < q[1] || (p[1] == q[1] && p[2] < q[2]))
end

function lcs_init(PS::AbsPolygon{2,T}) where {T}
    # PS should already be sorted
    pt = [PointTurn() for _ = 1:length(PS)]
    lct = LCSContext{T}(PS, pt, Int[])

    for i = length(PS):-1:1
        lct.turns[i] = pt_init(PS, i)
    end

    return lct
end

################################################################################
# Dynamic Programming Core
################################################################################

function get_max_chain(lct::LCSContext, i::Int, j::Int, target::Int)
    if (target > 0) && (j > target)
        return -1, -1
    end
    if (target > 0) && (j == target)
        return 1, 0
    end

    ps = lct.PS
    p = ps[i]
    q = ps[j]

    curr = lct.turns[j]
    tail = curr.tail
    r_max = length(tail)

    if r_max == 0
        return (target < 0) ? (1, 0) : (-1, -1)
    end

    # Binary search for the first point k in tail such that (p, q, ps[k]) is a right turn.
    # tail is sorted CCW around q.
    if is_right_turn(p, q, ps[tail[1]])
        return curr.tvals[1], 1
    end
    if !is_right_turn(p, q, ps[tail[r_max]])
        return (target < 0) ? (1, 0) : (-1, -1)
    end

    l, r = 1, r_max
    while l < r
        if (l + 1) == r
            return curr.tvals[r], r
        end
        m = (l + r) >> 1
        if is_right_turn(p, q, ps[tail[m]])
            r = m
        else
            l = m
        end
    end
    return -1, -1
end

function fill_vals_for_point!(lct::LCSContext, i::Int, target::Int)
    curr = lct.turns[i]
    if (target > 0) && (i > target)
        return
    end

    len = length(curr.tail)
    if len > 0
        curr.tvals[len] = -1
    end

    for j = len:-1:1
        t_j = curr.tail[j]
        mx, loc = get_max_chain(lct, i, t_j, target)

        next_pos = min(j + 1, len)
        t_prev = curr.tvals[next_pos]

        mx_next::Int = (mx < 0) ? mx : mx + 1

        if mx_next > t_prev
            curr.tvals[j] = mx_next
            curr.jump[j] = true
        else
            curr.tvals[j] = t_prev
            curr.jump[j] = false
        end
    end

    if all(==(false), curr.jump)
        curr.f_terminal = true
    end
end

function compute_dp_table!(lct::LCSContext, target::Int = -1)
    n = length(lct.PS)

    for i = n:-1:1
        fill_vals_for_point!(lct, i, target)
    end

    for i = 1:n
        curr = lct.turns[i]
        val = 1
        if length(curr.tvals) > 0
            val = first(curr.tvals)
        end
        if target > 0
            if i > target
                val = -2 * n
            elseif i == target
                val = 1
            end
        end
        push!(lct.Lens, val)
    end
end

################################################################################
# Solution Extraction
################################################################################

function extract_solution!(
    lct::LCSContext,
    sol::Vector{Int},
    loc::Int,
    prev::Int,
    target::Int = -1,
)
    push!(sol, loc)
    (loc == target) && return

    curr = lct.turns[loc]
    tvals = curr.tvals
    jump = curr.jump

    len_tvals = length(tvals)
    if prev == -1
        (len_tvals == 0) && return
        l = 1
        while (l < len_tvals) && (jump[l] == false)
            l += 1
        end
        extract_solution!(lct, sol, curr.tail[l], loc, target)
        return
    end

    mx, new_loc = get_max_chain(lct, prev, loc, target)
    if new_loc <= 0
        return
    end
    while (new_loc < length(curr.jump)) && (curr.jump[new_loc] == false)
        new_loc += 1
    end
    next_point_index = curr.tail[new_loc]
    extract_solution!(lct, sol, next_point_index, loc, target)
end

################################################################################
# Public API
################################################################################

"""
    longest_convex_chain(ps::Polygon{2, T})

Finds the longest convex chain in a set of points.
"""
function longest_convex_chain(ps::AbsPolygon{2,T}) where {T}
    n = length(ps)
    n == 0 && return Point{2,T}[]
    n < 3 && return ps.pnts

    # Points must be sorted
    ps_sorted = deepcopy(ps)
    sort_points!(ps_sorted)

    lct = lcs_init(ps_sorted)
    compute_dp_table!(lct)

    max_val, start_node = findmax(lct.Lens)
    sol = Int[]
    extract_solution!(lct, sol, start_node, -1)

    return [ps_sorted[i] for i in sol]
end

"""
    longest_concave_chain(ps::Polygon{2, T})

Finds the longest concave chain in a set of points by reflecting on the x-axis.
"""
function longest_concave_chain(ps::AbsPolygon{2,T}) where {T}
    ps_reflected = Polygon{2,T}([point(p[1], -p[2]) for p in ps])
    chain = longest_convex_chain(ps_reflected)
    return [point(p[1], -p[2]) for p in chain]
end

"""
    compute_largest_convex_subset(ps::Polygon{2, T})

Computes the largest convex subset of the given points.
"""
function compute_largest_convex_subset(ps::AbsPolygon{2,T}) where {T}
    n = length(ps)
    if n <= 2
        return deepcopy(ps)
    end

    # Work on a sorted copy
    work_ps = deepcopy(ps)
    sort_points!(work_ps)

    ps_reflected = Polygon{2,T}([point(p[1], -p[2]) for p in work_ps])
    # Note: reflection doesn't change x-order, so ps_reflected is also sorted by x.

    max_sz = -1
    max_sol = Int[]

    # We iterate over every possible rightmost point (target)
    for target = 1:n
        lct_upper = lcs_init(work_ps)
        compute_dp_table!(lct_upper, target)

        lct_lower = lcs_init(ps_reflected)
        compute_dp_table!(lct_lower, target)

        # Find best starting point for this target
        best_start = -1
        max_val = -1
        for i = 1:n
            if lct_upper.Lens[i] >= 0 && lct_lower.Lens[i] >= 0
                comb = lct_upper.Lens[i] + lct_lower.Lens[i]
                if comb > max_val
                    max_val = comb
                    best_start = i
                end
            end
        end

        if best_start != -1
            sol_upper = Int[]
            extract_solution!(lct_upper, sol_upper, best_start, -1, target)

            sol_lower = Int[]
            extract_solution!(lct_lower, sol_lower, best_start, -1, target)

            # Combine them
            if length(sol_lower) > 1
                # Remove start and end points from lower chain to avoid duplicates
                # sol_lower goes from start to target.
                # We want a loop: start -> ... -> target -> ... -> start
                rev_lower = reverse(sol_lower) # target -> ... -> start
                # sol_upper is start -> ... -> target
                # combined is sol_upper then rev_lower[2:end-1]
                combined = [sol_upper; rev_lower[2:(end-1)]]
            else
                combined = sol_upper
            end

            if length(combined) > max_sz
                max_sz = length(combined)
                max_sol = combined
            end
        end
    end

    return Polygon{2,T}([work_ps[i] for i in max_sol])
end

end # module LongestConvexSubset
