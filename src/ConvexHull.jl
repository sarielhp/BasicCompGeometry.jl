raw"""
    convex_hull(points)

Compute the 2D convex hull of a collection of points using the Monotone Chain algorithm.
Returns a PntSeq representing the convex hull. The vertices are returned in 
counter-clockwise order, starting from the point with the minimum x-coordinate.
The algorithm has $O(n \log n)$ time complexity due to the initial sorting step.
"""

function convex_hull(ps::AbsPntSeq{2,T}) where {T}
    return convex_hull(Points(ps))
end

function convex_hull(pnts::Vector{Point{2,T}}) where {T}
    n = length(pnts)
    if n <= 2
        return PntSeq(pnts)
    end

    # Sort points primarily by x-coordinate, secondarily by y-coordinate
    spnts = sort(pnts, by = p -> (p[1], p[2]))

    # Helper to build a half-hull
    function build_half_hull(points_seq)
        hull = Point{2,T}[]
        for p in points_seq
            while length(hull) >= 2 && !is_left_turn(hull[end-1], hull[end], p)
                pop!(hull)
            end
            push!(hull, p)
        end
        return hull
    end

    # Build lower and upper hulls
    lower = build_half_hull(spnts)
    upper = build_half_hull(reverse(spnts))

    # Concatenate hulls, removing the last point of each as it's the start of the other
    # This result is in counter-clockwise order.
    pop!(lower)
    pop!(upper)

    return PntSeq(vcat(lower, upper))
end

export convex_hull
