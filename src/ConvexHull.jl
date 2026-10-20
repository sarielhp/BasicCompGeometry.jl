raw"""
    convex_hull(points; kernel=FloatingKernel())

Compute the 2D convex hull of a collection of points using the Monotone Chain algorithm.
Returns a PntSeq representing the convex hull. The vertices are returned in 
counter-clockwise order, starting from the point with the minimum x-coordinate.
The algorithm has $O(n \log n)$ time complexity due to the initial sorting step.
The optional kernel is passed to both coordinate comparisons and turn predicates.
`ExactKernel()` computes an exact hull of the supplied coordinates without
changing the coordinate type of the returned vertices.
"""

function convex_hull(
    ps::AbsPntSeq{2,T};
    kernel::AbstractKernel = FloatingKernel(),
) where {T}
    return convex_hull(Points(ps); kernel = kernel)
end

function convex_hull(
    pnts::Vector{Point{2,T}};
    kernel::AbstractKernel = FloatingKernel(),
) where {T}
    return _convex_hull(kernel, pnts)
end

function _convex_hull(kernel::K, pnts::Vector{Point{2,T}}) where {K<:AbstractKernel,T}
    n = length(pnts)
    if n <= 2
        return PntSeq(pnts)
    end

    # Sort points primarily by x-coordinate, secondarily by y-coordinate
    spnts = sort(pnts, lt = (p, q) -> compare_xy(kernel, p, q) < 0)

    # Helper to build a half-hull
    function build_half_hull(points_seq)
        hull = Point{2,T}[]
        for p in points_seq
            while length(hull) >= 2 && !is_left_turn(kernel, hull[end-1], hull[end], p)
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
