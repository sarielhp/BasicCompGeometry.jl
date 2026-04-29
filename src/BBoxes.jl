######################################################################
### Bounding box

"""
    BBox{D, T}

Axis parallel bounding box in D dimensions. 
The fields `mini` and `maxi` store the lower and upper bounds of the box for each dimension.
"""
@with_kw_noshow mutable struct BBox{D,T}
    """Whether the box has been initialized with at least one point."""
    f_init::Bool = false
    """Lower bounds vector."""
    mini::MVector{D,T} = zeros(MVector{D,T})
    """Upper bounds vector."""
    maxi::MVector{D,T} = zeros(MVector{D,T})
end

# Constructors
"""
    BBox(p::Point{D, T}, q::Point{D, T})

Construct a bounding box that tightly encloses points `p` and `q`.
"""
BBox(p::Point{D,T}, q::Point{D,T}) where {D,T} = init!(BBox{D,T}(), p, q)

"""
    BBox(P::AbsPntSeq{D, T})

Construct a bounding box that tightly encloses all vertices of point sequence `P`.
"""
BBox(P::AbsPntSeq{D,T}) where {D,T} = bound!(BBox{D,T}(), P)

"""
    width(bb, dim=1)

Return the extent of the bounding box in the specified dimension `dim`.
Defaults to the first dimension (width in 2D).
"""
@inline width(bb::BBox{D,T}, dim::Int = 1) where {D,T} = bb.maxi[dim] - bb.mini[dim]

"""
    height(bb, dim=2)

Return the extent of the bounding box in the specified dimension `dim`.
Defaults to the second dimension (height in 2D).
"""
@inline height(bb::BBox{D,T}, dim::Int = 2) where {D,T} = bb.maxi[dim] - bb.mini[dim]

function Base.show(io::IO, bb::BBox{D,T}) where {D,T}
    if !bb.f_init
        print(io, "BBox{$D,$T} (uninitialized)")
    else
        for i = 1:D
            print(io, " [[", bb.mini[i], "..", bb.maxi[i], "]] ")
        end
    end
end

"""
    middle(bb, d)

Return the coordinate of the center of the bounding box in dimension `d`.
"""
@inline middle(bb::BBox{D,T}, d::Int) where {D,T} = (bb.mini[d] + bb.maxi[d]) / 2.0

"""
    d_min(bb, d)

Return the minimum coordinate of the box in dimension `d`.
"""
@inline d_min(bb::BBox{D,T}, d::Int) where {D,T} = bb.mini[d]

"""
    d_max(bb, d)

Return the maximum coordinate of the box in dimension `d`.
"""
@inline d_max(bb::BBox{D,T}, d::Int) where {D,T} = bb.maxi[d]

"""
    init!(bb, p, q)

Reset and initialize the bounding box `bb` in-place to tightly enclose points `p` and `q`.
Returns the modified `bb`.
"""
function init!(bb::BBox{D,T}, p::Point{D,T}, q::Point{D,T}) where {D,T}
    bb.f_init = true
    bb.mini .= min.(p, q)
    bb.maxi .= max.(p, q)
    return bb
end

"""
    diam(bb)

Return the length of the diagonal of the bounding box.
"""
function diam(bb::BBox{D,T}) where {D,T}
    !bb.f_init && return zero(T)
    return dist(Point{D,T}(bb.mini), Point{D,T}(bb.maxi))
end

"""
    dist_sq(p::Point, bb::BBox)

Return the minimum squared Euclidean distance between point `p` and bounding box `bb`.
Returns 0 if the point is inside the box.
"""
function dist_sq(p::Point{D,T1}, bb::BBox{D,T2}) where {D,T1,T2}
    !bb.f_init && return zero(promote_type(T1, T2))
    s = zero(promote_type(T1, T2))
    for i = 1:D
        d_i = if p[i] > bb.maxi[i]
            p[i] - bb.maxi[i]
        elseif p[i] < bb.mini[i]
            bb.mini[i] - p[i]
        else
            zero(promote_type(T1, T2))
        end
        s += d_i^2
    end
    return s
end

"""
    dist(p::Point, bb::BBox)

Return the minimum Euclidean distance between point `p` and bounding box `bb`.
Returns 0 if the point is inside the box.
"""
@inline dist(p::Point{D,T1}, bb::BBox{D,T2}) where {D,T1,T2} = sqrt(dist_sq(p, bb))

"""
    dist(bb1, bb2)

Return the minimum Euclidean distance between two bounding boxes.
Returns 0 if the boxes intersect.
"""
function dist(b::BBox{D,T}, c::BBox{D,T}) where {D,T}
    (!b.f_init || !c.f_init) && return zero(T)
    s = zero(T)
    for i = 1:D
        # Minimum distance in dimension i
        d_i = if b.mini[i] > c.maxi[i]
            b.mini[i] - c.maxi[i]
        elseif b.maxi[i] < c.mini[i]
            c.mini[i] - b.maxi[i]
        else
            zero(T)
        end
        s += d_i^2
    end
    return sqrt(s)
end

"""
    max_dist(bb1, bb2)

Return the maximum possible Euclidean distance between any point in `bb1` and any point in `bb2`.
This is the distance between the two most distant corners of the boxes.
"""
function max_dist(b::BBox{D,T}, c::BBox{D,T}) where {D,T}
    !b.f_init && return zero(T)
    !c.f_init && return zero(T)
    s = zero(T)
    for i = 1:D
        d_i = max(abs(b.maxi[i] - c.mini[i]), abs(c.maxi[i] - b.mini[i]))
        s += d_i^2
    end
    return sqrt(s)
end

"""
    max_dist_subspace(b, c, dir)

Return an upper bound on the maximum distance between any point in `b` and any point in `c`
after projection into the subspace orthogonal to `dir`.
"""
function max_dist_subspace(b::BBox{D,T}, c::BBox{D,T}, dir::Point{D,T}) where {D,T}
    !b.f_init && return zero(T)
    !c.f_init && return zero(T)
    
    # The maximum distance between two boxes b and c is achieved at their corners.
    # Specifically, it's the maximum value of ||p - q||_subspace where p is a corner of b
    # and q is a corner of c. This is equivalent to maximizing ||d||_subspace where
    # d is in the Minkowski difference box B_diff = [b.mini - c.maxi, b.maxi - c.mini].
    # Since the subspace norm is convex, the maximum is achieved at a vertex of B_diff.
    # There are 2^D such vertices.
    
    # Construct the Minkowski difference box
    diff_box = BBox{D,T}(true, MVector{D,T}(b.mini .- c.maxi), MVector{D,T}(b.maxi .- c.mini))
    
    max_d_sq = 0.0
    for i in 0:(2^D - 1)
        d = corner(diff_box, i)
        
        proj = dot(d, dir)
        d_sq = dot(d, d) - proj^2
        if d_sq > max_d_sq
            max_d_sq = d_sq
        end
    end
    
    return sqrt(max(0.0, max_d_sq))
end

"""
    bottom_left(bb)

Return the `mini` corner of the box as a `Point`.
"""
@inline bottom_left(bb::BBox{D,T}) where {D,T} = Point{D,T}(bb.mini)

"""
    corner(bb, i)

Return the i-th corner of the bounding box `bb`, where `i` is between 0 and 2^D - 1.
The bits of `i` encode which side (min or max) to take for each dimension.
"""
function corner(bb::BBox{D,T}, i::Integer) where {D,T}
    @assert 0 <= i < (1 << D) "Corner index $i must be in the range [0, $( (1 << D) - 1 )] for D=$D"
    return Point{D,T}(ntuple(j -> (i >> (j - 1)) & 1 == 1 ? bb.maxi[j] : bb.mini[j], D))
end

"""
    top_right(bb)

Return the `maxi` corner of the box as a `Point`.
"""
@inline top_right(bb::BBox{D,T}) where {D,T} = Point{D,T}(bb.maxi)

"""
    expand!(bb, factor)

Rescale the bounding box `bb` in-place by `factor` around its center.
`factor=2` doubles the side lengths.
"""
function expand!(bb::BBox{D,T}, factor::Real) where {D,T}
    center = (bb.maxi + bb.mini) / 2.0
    half_size = (bb.maxi - bb.mini) * (factor / 2.0)
    bb.mini .= center - half_size
    bb.maxi .= center + half_size
    return bb
end

"""
    expand_add!(bb, delta)

Padding operation: adds `delta` to each side of the box in-place.
Each dimension's length increases by `2 * delta`.
"""
function expand_add!(bb::BBox{D,T}, delta::Real) where {D,T}
    bb.mini .-= delta
    bb.maxi .+= delta
    return bb
end

"""
    bb + delta

Return a new bounding box expanded by padding `delta` on each side.
"""
Base.:+(b::BBox{D,T}, delta::Real) where {D,T} = expand_add!(deepcopy(b), delta)

"""
    is_inside(p, bb)

Return `true` if the point `p` is contained within or on the boundary of bounding box `bb`.
"""
function is_inside(p::Point{D,S}, bb::BBox{D,T}) where {D,T,S}
    !bb.f_init && return false
    @inbounds for i = 1:D
        if p[i] < bb.mini[i] || p[i] > bb.maxi[i]
            return false
        end
    end
    return true
end

"""
    bound!(bb, pnt::Point)

Update the bounding box `bb` in-place to include point `pnt`.
If the box was uninitialized, it becomes a zero-size box containing only `pnt`.
"""
function bound!(bb::BBox{D,T}, pnt::Point{D,T}) where {D,T}
    if !bb.f_init
        bb.f_init = true
        bb.mini .= pnt
        bb.maxi .= pnt
    else
        bb.mini .= min.(bb.mini, pnt)
        bb.maxi .= max.(bb.maxi, pnt)
    end
    return bb
end

"""
    bound!(bb, collection)

Update the bounding box `bb` in-place to include all points in the collection.
"""
function bound!(bb::BBox{D,T}, P) where {D,T}
    for p in P
        bound!(bb, p)
    end
    return bb
end

"""
    bound!(bb, P::AbsPntSeq)

Update the bounding box `bb` in-place to include all vertices of point sequence `P`.
"""
function bound!(bb::BBox{D,T}, P::AbsPntSeq{D,T}) where {D,T}
    for p in Points(P)
        bound!(bb, p)
    end
    return bb
end

"""
    BBox2F

Alias for a 2D bounding box with `Float64` coordinates.
"""
const BBox2F = BBox{2,Float64}

export BBox, BBox2F
export expand!, expand_add!
export init!, bound!, diam, dist, max_dist, is_inside
export width, height, middle, d_min, d_max
export bottom_left, top_right
