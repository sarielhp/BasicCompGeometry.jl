module MVBB

using ..BasicCompGeometry
using ..BasicCompGeometry.BBT
using ..BasicCompGeometry.VirtArray
using LinearAlgebra
using StaticArrays
using DataStructures

"""
    OBBox{D, T}

Arbitrarily oriented bounding box in D dimensions.
Represented by an orthonormal basis (axes) and the min/max extents along each axis.
"""
mutable struct OBBox{D,T}
    axes::SVector{D,Point{D,T}}
    low::MVector{D,T}
    high::MVector{D,T}
    f_init::Bool
end

function OBBox{D,T}() where {D,T}
    axes = SVector{D,Point{D,T}}([Point{D,T}(i == j ? 1 : 0 for j = 1:D) for i = 1:D])
    return OBBox{D,T}(axes, zeros(MVector{D,T}), zeros(MVector{D,T}), false)
end

function OBBox(axes::SVector{D,Point{D,T}}) where {D,T}
    return OBBox{D,T}(axes, zeros(MVector{D,T}), zeros(MVector{D,T}), false)
end

function volume(obb::OBBox{D,T}) where {D,T}
    !obb.f_init && return zero(T)
    vol = one(T)
    for i = 1:D
        vol *= (obb.high[i] - obb.low[i])
    end
    return vol
end

function bound!(obb::OBBox{D,T}, pnt::Point{D,T}) where {D,T}
    if !obb.f_init
        for i = 1:D
            val = dot(obb.axes[i], pnt)
            obb.low[i] = val
            obb.high[i] = val
        end
        obb.f_init = true
    else
        for i = 1:D
            val = dot(obb.axes[i], pnt)
            if val < obb.low[i]
                obb.low[i] = val
            end
            if val > obb.high[i]
                obb.high[i] = val
            end
        end
    end
    return obb
end

function bound!(obb::OBBox{D,T}, PS::AbsPntSeq{D,T}) where {D,T}
    for p in Points(PS)
        bound!(obb, p)
    end
    return obb
end

"""
    GPointPair{D, T}

A pair of points and their distance.
"""
mutable struct GPointPair{D,T}
    p::Point{D,T}
    q::Point{D,T}
    distance::T
end

function update_diam!(pp::GPointPair{D,T}, p::Point{D,T}, q::Point{D,T}) where {D,T}
    d = dist(p, q)
    if d > pp.distance
        pp.p = p
        pp.q = q
        pp.distance = d
    end
end

function update_diam_subspace!(pp::GPointPair{D,T}, p::Point{D,T}, q::Point{D,T}, dir::Point{D,T}) where {D,T}
    # Compute distance in the subspace orthogonal to dir
    # Project p-q onto dir
    diff = p - q
    proj_len = dot(diff, dir)
    d_sq = dot(diff, diff) - proj_len^2
    d = d_sq > 0 ? sqrt(d_sq) : zero(T)
    if d > pp.distance
        pp.p = p
        pp.q = q
        pp.distance = d
    end
end

"""
    generate_orthonormal_base(in_vec::Point{3, T})

Generate an orthonormal base (v1, v2, v3) where v1 is parallel to in_vec.
"""
function generate_orthonormal_base(in_vec::Point{3,T}) where {T}
    v1 = normalize(in_vec)
    if abs(v1[1]) < 1e-6 && abs(v1[2]) < 1e-6
        v2 = Point{3,T}(1, 0, 0)
    else
        v2 = normalize(Point{3,T}(-v1[2], v1[1], 0))
    end
    v3 = normalize(v1 × v2)
    return SVector{3,Point{3,T}}(v1, v2, v3)
end

# Rotating Calipers for 2D Minimum Area Rectangle
struct BBox2DInfo{T}
    vertices::SVector{4,Point{2,T}}
    area::T
end

function rotating_calipers_min_area(hull::AbsPntSeq{2,T}) where {T}
    n = cardin(hull)
    if n == 0
        return BBox2DInfo(SVector{4,Point{2,T}}(zeros(Point{2,T}, 4)), zero(T))
    elseif n == 1
        p = hull[1]
        return BBox2DInfo(SVector{4,Point{2,T}}(p, p, p, p), zero(T))
    elseif n == 2
        p1, p2 = hull[1], hull[2]
        return BBox2DInfo(SVector{4,Point{2,T}}(p1, p2, p2, p1), zero(T))
    end

    # Precompute edge angles
    angles = Vector{T}(undef, n)
    for i = 1:n
        p1 = hull[i]
        p2 = hull[i % n + 1]
        angles[i] = atan(p2[2] - p1[2], p2[1] - p1[1])
        if angles[i] < 0
            angles[i] += 2 * pi
        end
    end

    # Rotating calipers
    min_area = Inf
    best_bbox = nothing

    # Helper to find vertex in a certain direction
    function find_vertex(start_idx, target_angle)
        idx = start_idx
        while true
            curr_angle = angles[idx]
            prev_angle = angles[idx == 1 ? n : idx - 1]
            
            # Normalize angles to [0, 2pi) and handle wrap-around
            if prev_angle <= curr_angle
                if prev_angle < target_angle <= curr_angle
                    return idx
                end
            else
                if prev_angle < target_angle || target_angle <= curr_angle
                    return idx
                end
            end
            idx = idx % n + 1
            if idx == start_idx
                break
            end
        end
        return idx
    end

    # Initial directions
    s = 1
    v = 1
    t = 1

    for u = 1:n
        ang0 = angles[u]
        ang1 = mod(ang0 + pi / 2, 2 * pi)
        ang2 = mod(ang0 + pi, 2 * pi)
        ang3 = mod(ang0 + 3 * pi / 2, 2 * pi)

        s = find_vertex(s, ang1)
        v = find_vertex(v, ang2)
        t = find_vertex(t, ang3)

        # Intersection of lines
        function intersect(p1, a1, p2, a2)
            # Line 1: p1 + t1 * (cos a1, sin a1)
            # Line 2: p2 + t2 * (cos a2, sin a2)
            dx = p2[1] - p1[1]
            dy = p2[2] - p1[2]
            det = cos(a1) * sin(a2) - sin(a1) * cos(a2)
            if abs(det) < 1e-12
                return p1 # Parallel
            end
            t1 = (dx * sin(a2) - dy * cos(a2)) / det
            return Point{2,T}(p1[1] + t1 * cos(a1), p1[2] + t1 * sin(a1))
        end

        v0 = intersect(hull[u], ang0, hull[s], ang1)
        v1 = intersect(hull[s], ang1, hull[v], ang2)
        v2 = intersect(hull[v], ang2, hull[t], ang3)
        v3 = intersect(hull[t], ang3, hull[u], ang0)

        area = dist(v0, v1) * dist(v0, v3)
        if area < min_area
            min_area = area
            best_bbox = BBox2DInfo(SVector{4,Point{2,T}}(v0, v1, v2, v3), area)
        end
    end

    return best_bbox
end

"""
    approx_diam(PS::AbsPntSeq{D, T}, ε::Real) where {D, T}

Approximate the diameter of a point set using the Fair Split Tree algorithm.
"""
function approx_diam(PS::AbsPntSeq{D,T}, ε::Real) where {D,T}
    tree = Tree_init(PS)
    # Ensure it's fully expanded if we want to use it like the C++ code
    # though the C++ code splits nodes lazily.
    # For now, let's follow the C++ logic: start with root-root pair.
    
    # Initialize with an arbitrary pair
    pp = GPointPair(PS[1], PS[end], dist(PS[1], PS[end]))
    
    # Priority Queue for pairs of nodes
    # We want to maximize the maximum possible distance between nodes
    # Max possible distance is max_dist(bb1, bb2)
    pq = PriorityQueue{Tuple{Node,Node},T}(Base.Order.Reverse)
    
    root = tree.root
    pq[(root, root)] = max_dist(root.bb, root.bb)
    
    # Initial diameter estimate from longest dimension of bounding box
    dim_max = get_longest_dim(root.bb)
    # Find extreme points in this dimension
    p_min = PS[1]
    p_max = PS[1]
    val_min = p_min[dim_max]
    val_max = p_max[dim_max]
    for p in PS
        if p[dim_max] < val_min
            val_min = p[dim_max]
            p_min = p
        end
        if p[dim_max] > val_max
            val_max = p[dim_max]
            p_max = p
        end
    end
    update_diam!(pp, p_min, p_max)

    while !isempty(pq)
        (n1, n2), d_max = first(pq)
        popfirst!(pq)
        
        if d_max <= (1 + ε) * pp.distance
            break
        end
        
        # If both are leaves, just update diameter
        if n1.f_leaf && n2.f_leaf
            # In a fair split tree, each leaf has at least one point.
            # Usually multiple points if they are identical.
            p1 = tree.PS[first(n1.r)]
            p2 = tree.PS[first(n2.r)]
            update_diam!(pp, p1, p2)
            continue
        end
        
        # Split the pair
        # Ensure children exist
        if !n1.f_leaf && (isnothing(n1.left) || isnothing(n1.right))
            Tree_refine_node(tree, n1)
        end
        if !n2.f_leaf && (isnothing(n2.left) || isnothing(n2.right))
            Tree_refine_node(tree, n2)
        end
        
        # Try to update diameter with some points from the nodes
        p1 = tree.PS[first(n1.r)]
        p2 = tree.PS[first(n2.r)]
        update_diam!(pp, p1, p2)

        if !n1.f_leaf && !n2.f_leaf
            # Both splitable
            for c1 in (n1.left, n1.right), c2 in (n2.left, n2.right)
                if !isnothing(c1) && !isnothing(c2)
                    d = max_dist(c1.bb, c2.bb)
                    if d > (1 + ε) * pp.distance
                        pq[(c1, c2)] = d
                    end
                end
            end
        elseif !n1.f_leaf
            # Only n1 splitable
            for c1 in (n1.left, n1.right)
                if !isnothing(c1)
                    d = max_dist(c1.bb, n2.bb)
                    if d > (1 + ε) * pp.distance
                        pq[(c1, n2)] = d
                    end
                end
            end
        else
            # Only n2 splitable
            for c2 in (n2.left, n2.right)
                if !isnothing(c2)
                    d = max_dist(n1.bb, c2.bb)
                    if d > (1 + ε) * pp.distance
                        pq[(n1, c2)] = d
                    end
                end
            end
        end
    end
    
    return pp
end

function get_longest_dim(bb::BBox{D,T}) where {D,T}
    dim = 1
    max_w = width(bb, 1)
    for i = 2:D
        w = width(bb, i)
        if w > max_w
            max_w = w
            dim = i
        end
    end
    return dim
end

"""
    approx_diam_subspace(PS::AbsPntSeq{D, T}, ε::Real, dir::Point{D, T}) where {D, T}

Approximate the diameter of a point set projected into the subspace orthogonal to `dir`.
"""
function approx_diam_subspace(PS::AbsPntSeq{D,T}, ε::Real, dir::Point{D,T}) where {D,T}
    tree = Tree_init(PS)
    pp = GPointPair(PS[1], PS[end], 0.0)
    update_diam_subspace!(pp, PS[1], PS[end], dir)
    
    pq = PriorityQueue{Tuple{Node,Node},T}(Base.Order.Reverse)
    root = tree.root
    pq[(root, root)] = max_dist_subspace(root.bb, root.bb, dir)
    
    while !isempty(pq)
        (n1, n2), d_max_ub = first(pq)
        popfirst!(pq)
        
        if d_max_ub <= (1 + ε) * pp.distance
            break
        end
        
        if n1.f_leaf && n2.f_leaf
            p1 = tree.PS[first(n1.r)]
            p2 = tree.PS[first(n2.r)]
            update_diam_subspace!(pp, p1, p2, dir)
            continue
        end

        if !n1.f_leaf && (isnothing(n1.left) || isnothing(n1.right))
            Tree_refine_node(tree, n1)
        end
        if !n2.f_leaf && (isnothing(n2.left) || isnothing(n2.right))
            Tree_refine_node(tree, n2)
        end

        # Sample points to update lower bound
        p1 = tree.PS[first(n1.r)]
        p2 = tree.PS[first(n2.r)]
        update_diam_subspace!(pp, p1, p2, dir)

        if !n1.f_leaf && !n2.f_leaf
            for c1 in (n1.left, n1.right), c2 in (n2.left, n2.right)
                if !isnothing(c1) && !isnothing(c2)
                    d = max_dist_subspace(c1.bb, c2.bb, dir)
                    if d > (1 + ε) * pp.distance
                        pq[(c1, c2)] = d
                    end
                end
            end
        elseif !n1.f_leaf
            for c1 in (n1.left, n1.right)
                if !isnothing(c1)
                    d = max_dist_subspace(c1.bb, n2.bb, dir)
                    if d > (1 + ε) * pp.distance
                        pq[(c1, n2)] = d
                    end
                end
            end
        else
            for c2 in (n2.left, n2.right)
                if !isnothing(c2)
                    d = max_dist_subspace(n1.bb, c2.bb, dir)
                    if d > (1 + ε) * pp.distance
                        pq[(n1, c2)] = d
                    end
                end
            end
        end
    end
    return pp
end

"""
    solve_2d_subspace(PS::AbsPntSeq{3, T}, dir::Point{3, T}) where {T}

Find the best OBBox for PS where one of the axes is parallel to `dir`.
"""
function solve_2d_subspace(PS::AbsPntSeq{3,T}, dir::Point{3,T}) where {T}
    norm_dir = normalize(dir)
    base = generate_orthonormal_base(norm_dir)
    # Project points to 2D
    pnts2d = [Point{2,T}(dot(p, base[2]), dot(p, base[3])) for p in Points(PS)]
    hull2d = convex_hull(PntSeq(pnts2d))
    bbox2d = rotating_calipers_min_area(hull2d)
    
    # Vertices of 2D box to get directions
    v0 = bbox2d.vertices[1]
    v1 = bbox2d.vertices[2]
    v3 = bbox2d.vertices[4]
    
    dir2_2d = normalize(v1 - v0)
    dir3_2d = normalize(v3 - v0)
    
    # Lift back to 3D
    dir2 = normalize(dir2_2d[1] * base[2] + dir2_2d[2] * base[3])
    dir3 = normalize(dir3_2d[1] * base[2] + dir3_2d[2] * base[3])
    
    obb = OBBox(SVector{3,Point{3,T}}(norm_dir, dir2, dir3))
    bound!(obb, PS)
    return obb
end

"""
    optimize(PS::AbsPntSeq{3, T}, obb::OBBox{3, T}, iterations::Int) where {T}

Iteratively refine the OBBox by rotating the subspace.
"""
function optimize(PS::AbsPntSeq{3,T}, obb::OBBox{3,T}, iterations::Int) where {T}
    best_obb = obb
    best_vol = volume(best_obb)
    
    for i = 1:iterations
        # Cycle through axes as projection normals
        axis_idx = (i - 1) % 3 + 1
        new_obb = solve_2d_subspace(PS, best_obb.axes[axis_idx])
        new_vol = volume(new_obb)
        
        if new_vol < best_vol
            best_vol = new_vol
            best_obb = new_obb
        end
    end
    return best_obb
end

"""
    approx_mvbb(PS::AbsPntSeq{3, T}, ε::Real) where {T}

Approximate the minimum volume bounding box of a 3D point set.
"""
function approx_mvbb(PS::AbsPntSeq{3,T}, ε::Real) where {T}
    # 1. Start with approximate diameter as first axis
    pp1 = approx_diam(PS, ε)
    dir1 = normalize(pp1.p - pp1.q)
    
    obb_init = solve_2d_subspace(PS, dir1)
    
    # 2. Iteratively refine
    obb_refined = optimize(PS, obb_init, 10)
    
    # 3. Also consider axis-aligned box
    aabb = BBox(PS)
    vol_obb = volume(obb_refined)
    vol_aabb = volume(OBBox(SVector{3,Point{3,T}}(
        Point{3,T}(1,0,0), Point{3,T}(0,1,0), Point{3,T}(0,0,1)),
        MVector{3,T}(aabb.mini), MVector{3,T}(aabb.maxi), true))
        
    if vol_aabb < vol_obb
        return OBBox(SVector{3,Point{3,T}}(
            Point{3,T}(1,0,0), Point{3,T}(0,1,0), Point{3,T}(0,0,1)),
            MVector{3,T}(aabb.mini), MVector{3,T}(aabb.maxi), true)
    end
    
    return obb_refined
end

export OBBox, volume, approx_diam, approx_mvbb, rotating_calipers_min_area

end # module
