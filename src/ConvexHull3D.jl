"""
    ConvexHull3D

This module implements the Randomized Incremental Algorithm for computing the 
convex hull of a set of points in 3D. 

The algorithm has an expected running time of O(n log n).
"""
module ConvexHull3D

using ..BasicCompGeometry
using LinearAlgebra
using Random

import ..BasicCompGeometry: dist, Point3F, AbsPntSeq

export convex_hull_3d, verify_convex_hull_3d

################################################################################
# Data Structures
################################################################################

"""
    Face

A triangular face on the 3D convex hull.
- `v`: Indices of the 3 vertices (counter-clockwise when viewed from outside).
- `normal`: Outward-facing unit normal vector.
- `offset`: Plane offset such that dot(normal, x) = offset for points on the face.
- `visible_points`: Indices of points that can "see" this face (in the conflict graph).
"""
mutable struct Face
    v::Vector{Int}
    normal::Point3F
    offset::Float64
    visible_points::Vector{Int}
    active::Bool

    function Face(v1::Int, v2::Int, v3::Int, pnts::Vector{Point3F})
        v = [v1, v2, v3]
        # Compute normal: (p2-p1) x (p3-p1)
        p1, p2, p3 = pnts[v1], pnts[v2], pnts[v3]
        n = cross(p2 - p1, p3 - p1)
        mag = norm(n)
        if mag < 1e-12
            # Degenerate face
            return new(v, point(0.0, 0.0, 0.0), 0.0, Int[], false)
        end
        n /= mag
        offset = dot(n, p1)
        return new(v, n, offset, Int[], true)
    end
end

"""
    is_visible(face, p)

Returns true if the point `p` is in the open half-space defined by the face 
and its outward normal (i.e., the point can "see" the face).
"""
@inline function is_visible(face::Face, p::Point3F, eps::Float64 = 1e-9)
    return dot(face.normal, p) > face.offset + eps
end

################################################################################
# Algorithm Helpers
################################################################################

"""
    find_initial_tetrahedron(pnts)

Finds 4 points that form a non-degenerate tetrahedron.
"""
function find_initial_tetrahedron(pnts::Vector{Point3F})
    n = length(pnts)
    n < 4 && return nothing

    # 1. Find two distinct points
    i1 = 1
    i2 = 0
    for i in 2:n
        if dist(pnts[i1], pnts[i]) > 1e-9
            i2 = i
            break
        end
    end
    i2 == 0 && return nothing

    # 2. Find a third point not collinear with the first two
    i3 = 0
    line_dir = pnts[i2] - pnts[i1]
    for i in 1:n
        (i == i1 || i == i2) && continue
        # Distance from point to line
        v = pnts[i] - pnts[i1]
        area = norm(cross(line_dir, v))
        if area > 1e-9
            i3 = i
            break
        end
    end
    i3 == 0 && return nothing

    # 3. Find a fourth point not coplanar with the first three
    i4 = 0
    p1, p2, p3 = pnts[i1], pnts[i2], pnts[i3]
    normal = cross(p2 - p1, p3 - p1)
    for i in 1:n
        (i == i1 || i == i2 || i == i3) && continue
        if abs(dot(normal, pnts[i] - p1)) > 1e-9
            i4 = i
            break
        end
    end
    i4 == 0 && return nothing

    return [i1, i2, i3, i4]
end

"""
    orient_tetrahedron!(indices, pnts)

Ensures the tetrahedron indices are oriented such that the first face (1,2,3) 
points away from the fourth point.
"""
function orient_tetrahedron!(idx, pnts)
    p1, p2, p3, p4 = pnts[idx[1]], pnts[idx[2]], pnts[idx[3]], pnts[idx[4]]
    n = cross(p2 - p1, p3 - p1)
    if dot(n, p4 - p1) > 0
        # Point 4 is in front of face (1,2,3), swap 2 and 3 to flip normal
        idx[2], idx[3] = idx[3], idx[2]
    end
end

################################################################################
# Main Algorithm
################################################################################

"""
    convex_hull_3d(points)

Computes the 3D convex hull of a set of points.
Returns a pair `(vertices, faces)` where `vertices` is a vector of indices 
into the original point set, and `faces` is a vector of triples of indices into 
the `vertices` vector.
"""
function convex_hull_3d(ps::AbsPntSeq{3,T}) where {T}
    return convex_hull_3d(Points(ps))
end

function convex_hull_3d(pnts_in::Vector{Point3F})
    n_in = length(pnts_in)
    if n_in < 4
        return collect(1:n_in), Vector{Int}[]
    end

    # 1. Initial tetrahedron
    tetra_idx = find_initial_tetrahedron(pnts_in)
    if tetra_idx === nothing
        # All points are coplanar (or worse)
        # TODO: Handle coplanar 2D hull in 3D
        return Int[], Vector{Int}[]
    end
    orient_tetrahedron!(tetra_idx, pnts_in)

    # 2. Setup internal structures
    faces = Vector{Face}()
    # Maps oriented edge (u, v) to face index
    edge_to_face = Dict{Tuple{Int, Int}, Int}()
    
    # Create 4 faces of the tetrahedron
    i1, i2, i3, i4 = tetra_idx
    push!(faces, Face(i1, i2, i3, pnts_in))
    push!(faces, Face(i1, i3, i4, pnts_in))
    push!(faces, Face(i1, i4, i2, pnts_in))
    push!(faces, Face(i2, i4, i3, pnts_in))

    for (i, f) in enumerate(faces)
        v = f.v
        edge_to_face[(v[1], v[2])] = i
        edge_to_face[(v[2], v[3])] = i
        edge_to_face[(v[3], v[1])] = i
    end

    # Remaining points
    claimed = fill(false, n_in)
    for i in tetra_idx; claimed[i] = true; end
    
    unclaimed = [i for i in 1:n_in if !claimed[i]]
    shuffle!(unclaimed)

    # Initial conflict graph
    for (f_idx, f) in enumerate(faces)
        for p_idx in unclaimed
            if is_visible(f, pnts_in[p_idx])
                push!(f.visible_points, p_idx)
            end
        end
    end

    # 3. Incremental add
    for p_idx in unclaimed
        # Find faces visible to this point
        visible_faces_indices = Int[]
        for (f_idx, f) in enumerate(faces)
            if f.active && p_idx in f.visible_points
                push!(visible_faces_indices, f_idx)
            end
        end

        if isempty(visible_faces_indices)
            continue
        end

        # Find the horizon edges
        edges_count = Dict{Tuple{Int, Int}, Int}()
        for f_idx in visible_faces_indices
            f = faces[f_idx]
            v = f.v
            e_list = [(v[1], v[2]), (v[2], v[3]), (v[3], v[1])]
            for e in e_list
                e_canon = e[1] < e[2] ? (e[1], e[2]) : (e[2], e[1])
                edges_count[e_canon] = get(edges_count, e_canon, 0) + 1
            end
        end

        # New faces to be added
        new_faces_info = [] # (nf_idx, u, v)

        for f_idx in visible_faces_indices
            f = faces[f_idx]
            f.active = false # Mark as inactive
            v = f.v
            e_list = [(v[1], v[2]), (v[2], v[3]), (v[3], v[1])]
            for e in e_list
                e_canon = e[1] < e[2] ? (e[1], e[2]) : (e[2], e[1])
                if edges_count[e_canon] == 1
                    # Horizon edge (u, v_e) from a visible face. 
                    # The non-visible face sharing this edge has it as (v_e, u).
                    u, v_e = e
                    f_nonvis_idx = edge_to_face[(v_e, u)]
                    f_nonvis = faces[f_nonvis_idx]
                    
                    # Create new face (u, v_e, p_idx)
                    new_f = Face(u, v_e, p_idx, pnts_in)
                    if new_f.active
                        push!(faces, new_f)
                        nf_idx = length(faces)
                        
                        # Update conflict list: candidates are points from 
                        # the two faces sharing the horizon edge.
                        candidates = [f.visible_points; f_nonvis.visible_points]
                        unique!(candidates)
                        for pt in candidates
                            if pt != p_idx && is_visible(new_f, pnts_in[pt])
                                push!(new_f.visible_points, pt)
                            end
                        end
                        push!(new_faces_info, (nf_idx, u, v_e))
                    end
                end
            end
        end

        # Update edge_to_face map
        # 1. Remove entries for deleted faces
        for f_idx in visible_faces_indices
            v = faces[f_idx].v
            delete!(edge_to_face, (v[1], v[2]))
            delete!(edge_to_face, (v[2], v[3]))
            delete!(edge_to_face, (v[3], v[1]))
        end
        # 2. Add entries for new faces
        for (nf_idx, u, v) in new_faces_info
            edge_to_face[(u, v)] = nf_idx
            edge_to_face[(v, p_idx)] = nf_idx
            edge_to_face[(p_idx, u)] = nf_idx
        end
    end

    # 4. Finalize result
    final_faces = Vector{Vector{Int}}()
    used_verts = Set{Int}()
    for f in faces
        if f.active
            push!(final_faces, f.v)
            for v in f.v
                push!(used_verts, v)
            end
        end
    end

    # Mapping from original indices to 1:m indices
    vert_list = sort(collect(used_verts))
    idx_map = Dict(orig => i for (i, orig) in enumerate(vert_list))
    
    mapped_faces = [ [idx_map[v] for v in f] for f in final_faces]

    return vert_list, mapped_faces
end

"""
    verify_convex_hull_3d(pnts, verts_idx, faces)

Verifies the geometric and topological properties of a 3D convex hull.
Checks:
1. Every face is a triangle.
2. Convexity: All points lie on one side of each face's plane.
3. Manifold: Every edge is shared by exactly two faces with opposite orientations.
4. Correctness: Every face has 3 edges.

Returns `(success::Bool, message::String)`.
"""
function verify_convex_hull_3d(pnts::Vector{Point3F}, verts_idx::Vector{Int}, faces::Vector{Vector{Int}})
    hull_pnts = pnts[verts_idx]
    
    # 1. Basic Triangle Check
    for (i, f) in enumerate(faces)
        if length(f) != 3
            return false, "Face $i is not a triangle (length: $(length(f)))"
        end
    end

    # 2. Convexity Check
    for (i, f) in enumerate(faces)
        p1, p2, p3 = hull_pnts[f[1]], hull_pnts[f[2]], hull_pnts[f[3]]
        # Normal pointing outward (CCW)
        n = cross(p2 - p1, p3 - p1)
        mag = norm(n)
        if mag < 1e-12
            return false, "Face $i is degenerate"
        end
        n /= mag
        offset = dot(n, p1)

        for (j, p) in enumerate(pnts)
            # Use a slightly relaxed epsilon for float precision
            if dot(n, p) > offset + 1e-7
                return false, "Convexity violated: Point $j is outside face $i by $(dot(n, p) - offset)"
            end
        end
    end

    # 3. Topological Check (Manifold & Orientation)
    # Map from oriented edge (u, v) to the face index that contains it
    edge_map = Dict{Tuple{Int, Int}, Int}()
    for (i, f) in enumerate(faces)
        edges = [(f[1], f[2]), (f[2], f[3]), (f[3], f[1])]
        for e in edges
            if haskey(edge_map, e)
                return false, "Edge $e is duplicated (shared by faces $(edge_map[e]) and $i with same orientation)"
            end
            edge_map[e] = i
        end
    end

    for (e, f_idx) in edge_map
        rev_e = (e[2], e[1])
        if !haskey(edge_map, rev_e)
            return false, "Edge $e (face $f_idx) has no matching reverse edge (hole in hull or inconsistent orientation)"
        end
    end

    return true, "Hull is valid"
end

end # module ConvexHull3D
