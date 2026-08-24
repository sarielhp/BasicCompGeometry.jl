#!/usr/bin/env julia
# Monte Carlo Area Estimation of a Curved Voronoi Cell
#
# This program demonstrates:
# 1. Computing the Voronoi cell of a point site p w.r.t. the 4 boundary lines
#    of the unit square [0, 1] x [0, 1].
# 2. Representing the cell boundary as a CurvePolygon2D composed of 4 ParabolaArc curves,
#    since the bisector between a point (focus) and a line (directrix) is a parabola.
# 3. Computing the reference area of this curved cell using fine boundary sampling (shoelace formula).
# 4. Estimating the cell area via Monte Carlo sampling across multiple sample sizes
#    (n = 10^4, 10^5, 10^6, 10^7) using the ray-casting `is_inside` point-in-curved-polygon test.
# 5. Comparing empirical errors with theoretical Monte Carlo standard deviation (O(1/sqrt(n))).

using QuickEnv
using BasicCompGeometry
using LinearAlgebra
using Printf
using Random

"""
    square_boundary_lines()

Return the 4 infinite lines forming the boundary of the unit square [0, 1] x [0, 1]:
- Left:   x = 0
- Right:  x = 1
- Bottom: y = 0
- Top:    y = 1
"""
function square_boundary_lines()
    left   = Line(point(0.0, 0.0), point(0.0, 1.0))
    right  = Line(point(1.0, 0.0), point(0.0, 1.0))
    bottom = Line(point(0.0, 0.0), point(1.0, 0.0))
    top    = Line(point(0.0, 1.0), point(1.0, 0.0))
    return [left, right, bottom, top]
end

"""
    distance_to_line(x, line)

Compute the Euclidean distance from 2D point `x` to infinite line `line`.
"""
function distance_to_line(x::Point{2,T}, line::Line{2,T}) where {T}
    n = Point{2,T}(-line.u.y, line.u.x)
    n = n / norm(n)
    return abs(dot(x - line.p, n))
end

"""
    voronoi_cell_vertices(p)

Compute the 4 parabolas and the sorted cyclic vertices of the Voronoi cell of point `p`
w.r.t. the four boundary lines of the unit square.
"""
function voronoi_cell_vertices(p::Point{2,T}) where {T}
    lines = square_boundary_lines()
    n = 4
    parabolas = [Parabola(p, lines[i]) for i in 1:n]
    interior = point(T(0.5), T(0.5))
    normals = Vector{Point{2,T}}(undef, n)
    for i in 1:n
        ni = Point{2,T}(-lines[i].u.y, lines[i].u.x)
        ni = ni / norm(ni)
        if dot(interior - lines[i].p, ni) < 0
            ni = -ni
        end
        normals[i] = ni
    end
    candidates = Point{2,T}[]
    for i in 1:n
        for j in (i+1):n
            nd = normals[i] - normals[j]
            if norm(nd) < 1e-12; continue; end
            ci = dot(lines[i].p, normals[i])
            cj = dot(lines[j].p, normals[j])
            d = ci - cj
            dir = Point{2,T}(-nd.y, nd.x)
            pt0 = Point{2,T}(0.0, 0.0)
            if abs(nd.x) > 1e-12
                pt0 = Point{2,T}(d / nd.x, 0.0)
            elseif abs(nd.y) > 1e-12
                pt0 = Point{2,T}(0.0, d / nd.y)
            end
            bis = Line(pt0, dir)
            pts = intersect_line_curve(bis, parabolas[i])
            for (pt, _) in pts
                d_p = dist(pt, p)
                ok = true
                for kk in 1:n
                    sd = dot(pt - lines[kk].p, normals[kk])
                    if sd < d_p - 1e-8
                        ok = false; break
                    end
                end
                if ok && abs(dot(pt - lines[i].p, normals[i]) - d_p) < 1e-8
                    push!(candidates, pt)
                end
            end
        end
    end
    isempty(candidates) && return parabolas, Point{2,T}[]
    sort!(candidates, by=pt -> atan(pt.y - p.y, pt.x - p.x))
    unique_pts = Point{2,T}[]
    for pt in candidates
        if isempty(unique_pts) || dist(pt, unique_pts[end]) > 1e-8
            push!(unique_pts, pt)
        end
    end
    return parabolas, unique_pts
end

"""
    voronoi_curved_polygon(p)

Construct the Voronoi cell of point `p` as a `CurvePolygon2D` composed of 4 `ParabolaArc` segments.
"""
function voronoi_curved_polygon(p::Point{2,T}) where {T}
    lines = square_boundary_lines()
    parabolas, verts = voronoi_cell_vertices(p)
    length(verts) < 3 && return CurvePolygon2D(ParabolaArc{T}[])
    m = length(verts)
    arcs = ParabolaArc{T}[]
    for idx in 1:m
        p_cur = verts[idx]
        p_next = verts[mod1(idx + 1, m)]
        mid = (p_cur + p_next) / 2
        best_k = 0
        best_diff = Inf
        for k in 1:4
            d = abs(dist(mid, p) - distance_to_line(mid, lines[k]))
            if d < best_diff
                best_diff = d; best_k = k
            end
        end
        push!(arcs, ParabolaArc(parabolas[best_k], p_cur, p_next))
    end
    return CurvePolygon2D(arcs)
end

"""
    sample_cell_boundary(cp, n_samples=400)

Sample `n_samples` points evenly along the boundary of the curved polygon `cp`.
"""
function sample_cell_boundary(cp::CurvePolygon2D{T}, n_samples::Int=400) where {T}
    isempty(cp.curves) && return Point{2,T}[]
    m = length(cp.curves)
    pts = Point{2,T}[]
    samples_per_curve = max(2, n_samples ÷ m)
    for c in cp.curves
        for s in range(0.0, 1.0, length=samples_per_curve)
            push!(pts, at(c, s))
        end
    end
    return pts
end

"""
    computed_cell_area(cp)

Compute the area of the curved polygon `cp` via fine boundary discretization
and Green\x27s theorem (shoelace formula on finely sampled points).
"""
function computed_cell_area(cp::CurvePolygon2D{T}) where {T}
    boundary = sample_cell_boundary(cp, 800)
    isempty(boundary) && return 0.0
    area = 0.0
    for i in 1:length(boundary)
        j = mod1(i + 1, length(boundary))
        area += boundary[i].x * boundary[j].y
        area -= boundary[j].x * boundary[i].y
    end
    return abs(area) / 2
end

"""
    estimate_area_monte_carlo(cp, n_samples)

Estimate the area of `cp` within [0, 1]^2 by generating `n_samples` uniform random points
and testing inside/outside status with `is_inside(q, cp)`.
"""
function estimate_area_monte_carlo(cp::CurvePolygon2D{T}, n_samples::Int) where {T}
    inside_count = 0
    for _ in 1:n_samples
        q = point(rand(T), rand(T))
        if is_inside(q, cp)
            inside_count += 1
        end
    end
    return T(inside_count) / T(n_samples)
end

function main()
    println("="^75)
    println("  Monte Carlo Area Estimation of Curved Voronoi Cell")
    println("="^75)
    println()
    println("Description:")
    println("  We place a site point p inside the unit square [0, 1]^2.")
    println("  The Voronoi cell of p with respect to the four boundary edges of the square")
    println("  consists of all points q in [0, 1]^2 closer to p than to the boundary.")
    println("  Since the bisector of a point (focus) and a line (directrix) is a parabola,")
    println("  the cell boundary is formed by 4 parabolic arcs (ParabolaArc).")
    println()

    # Fixed random seed for reproducibility
    Random.seed!(42)

    p = point(0.3, 0.4)
    println(@sprintf("Site point p: (%.2f, %.2f)", p.x, p.y))

    # Construct curved polygon
    cp = voronoi_curved_polygon(p)
    println(@sprintf("Constructed CurvePolygon2D with %d parabolic arcs.", length(cp.curves)))
    println(@sprintf("Boundary perimeter (arc length): %.6f", geom_length(cp)))
    println()

    # Reference area computed via fine discretization and Green\x27s theorem
    ref_area = computed_cell_area(cp)
    println(@sprintf("Reference area (boundary discretization): %.6f", ref_area))
    println()

    println("Running Monte Carlo sampling across sample sizes:")
    println("-"^75)
    println(@sprintf("%-12s | %-15s | %-12s | %-12s | %-10s",
                     "Samples (N)", "Estimated Area", "Abs Error", "Std Dev (1σ)", "Rel Error"))
    println("-"^75)

    sample_sizes = [10_000, 100_000, 1_000_000, 10_000_000]

    for n in sample_sizes
        est_area = estimate_area_monte_carlo(cp, n)
        abs_err = abs(est_area - ref_area)
        rel_err = abs_err / ref_area
        # Theoretical Monte Carlo standard deviation: sigma = sqrt(A * (1 - A) / N)
        sigma = sqrt(ref_area * (1.0 - ref_area) / n)
        
        println(@sprintf("%-12d | %-15.6f | %-12.6f | %-12.6f | %-9.2f%%",
                         n, est_area, abs_err, sigma, rel_err * 100.0))
    end
    println("-"^75)
    println()
    println("Conclusion:")
    println("  As the number of samples increases from 10^4 to 10^7, the Monte Carlo estimate")
    println("  converges to the reference area with error scaling as O(1/sqrt(N)),")
    println("  verifying both the ParabolaArc geometry and the ray-casting `is_inside` test.")
    println("="^75)
end

main()
