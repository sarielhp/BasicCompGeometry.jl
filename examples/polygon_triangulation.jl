#!/usr/bin/env julia

# Compare O(n log n) plane-sweep, expected O(n log* n) Seidel, :ear_cut_z_order, and O(n r)
# ear-clipping polygon triangulation, and draw a side-by-side partition + triangulation Ipe figure.
include(joinpath(@__DIR__, "support/runtime.jl"))

using QuickEnv
using BasicCompGeometry
using BasicCompGeometry.IpeDraw
using Printf
using Random

"""Build a CCW non-monotone horizontal comb polygon with `4m + 4` vertices and `2m` reflex notches."""
function make_comb_polygon(m_prongs::Int)
    pts = Point2F[point(0.0, 0.0)]
    for i in 1:m_prongs
        y_bot = 2.0 * (i - 1)
        x_notch = 1.2 + 1.4 * (i % 2) + 0.5 * sin(pi * i / (m_prongs + 1))
        push!(pts, point(10.0, y_bot))
        push!(pts, point(10.0, y_bot + 1.0))
        push!(pts, point(x_notch, y_bot + 1.2))
        push!(pts, point(x_notch, y_bot + 1.8))
    end
    push!(pts, point(10.0, 2.0 * m_prongs))
    push!(pts, point(10.0, 2.0 * m_prongs + 1.0))
    push!(pts, point(0.0, 2.0 * m_prongs + 1.0))
    return pts
end

"""Build a 22-vertex non-monotone fjord polygon in generic position for visual comparison."""
function make_fjord_polygon()
    return Point2F[
        point(0.5, 0.4),
        point(3.4, 0.2),
        point(2.5, 2.7),
        point(4.6, 1.5),
        point(5.5, 3.6),
        point(6.6, 0.6),
        point(9.6, 0.3),
        point(9.8, 3.2),
        point(7.5, 2.9),
        point(7.9, 5.1),
        point(9.7, 4.9),
        point(9.4, 7.6),
        point(7.0, 7.8),
        point(6.1, 5.4),
        point(4.9, 7.2),
        point(3.8, 4.7),
        point(2.9, 7.7),
        point(0.7, 7.4),
        point(0.3, 5.2),
        point(2.6, 5.8),
        point(1.4, 3.6),
        point(0.4, 2.3),
    ]
end

norm_diag(a::Int, b::Int) = a < b ? (a, b) : (b, a)

# Check Euler counts, triangle orientation and area conservation; unlike `@assert`, these
# checks are never compiled away.
function verify_triangulation(poly, res)
    n = length(poly)
    length(res.diagonals) == n - 3 || error("expected $(n - 3) diagonals, got $(length(res.diagonals))")
    length(res.triangles) == n - 2 || error("expected $(n - 2) triangles, got $(length(res.triangles))")
    poly_area = twice_area(poly)
    poly_area > 0 || error("polygon is not counterclockwise")
    tri_area = 0.0
    for (i, j, k) in res.triangles
        a = twice_area([poly[i], poly[j], poly[k]])
        a > 0 || error("triangle ($i, $j, $k) is not counterclockwise")
        tri_area += a
    end
    isapprox(tri_area, poly_area; rtol = 1e-8) || error("triangle areas do not sum to the polygon area")
    return nothing
end

function run_benchmarks()
    prongs_list = example_size([50, 200, 600, 1200], [8, 20, 40])
    # Warm up all six algorithms (with > 16 reflex vertices to exercise the spatial indices)
    warm = make_comb_polygon(12)
    for alg in (:sweep, :seidel, :ear_cut_z_order, :ear_cut_kd_tree, :ear_cut_u_kd_tree, :ear_cut)
        triangulate_polygon(warm; algorithm = alg, rng = MersenneTwister(1))
    end

    println("Polygon Triangulation Benchmark (:sweep, :seidel, :ear_cut_z_order, :ear_cut_kd_tree, :ear_cut_u_kd_tree, :ear_cut)")
    println("  Monotone?: false for all test combs (2m reflex vertices)")
    @printf(
        "%-7s | %-10s | %-10s | %-12s | %-11s | %-11s | %-11s | %-7s | %-7s\n",
        "n",
        ":sweep",
        ":seidel",
        ":z_order",
        ":kd_tree",
        ":u_kd_tree",
        ":ear_cut",
        "z/u_kd",
        "Area OK",
    )
    @printf(
        "%-7s | %10s | %10s | %12s | %11s | %11s | %11s | %7s | %-7s\n",
        "",
        "(ms)",
        "(ms)",
        "(ms)",
        "(ms)",
        "(ms)",
        "(ms)",
        "speedup",
        "",
    )
    println("-"^113)

    for m in prongs_list
        poly = make_comb_polygon(m)
        @assert !is_x_monotone(poly) && !is_y_monotone(poly)
        n = length(poly)

        t_sweep = @elapsed res_sweep = triangulate_polygon(poly)
        t_seidel = @elapsed res_seidel = triangulate_polygon(poly; algorithm = :seidel, rng = MersenneTwister(42))
        t_zorder = @elapsed res_zorder = triangulate_polygon(poly; algorithm = :ear_cut_z_order)
        t_kdtree = @elapsed res_kdtree = triangulate_polygon(poly; algorithm = :ear_cut_kd_tree)
        t_ukdtree = @elapsed res_ukdtree = triangulate_polygon(poly; algorithm = :ear_cut_u_kd_tree)
        t_ear = @elapsed res_ear = triangulate_polygon(poly; algorithm = :ear_cut)

        verify_triangulation(poly, res_sweep)
        verify_triangulation(poly, res_seidel)
        verify_triangulation(poly, res_zorder)
        verify_triangulation(poly, res_kdtree)
        verify_triangulation(poly, res_ukdtree)
        verify_triangulation(poly, res_ear)

        speedup = t_zorder / max(t_ukdtree, 1e-9)
        @printf(
            "%-7d | %10.3f | %10.3f | %12.3f | %11.3f | %11.3f | %11.3f | %6.2fx | %-7s\n",
            n,
            1e3 * t_sweep,
            1e3 * t_seidel,
            1e3 * t_zorder,
            1e3 * t_kdtree,
            1e3 * t_ukdtree,
            1e3 * t_ear,
            speedup,
            "pass",
        )
    end
    println("Verified: every algorithm passed the Euler-count, orientation, and area checks")
end

function draw_panel!(fig, poly, part_diags, all_diags, theme, shift_x::Float64, caption_tex::String)
    shifted = [point(p.x + shift_x, p.y) for p in poly]
    part_set = Set(norm_diag(a, b) for (a, b) in part_diags)
    all_set = Set(norm_diag(a, b) for (a, b) in all_diags)
    rest_diags = sort!(collect(setdiff(all_set, part_set)))

    draw!(fig, shifted; style = theme.polygon)
    for (a, b) in rest_diags
        draw!(fig, Segment(shifted[a], shifted[b]); style = theme.secondary_diag)
    end
    for (a, b) in sort!(collect(part_set))
        draw!(fig, Segment(shifted[a], shifted[b]); style = theme.partition_diag)
    end
    mark!(fig, shifted; stroke = :navy, fill = :navy, shape = :disk, size = :small)
    label!(fig, point(shift_x + 5.0, -0.55), caption_tex; anchor = :center, style = :normal, size = :small)
end

function build_figure()
    fig_poly = make_fjord_polygon()
    @assert !is_x_monotone(fig_poly) && !is_y_monotone(fig_poly)

    # Default :sweep triangulation + its monotone-partition diagonals
    res_sweep = triangulate_polygon(fig_poly)
    sweep_part = BasicCompGeometry._sweep_partition_diagonals(FloatingKernel(), fig_poly)
    verify_triangulation(fig_poly, res_sweep)

    # Seidel triangulation + its trapezoidal-partition diagonals with shared RNG seed
    rng_seidel = MersenneTwister(17)
    trap_map = BasicCompGeometry._seidel_build_trapezoidal_map(FloatingKernel(), fig_poly, copy(rng_seidel))
    seidel_part = BasicCompGeometry._seidel_extract_diagonals(trap_map)
    res_seidel = triangulate_seidel(fig_poly; rng = rng_seidel)
    verify_triangulation(fig_poly, res_seidel)

    output_dir = example_output()
    mkpath(output_dir)
    pdf_path = joinpath(output_dir, "polygon_triangulation.pdf")

    dx = 11.4
    fit_box = BBox(point(-0.2, -1.2), point(dx + 10.2, 9.2))
    theme = Theme(
        polygon = Style(fill = :lightblue, fill_opacity = 0.2, stroke = :navy, pen = :heavier),
        partition_diag = Style(stroke = :red, pen = :heavier, dash = :dashed),
        secondary_diag = Style(stroke = :seagreen, pen = :thin),
    )

    figure(pdf_path; fit = fit_box, margin = 22, preview = "--preview" in ARGS) do fig
        draw_panel!(
            fig,
            fig_poly,
            sweep_part,
            res_sweep.diagonals,
            theme,
            0.0,
            raw"Plane-Sweep ($O(n \log n)$): monotone partition + stack",
        )
        draw_panel!(
            fig,
            fig_poly,
            seidel_part,
            res_seidel.diagonals,
            theme,
            dx,
            raw"Seidel ($O(n \log^* n)$ exp.): trapezoidal split + mountain",
        )
        legend!(
            fig,
            [
                raw"Simple Polygon $\partial P$" => theme.polygon,
                raw"Partition Diagonals" => theme.partition_diag,
                raw"Completion Diagonals" => theme.secondary_diag,
            ];
            position = :northwest,
            width = 165,
        )
    end

    println("Generated: ", pdf_path)
end

run_benchmarks()
build_figure()
