# Runtime-dispatch check of the hot entry points with JET.report_opt.
#
# Run by `tools/gate --full` (not by Pkg.test, so CI and users never install
# JET). Any dynamic dispatch in BasicCompGeometry code reachable from these
# calls fails the check: an abstract field type, an abstract container element
# type, or a type-unstable local makes a hot loop dispatch at run time.
#
#   julia --project=test/jet -e 'using Pkg; Pkg.instantiate()'
#   julia --project=test/jet test/jet/jet_opt.jl

using BasicCompGeometry, JET, Random
using BasicCompGeometry: BBT, WSPD, MVBB
using BasicCompGeometry.ConvexHull3D,
    BasicCompGeometry.LongestConvexSubset, BasicCompGeometry.PntSeqHierarchy

const B = BasicCompGeometry

# Error-message construction (string interpolation of geometry before a throw)
# is not a hot path; do not descend into it.
const _COLD = (Base.string, Base.print, Base.show, Base.print_to_string, Base.sprint)
_hot(@nospecialize(f)) = !any(c -> f === typeof(c), _COLD)
const _COLD_NAMES = (:string, :print, :show, :print_to_string, :sprint)
_cold_frame(fr) = (m = fr.linfo.def; m isa Method && m.name in _COLD_NAMES)
_hot_reports(result) = filter(r -> !any(_cold_frame, r.vst), JET.get_reports(result))

function star(n, rng)
    θ = sort(rand(rng, n)) .* 2π
    return [(r = 0.5 + rand(rng); point(r * cos(t), r * sin(t))) for t in θ]
end

function main()
    rng = Xoshiro(1)
    P2 = PntSeq([point(rand(rng), rand(rng)) for _ in 1:200])
    P3 = [point(rand(rng), rand(rng), rand(rng)) for _ in 1:200]
    poly = star(200, rng)
    sq = [point(0.0, 0.0), point(1.0, 0.0), point(1.0, 1.0), point(0.0, 1.0)]
    lines = [Line(point(rand(rng), 0.0), point(0.0, rand(rng))) for _ in 1:10]
    view = BBox(point(0.0, 0.0), point(1.0, 1.0))
    W = WSPD.init(P2, 2.0)
    keys = rand(rng, 1000)
    h = hierarchy_init(P2)

    cases = [
        (
            "convex_hull",
            () -> @report_opt target_modules = (B,) function_filter = _hot convex_hull(P2)
        ),
        (
            "convex_hull_3d",
            () ->
                @report_opt target_modules = (B,) function_filter = _hot convex_hull_3d(P3)
        ),
        (
            "BBT.Tree_init",
            () ->
                @report_opt target_modules = (B,) function_filter = _hot BBT.Tree_init(P2)
        ),
        (
            "WSPD.init",
            () -> @report_opt target_modules = (B,) function_filter = _hot WSPD.init(
                P2,
                2.0,
            )
        ),
        (
            "WSPD.expand!",
            () ->
                @report_opt target_modules = (B,) function_filter = _hot WSPD.expand!(W)
        ),
        (
            "approx_diameter",
            () ->
                @report_opt target_modules = (B,) function_filter = _hot approx_diameter(
                    P2,
                    0.1,
                )
        ),
        (
            "MVBB.approx_mvbb",
            () ->
                @report_opt target_modules = (B,) function_filter = _hot MVBB.approx_mvbb(
                    PntSeq(P3),
                    0.1,
                )
        ),
        (
            "compute_largest_convex_subset",
            () ->
                @report_opt target_modules = (B,) function_filter = _hot compute_largest_convex_subset(
                    P2,
                )
        ),
        (
            "triangulate :ear_cut",
            () ->
                @report_opt target_modules = (B,) function_filter = _hot triangulate_polygon(
                    poly;
                    algorithm = :ear_cut,
                )
        ),
        (
            "triangulate :sweep",
            () ->
                @report_opt target_modules = (B,) function_filter = _hot triangulate_polygon(
                    poly;
                    algorithm = :sweep,
                )
        ),
        (
            "triangulate :seidel",
            () ->
                @report_opt target_modules = (B,) function_filter = _hot triangulate_polygon(
                    poly;
                    algorithm = :seidel,
                )
        ),
        (
            "triangulate_monotone",
            () ->
                @report_opt target_modules = (B,) function_filter = _hot triangulate_monotone(
                    poly,
                )
        ),
        (
            "triangulate_ear_cut small",
            () ->
                @report_opt target_modules = (B,) function_filter = _hot triangulate_ear_cut(
                    poly[1:4:end],
                )
        ),
        (
            "triangulate_ear_cut_u_kd_tree",
            () ->
                @report_opt target_modules = (B,) function_filter = _hot triangulate_ear_cut_u_kd_tree(
                    poly,
                )
        ),
        (
            "minkowski_sum",
            () ->
                @report_opt target_modules = (B,) function_filter = _hot minkowski_sum(
                    sq,
                    sq,
                )
        ),
        (
            "minkowski_sum_decomposition",
            () ->
                @report_opt target_modules = (B,) function_filter = _hot minkowski_sum_decomposition(
                    poly,
                    sq,
                )
        ),
        (
            "visibility_polygon",
            () ->
                @report_opt target_modules = (B,) function_filter = _hot visibility_polygon(
                    poly,
                    point(0.0, 0.0),
                )
        ),
        (
            "polygon_kernel",
            () ->
                @report_opt target_modules = (B,) function_filter = _hot polygon_kernel(
                    poly,
                )
        ),
        (
            "vertical_decomposition",
            () ->
                @report_opt target_modules = (B,) function_filter = _hot vertical_decomposition(
                    lines,
                    view,
                )
        ),
        (
            "build_rmq",
            () -> @report_opt target_modules = (B,) function_filter = _hot build_rmq(keys)
        ),
        (
            "hierarchy_build!",
            () ->
                @report_opt target_modules = (B,) function_filter = _hot hierarchy_build!(
                    h,
                )
        ),
        (
            "sample_uniformly",
            () ->
                @report_opt target_modules = (B,) function_filter = _hot sample_uniformly(
                    P2,
                    500,
                )
        ),
    ]

    failed = 0
    for (name, run) in cases
        result = run()
        reports = _hot_reports(result)
        n = length(reports)
        println(rpad(name, 32), n == 0 ? "ok" : "$n runtime-dispatch report(s)")
        if n > 0
            failed += 1
            foreach(r -> (JET.print_report(stdout, r); println()), reports)
        end
    end
    failed == 0 || (
        println(
            "JET report_opt: $failed of $(length(cases)) entry points dispatch at run time",
        );
        exit(1)
    )
    println("JET report_opt: all $(length(cases)) entry points free of runtime dispatch")
end

main()
