# Allocation-regression check of the hot entry points.
#
# Wall-clock benchmarks are too noisy for a gate, but the number of bytes a
# warmed-up call allocates is deterministic for a given Julia version. Each
# entry point below is measured (min of 3 runs of `@allocated` after warm-up,
# inside a function) on deterministic inputs, mostly at two sizes n and 4n so
# that a complexity blowup shows up as well, and compared against the
# committed baseline `test/alloc/baseline.toml`.
#
#   julia --project=. test/alloc/alloc_check.jl            # check
#   julia --project=. test/alloc/alloc_check.jl --update   # rewrite baseline
#
# A measurement fails when it exceeds baseline * 1.25 + 1024 bytes; one that
# dropped below baseline * 0.8 is reported (not failed) as a candidate for
# `--update`. The baseline records the Julia major.minor version; on another
# version the comparison is skipped (exit 0), since allocation counts change
# between Julia releases. Entry points that are missing or fail to run (e.g.
# when checking an older tree of the package) are reported as skipped.

using BasicCompGeometry, Random, Printf, TOML
using BasicCompGeometry: BBT, WSPD, MVBB

const B = BasicCompGeometry
const BASELINE = joinpath(@__DIR__, "baseline.toml")
const FAIL_RATIO = 1.25
const FAIL_SLACK = 1024
const DROP_RATIO = 0.8
const RUNS = 3

# Resolve `mod.sym` lazily, so a missing entry point is skipped, not a crash.
fn(mod::Module, sym::Symbol) =
    isdefined(mod, sym) ? getfield(mod, sym) : error("$mod.$sym is not defined")
fn(sub::Symbol, sym::Symbol) =
    isdefined(B, sub) ? fn(getfield(B, sub), sym) : error("$B.$sub is not defined")

# ---------------------------------------------------------------------------
# Deterministic inputs
# ---------------------------------------------------------------------------

rv2(n, s) = (r = Xoshiro(s); [point(rand(r), rand(r)) for _ in 1:n])
rp2(n, s) = PntSeq(rv2(n, s))
rp3(n, s) = (r = Xoshiro(s); PntSeq([point(rand(r), rand(r), rand(r)) for _ in 1:n]))

# Simple counterclockwise star-shaped polygon around the origin.
function star(n, s)
    r = Xoshiro(s)
    θ = sort(rand(r, n)) .* 2π
    return [(ρ = 0.5 + rand(r); point(ρ * cos(t), ρ * sin(t))) for t in θ]
end

# Counterclockwise x-monotone polygon: lower chain left to right, upper back.
function xmonotone(n, s)
    r = Xoshiro(s)
    k = n ÷ 2
    xs_lo = sort(rand(r, k))
    xs_hi = sort(rand(r, n - k); rev = true)
    lo = [point(-1.0 + 2x, -0.1 - rand(r)) for x in xs_lo]
    hi = [point(-1.0 + 2x, 0.1 + rand(r)) for x in xs_hi]
    return vcat([point(-1.5, 0.0)], lo, [point(1.5, 0.0)], hi)
end

# Regular convex k-gon of radius ρ.
cvx(k, ρ) = [point(ρ * cos(2π * i / k), ρ * sin(2π * i / k)) for i in 0:(k-1)]

# ---------------------------------------------------------------------------
# Entry points: name => (sizes, size -> zero-argument closure to measure)
# ---------------------------------------------------------------------------

# A loop of `m` calls of `f(canvas, args(i)...; kw...)` on a fresh canvas.
function ipe_loop(new::C, f::F, m, args::A, kw::NamedTuple = (;)) where {C,F,A}
    canvas = new()
    for i in 1:m
        f(canvas, args(i)...; kw...)
    end
    return canvas
end

function entries()
    E = Tuple{String,Vector{Int},Function}[]
    add!(make, name, sizes) = push!(E, (name, collect(Int, sizes), make))

    for alg in (:ear_cut, :sweep, :seidel)
        add!("triangulate_polygon :$alg", (1000, 4000)) do n
            f = fn(B, :triangulate_polygon)
            P = star(n, 1)
            () -> (Random.seed!(1); f(P; algorithm = alg))
        end
    end
    add!("triangulate_monotone", (1000, 4000)) do n
        f = fn(B, :triangulate_monotone)
        P = xmonotone(n, 2)
        () -> f(P)
    end
    add!("triangulate_ear_cut_u_kd_tree", (1000, 4000)) do n
        f = fn(B, :triangulate_ear_cut_u_kd_tree)
        P = star(n, 3)
        () -> f(P)
    end
    add!("convex_hull", (10_000, 40_000)) do n
        f = fn(B, :convex_hull)
        P = rp2(n, 4)
        () -> f(P)
    end
    add!("convex_hull_3d", (2000, 8000)) do n
        f = fn(:ConvexHull3D, :convex_hull_3d)
        P = rp3(n, 5)
        () -> (Random.seed!(1); f(P))
    end
    add!("BBT.Tree_init", (10_000, 40_000)) do n
        f = fn(BBT, :Tree_init)
        P = rp2(n, 6)
        () -> f(P)
    end
    add!("WSPD.init+expand!", (1000, 4000)) do n
        init = fn(WSPD, :init)
        expand! = fn(WSPD, :expand!)
        P = rp2(n, 7)
        () -> expand!(init(P, 2.0))
    end
    add!("approx_diameter", (10_000, 40_000)) do n
        f = fn(B, :approx_diameter)
        P = rp2(n, 8)
        () -> f(P, 0.05)
    end
    add!("MVBB.approx_mvbb", (5000, 20_000)) do n
        f = fn(MVBB, :approx_mvbb)
        P = rp3(n, 9)
        () -> (Random.seed!(1); f(P, 0.1))
    end
    add!("compute_largest_convex_subset", (25, 50)) do n
        f = fn(:LongestConvexSubset, :compute_largest_convex_subset)
        P = rp2(n, 10)
        () -> f(P)
    end
    add!("minkowski_sum", (500, 2000)) do n
        f = fn(B, :minkowski_sum)
        P = cvx(n, 1.0)
        Q = cvx(n ÷ 2 + 1, 0.3)
        () -> f(P, Q)
    end
    add!("minkowski_sum_decomposition", (100, 400)) do n
        f = fn(B, :minkowski_sum_decomposition)
        P = star(n, 12)
        Q = cvx(4, 0.05)
        () -> f(P, Q)
    end
    add!("visibility_polygon", (1000, 4000)) do n
        f = fn(B, :visibility_polygon)
        P = star(n, 13)
        () -> f(P, point(0.0, 0.0))
    end
    add!("polygon_kernel", (1000, 4000)) do n
        f = fn(B, :polygon_kernel)
        P = star(n, 14)
        () -> f(P)
    end
    add!("vertical_decomposition", (10, 40)) do n
        f = fn(B, :vertical_decomposition)
        r = Xoshiro(15)
        L = [Line(point(rand(r), 0.0), point(0.0, rand(r))) for _ in 1:n]
        view = BBox(point(0.0, 0.0), point(1.0, 1.0))
        () -> f(L, view)
    end
    add!("build_rmq", (10_000, 40_000)) do n
        f = fn(B, :build_rmq)
        k = rand(Xoshiro(16), n)
        () -> f(k)
    end
    add!("hierarchy_init+hierarchy_build!", (500, 2000)) do n
        init = fn(:PntSeqHierarchy, :hierarchy_init)
        build! = fn(:PntSeqHierarchy, :hierarchy_build!)
        r = Xoshiro(17)
        P = PntSeq(cumsum([point(randn(r), randn(r)) for _ in 1:n]))
        () -> build!(init(P))
    end
    add!("sample_uniformly", (10_000, 40_000)) do n
        f = fn(B, :sample_uniformly)
        P = rp2(n ÷ 4, 18)
        () -> f(P, n)
    end
    add!("union_bbox", (10_000, 40_000)) do n
        f = fn(B, :union_bbox)
        P = rp2(n + 1, 19)
        S = [Segment(P[i], P[i+1]) for i in 1:n]
        () -> f(S)
    end
    add!("dist_segment_segment x n", (10_000, 40_000)) do n
        f = fn(B, :dist_segment_segment)
        P = rv2(n + 3, 20)
        () -> begin
            s = 0.0
            for i in 1:n
                s += f(P[i], P[i+1], P[i+2], P[i+3])
            end
            s
        end
    end

    # IpeDraw primitives: a loop of m calls on a fresh canvas.
    m = 2000
    pts = rv2(m + 4, 21)
    add!("IpeDraw.draw_point! x$m", (m,)) do _
        f = fn(:IpeDraw, :draw_point!)
        new = fn(:IpeDraw, :IpeCanvas)
        () -> ipe_loop(new, f, m, i -> (pts[i],))
    end
    add!("IpeDraw.draw_segment! x$m", (m,)) do _
        f = fn(:IpeDraw, :draw_segment!)
        new = fn(:IpeDraw, :IpeCanvas)
        () -> ipe_loop(new, f, m, i -> (pts[i], pts[i+1]))
    end
    add!("IpeDraw.draw_circle! x$m", (m,)) do _
        f = fn(:IpeDraw, :draw_circle!)
        new = fn(:IpeDraw, :IpeCanvas)
        () -> ipe_loop(new, f, m, i -> (pts[i].x, pts[i].y, 0.1))
    end
    add!("IpeDraw.draw_box! x$m", (m,)) do _
        f = fn(:IpeDraw, :draw_box!)
        new = fn(:IpeDraw, :IpeCanvas)
        # Explicit attributes, so attribute formatting is not dwarfed by @sprintf.
        kw = (stroke = :red, fill = :blue, pen = :fat, dash = :dashed, opacity = :op50)
        () -> ipe_loop(new, f, m, i -> (BBox(pts[i], pts[i+1]),), kw)
    end
    add!("IpeDraw.draw_bezier! x$m", (m,)) do _
        f = fn(:IpeDraw, :draw_bezier!)
        new = fn(:IpeDraw, :IpeCanvas)
        C = fn(B, :CubicBezier)
        bz = [C(pts[i], pts[i+1], pts[i+2], pts[i+3]) for i in 1:m]
        kw = (
            stroke = :red,
            fill = :blue,
            pen = :fat,
            dash = :dashed,
            arrow = :normal,
            rarrow = :pointed,
        )
        () -> ipe_loop(new, f, m, i -> (bz[i],), kw)
    end
    add!("IpeDraw.mark! x$m", (m,)) do _
        f = fn(:IpeDraw, :mark!)
        new = fn(:IpeDraw, :IpeCanvas)
        () -> ipe_loop(new, f, m, i -> (pts[i],))
    end
    add!("IpeDraw.label! x$m", (m,)) do _
        f = fn(:IpeDraw, :label!)
        new = fn(:IpeDraw, :IpeCanvas)
        txt = ["p_{$i}" for i in 1:m]
        () -> ipe_loop(new, f, m, i -> (pts[i], txt[i]))
    end
    add!("IpeDraw.draw_polygon!", (2500, 10_000)) do n
        f = fn(:IpeDraw, :draw_polygon!)
        new = fn(:IpeDraw, :IpeCanvas)
        P = rv2(n, 22)
        () -> f(new(), P)
    end
    add!("IpeDraw.draw_spline!", (500, 2000)) do n
        f = fn(:IpeDraw, :draw_spline!)
        new = fn(:IpeDraw, :IpeCanvas)
        sp = fn(B, :interpolate_catmull_rom)(rp2(n, 23))
        () -> f(new(), sp)
    end
    add!("IpeDraw.draw_bspline!", (2500, 10_000)) do n
        f = fn(:IpeDraw, :draw_bspline!)
        new = fn(:IpeDraw, :IpeCanvas)
        P = rp2(n, 24)
        () -> f(new(), P)
    end
    return E
end

# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------

function measure(f::F) where {F}
    f()
    f()
    best = typemax(Int)
    for _ in 1:RUNS
        best = min(best, @allocated f())
    end
    return best
end

key(name, n) = "$name @ n=$n"

function measure_all()
    results = Dict{String,Int}()
    skipped = String[]
    for (name, sizes, make) in entries()
        for n in sizes
            k = key(name, n)
            try
                results[k] = measure(make(n))
            catch err
                msg = sprint(showerror, err; context = :limit => true)
                push!(skipped, "$k: " * first(split(msg, '\n')))
            end
        end
    end
    return results, skipped
end

julia_mm() = "$(VERSION.major).$(VERSION.minor)"

fmt_bytes(b) =
    b < 1024 ? "$(b) B" :
    b < 1024^2 ? @sprintf("%.1f KiB", b / 1024) : @sprintf("%.2f MiB", b / 1024^2)

function main(args)
    update = "--update" in args
    t0 = time()
    results, skipped = measure_all()
    elapsed = time() - t0

    if update
        if !isempty(skipped)
            println("alloc_check: cannot update the baseline; entry points failed:")
            foreach(s -> println("  ", s), skipped)
            exit(1)
        end
        data = Dict(
            "julia" => julia_mm(),
            "allocations" =>
                Dict(k => results[k] for k in sort(collect(keys(results)))),
        )
        open(BASELINE, "w") do io
            println(
                io,
                "# Allocation baseline for test/alloc/alloc_check.jl (bytes per call).",
            )
            println(
                io,
                "# Regenerate with: julia --project=. test/alloc/alloc_check.jl --update",
            )
            TOML.print(io, data; sorted = true)
        end
        @printf(
            "alloc_check: wrote %d measurements to %s (%.1f s)\n",
            length(results),
            BASELINE,
            elapsed
        )
        return
    end

    isfile(BASELINE) ||
        (println("alloc_check: no baseline at $BASELINE; run with --update"); exit(1))
    base = TOML.parsefile(BASELINE)
    if get(base, "julia", "") != julia_mm()
        println(
            "alloc_check: baseline was recorded on Julia $(get(base, "julia", "?")), " *
            "this is Julia $(julia_mm()); skipping the comparison " *
            "(regenerate with --update on this version).",
        )
        return
    end
    ref = base["allocations"]

    failed = Tuple{String,Int,Int}[]
    dropped = Tuple{String,Int,Int}[]
    w = maximum(length, keys(ref); init = 20) + 2
    @printf("%-*s %14s %14s %8s\n", w, "entry @ size", "baseline", "now", "ratio")
    for k in sort(collect(keys(ref)))
        b = Int(ref[k])
        if !haskey(results, k)
            @printf("%-*s %14s %14s %8s\n", w, k, fmt_bytes(b), "-", "skipped")
            continue
        end
        a = results[k]
        ratio = b == 0 ? (a == 0 ? 1.0 : Inf) : a / b
        status = ""
        if a > b * FAIL_RATIO + FAIL_SLACK
            push!(failed, (k, b, a))
            status = "  FAIL"
        elseif a < b * DROP_RATIO
            push!(dropped, (k, b, a))
            status = "  dropped"
        end
        @printf("%-*s %14s %14s %8.2f%s\n", w, k, fmt_bytes(b), fmt_bytes(a), ratio, status)
    end
    for k in sort(collect(keys(results)))
        haskey(ref, k) || println("note: $k is not in the baseline (run with --update)")
    end
    isempty(skipped) || println("\nskipped (missing or failing in this tree):")
    foreach(s -> println("  ", s), skipped)
    @printf("\nmeasured %d entries in %.1f s\n", length(results), elapsed)

    if !isempty(dropped)
        println(
            "\n$(length(dropped)) measurement(s) dropped below $(DROP_RATIO)x baseline; " *
            "consider tightening with --update:",
        )
        for (k, b, a) in dropped
            @printf("  %s: %s -> %s\n", k, fmt_bytes(b), fmt_bytes(a))
        end
    end
    if !isempty(failed)
        println(
            "\nalloc_check FAILED: $(length(failed)) measurement(s) exceed " *
            "baseline x $(FAIL_RATIO) + $(FAIL_SLACK) B:",
        )
        for (k, b, a) in failed
            @printf(
                "  %-*s %14s -> %14s  (%.2fx)\n",
                w,
                k,
                fmt_bytes(b),
                fmt_bytes(a),
                b == 0 ? Inf : a / b
            )
        end
        exit(1)
    end
    println(
        "alloc_check: all measurements within baseline x $(FAIL_RATIO) + $(FAIL_SLACK) B",
    )
end

main(ARGS)
