# Counters and timings for :diagonal_split vs :sweep and :ear_cut_u_kd_tree.
# Usage: julia --project=../.. bench.jl [log2 sizes...]   (default: 16 18)
include(joinpath(@__DIR__, "..", "improve_u_kd", "polys.jl"))
const B = BasicCompGeometry

function hcomb(m)
    pts = [point(0.0, 0.0)]
    for i in 1:m
        yb = 2.0 * (i - 1); xn = 1.2 + 1.4 * (i % 2) + 0.5 * sin(pi * i / (m + 1))
        push!(pts, point(10.0, yb), point(10.0, yb + 1.0), point(xn, yb + 1.2), point(xn, yb + 1.8))
    end
    push!(pts, point(10.0, 2.0 * m), point(10.0, 2.0 * m + 1.0), point(0.0, 2.0 * m + 1.0)); pts
end

exps = isempty(ARGS) ? [16, 18] : parse.(Int, ARGS)
for e in exps
    n = 2^e
    for (nm, P) in (("star", star(n)), ("altstar", altstar(n)), ("coast", coast(n)),
                    ("spiral", spiral(n)), ("comb", hcomb(n ÷ 4)))
        nP = length(P)
        _, s = B._triangulate_diagonal_split_impl(P, FloatingKernel(), Xoshiro(1))
        println(rpad("2^$e $nm", 14), " tests/n=", rpad(round(s[1] / nP, digits=3), 6),
                " fail/n=", rpad(round(s[2] / nP, digits=3), 6),
                " nodes/test=", rpad(round(s[3] / max(s[1], 1), digits=1), 6),
                " edges/test=", rpad(round(s[4] / max(s[1], 1), digits=1), 6),
                " splits/n=", rpad(round(s[5] / nP, digits=3), 6),
                " gallop/n=", rpad(round(s[8] / nP, digits=3), 6), " fans=", s[7], " fallbacks=", s[6])
        print(" "^14)
        for a in (:sweep, :ear_cut_u_kd_tree, :diagonal_split)
            triangulate_polygon(P; algorithm=a, rng=Xoshiro(1))
            t = minimum(@elapsed(triangulate_polygon(P; algorithm=a, rng=Xoshiro(1))) for _ in 1:3)
            print(rpad(" $a=$(round(1e3t, digits=1))ms", 30))
        end
        println()
    end
end
