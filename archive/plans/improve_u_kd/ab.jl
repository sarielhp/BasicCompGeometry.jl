include("polys.jl")
function hcomb(m)
    pts = [point(0.0, 0.0)]
    for i in 1:m
        yb = 2.0 * (i - 1)
        xn = 1.2 + 1.4 * (i % 2) + 0.5 * sin(pi * i / (m + 1))
        push!(
            pts,
            point(10.0, yb),
            point(10.0, yb + 1.0),
            point(xn, yb + 1.2),
            point(xn, yb + 1.8),
        )
    end
    push!(pts, point(10.0, 2.0 * m), point(10.0, 2.0 * m + 1.0), point(0.0, 2.0 * m + 1.0))
    pts
end
Base.include(BasicCompGeometry, joinpath(@__DIR__, "old_run.jl"))
const B = BasicCompGeometry
function run_with(P, idxkind, old::Bool)
    n = length(P)
    k = FloatingKernel()
    prev_v, next_v, is_active = B._init_ear_clip_ring(n)
    is_reflex = falses(n)
    rc = 0
    for i in 1:n
        if orient_sign(k, P[prev_v[i]], P[i], P[next_v[i]]) <= 0
            is_reflex[i] = true
            rc += 1
        end
    end
    idx = if idxkind == :u
        B._build_kdtree_reflex_index(P, Vector{Bool}(is_reflex), rc, Val(true))
    elseif idxkind == :kd
        B._build_kdtree_reflex_index(P, Vector{Bool}(is_reflex), rc, Val(false))
    else
        ok, a, b, c, d = B._compute_z_bounds(P)
        B._build_zorder_reflex_index(P, is_reflex, rc, a, b, c, d)
    end
    st = B._EarClipState(k, P, prev_v, next_v, is_active, idx)
    return old ? B._ear_clip_run_old!(st) : B._ear_clip_run!(st)
end
cases = [
    ("star", star(2^17)),
    ("altstar", altstar(2^18)),
    ("coast", coast(2^18)),
    ("spiral", spiral(2^18)),
    ("comb", hcomb(2^12)),
]
for kind in (:u, :kd, :z)
    for (nm, P) in cases
        (kind == :z && nm == "star") && continue
        run_with(P, kind, true)
        run_with(P, kind, false)
        to = Inf
        tn = Inf
        for _ in 1:4
            to = min(to, @elapsed run_with(P, kind, true))
            tn = min(tn, @elapsed run_with(P, kind, false))
        end
        println(
            rpad("$kind $nm", 14),
            " old=",
            rpad(round(1e3to, digits = 1), 9),
            " new=",
            rpad(round(1e3tn, digits = 1), 9),
            " ratio=",
            round(tn/to, digits = 2),
        )
    end
end
