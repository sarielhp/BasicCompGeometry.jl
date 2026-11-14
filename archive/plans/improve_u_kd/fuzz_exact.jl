include("lib.jl")
function randpoly(rng, n, G)
    for _ = 1:500
        pts = unique([
            point(Rational{BigInt}(rand(rng, 0:G)), Rational{BigInt}(rand(rng, 0:G)))
            for _ = 1:n
        ])
        length(pts)<3 && continue
        P=copy(pts)
        shuffle!(rng, P)
        m=length(P)
        for it = 1:2000
            changed=false
            for i = 1:m, j = (i+2):m
                (i==1 && j==m) && continue
                if seg_intersect(P[i], P[i+1], P[j], P[mod1(j+1, m)])
                    reverse!(P, i+1, j)
                    changed=true
                end
            end
            changed || break
        end
        is_simple(P) || continue
        area2(P)<0 && reverse!(P)
        area2(P)==0 && continue
        return P
    end
    nothing
end
rng=Xoshiro(parse(Int, ARGS[1]));
fails=Dict{Symbol,Int}();
tried=0
algs =
    length(ARGS) >= 3 ? Tuple(Symbol.(split(ARGS[3], ","))) :
    (:sweep, :seidel, :ear_cut, :ear_cut_z_order, :ear_cut_kd_tree, :ear_cut_u_kd_tree)
for trial = 1:parse(Int, ARGS[2])
    n=rand(rng, 4:60)
    G=rand(rng, [4, 8, 15, 30])
    P=randpoly(rng, n, G)
    P===nothing && continue
    global tried+=1
    for a in algs
        r = try
            verify(
                P,
                triangulate_polygon(
                    P;
                    kernel = ExactKernel(),
                    algorithm = a,
                    rng = Xoshiro(trial),
                ),
            )
        catch e
            "EXC: "*sprint(showerror, e)[1:min(end, 100)]
        end
        if r!==nothing
            fails[a]=get(fails, a, 0)+1
            fails[a]<=2 &&
                println("FAIL $a $r P=$([(Float64(p.x),Float64(p.y)) for p in P])")
        end
    end
end
println("tried=$tried fails=$fails")
