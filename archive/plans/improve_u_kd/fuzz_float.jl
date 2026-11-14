include("lib.jl")
function randpoly(rng, n, G, scale)
    for _ = 1:500
        pts = unique([
            point(
                Rational{BigInt}(Float64(rand(rng, 0:G)*scale)),
                Rational{BigInt}(Float64(rand(rng, 0:G)*scale)),
            ) for _ = 1:n
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
fails=0;
tried=0;
argerr=0
log = open(
    joinpath(
        mkpath(joinpath(@__DIR__, "..", "..", "..", "output", "tmp")),
        "fuzz_float_last_$(ARGS[1]).txt",
    ),
    "w",
)
for trial = 1:parse(Int, ARGS[2])
    n=rand(rng, 4:25)
    G=rand(rng, [4, 6, 10])
    sc=rand(rng, [0.1, 0.3, 1/3, 0.7, 1e-9*0.1])
    P=randpoly(rng, n, G, sc)
    P===nothing && continue
    Pf=[point(Float64(p.x), Float64(p.y)) for p in P]
    for s = 1:3
        global tried+=1
        seekstart(log)
        println(log, "seed=$s Pf=", Pf)
        flush(log)
        r = try
            verify(P, triangulate_polygon(Pf; algorithm = Symbol(ARGS[3]), rng = Xoshiro(s)))
        catch e
            e isa ArgumentError ? (global argerr+=1; nothing) :
            "EXC: "*sprint(showerror, e)[1:min(end, 100)]
        end
        if r!==nothing
            global fails+=1
            fails<=4 && println("FAIL seed=$s $r\n  Pf=$Pf")
            break
        end
    end
end
println(ARGS[3], " tried=$tried argerr=$argerr fails(wrong or non-Arg exc)=$fails")
