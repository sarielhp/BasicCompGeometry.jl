include("polys.jl")
const B = BasicCompGeometry
cases = [
    ("star", star(2^17)),
    ("altstar", altstar(2^18)),
    ("coast", coast(2^18)),
    ("spiral", spiral(2^18)),
]
res = Dict{String,Vector{Float64}}()
for (nm, P) in cases
    triangulate_ear_cut_u_kd_tree(P)
    res[nm] = [minimum(@elapsed(triangulate_ear_cut_u_kd_tree(P)) for _ in 1:4)]
end
@eval B @inline _leaf_may_block(
    q::_FlatTriQueryF64,
    px::Vector{Float64},
    py::Vector{Float64},
    rs::Int,
    le::Int,
)::Bool = true
for (nm, P) in cases
    triangulate_ear_cut_u_kd_tree(P)
    push!(res[nm], minimum(@elapsed(triangulate_ear_cut_u_kd_tree(P)) for _ in 1:4))
end
for (nm, _) in cases
    a, b = res[nm]
    println(
        rpad(nm, 10),
        " with prescan=",
        round(1e3a, digits = 1),
        " without=",
        round(1e3b, digits = 1),
        " ratio=",
        round(a/b, digits = 2),
    )
end
