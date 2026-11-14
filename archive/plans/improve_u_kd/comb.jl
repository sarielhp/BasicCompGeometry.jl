using BasicCompGeometry
const C = zeros(Int, 4)
@eval BasicCompGeometry begin
    @inline function _is_ear_vertex(st::_EarClipState, i::Int)::Bool
        Main.C[1] += 1
        _is_reflex_vertex(st.rindex, i) && return false
        Main.C[2] += 1
        pr = st.prev_v[i]
        nx = st.next_v[i]
        poly = st.poly
        r =
            !_ear_blocked(
                st.rindex,
                poly,
                st.kernel,
                pr,
                i,
                nx,
                poly[pr],
                poly[i],
                poly[nx],
            )
        r && (Main.C[3] += 1)
        return r
    end
    function _find_fallback_ear(st::_EarClipState, hint::Int)::Int
        Main.C[4] += 1
        start_v = st.is_active[hint] ? hint : findfirst(st.is_active)
        start_v === nothing && return 0
        curr = start_v
        while true
            _is_ear_vertex(st, curr) && return curr
            curr = st.next_v[curr]
            curr == start_v && return 0
        end
    end
end
function _hcomb(m_prongs::Int)
    pts = [point(0.0, 0.0)]
    for i in 1:m_prongs
        yb = 2.0 * (i - 1)
        xn = 1.2 + 1.4 * (i % 2) + 0.5 * sin(pi * i / (m_prongs + 1))
        push!(
            pts,
            point(10.0, yb),
            point(10.0, yb + 1.0),
            point(xn, yb + 1.2),
            point(xn, yb + 1.8),
        )
    end
    push!(
        pts,
        point(10.0, 2.0 * m_prongs),
        point(10.0, 2.0 * m_prongs + 1.0),
        point(0.0, 2.0 * m_prongs + 1.0),
    )
    return pts
end
P = _hcomb(1200);
n = length(P)
for alg in (:ear_cut, :ear_cut_u_kd_tree, :ear_cut_z_order)
    fill!(C, 0)
    triangulate_polygon(P; algorithm = alg)
    println(
        rpad(alg, 20),
        " calls/n=",
        round(C[1]/n, digits = 2),
        " convex-tests/n=",
        round(C[2]/n, digits = 2),
        " ears found/n=",
        round(C[3]/n, digits = 2),
        " fallbacks=",
        C[4],
    )
end
