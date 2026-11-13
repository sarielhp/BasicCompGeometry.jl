function _ear_clip_run_old!(st::_EarClipState)
    n = length(st.poly)
    ear_queue = Int[]
    sizehint!(ear_queue, 3n)
    vert_ver = ones(Int, n)
    in_queue_ver = zeros(Int, n)
    for i in 1:n
        if _is_ear_vertex(st, i)
            push!(ear_queue, i)
            in_queue_ver[i] = 1
        end
    end

    diagonals = Tuple{Int,Int}[]
    sizehint!(diagonals, n - 3)
    triangles = NTuple{3,Int}[]
    sizehint!(triangles, n - 2)

    active_count = n
    active_hint = 1
    q_pos = 1

    while active_count > 3
        if q_pos > length(ear_queue)
            fb = _find_fallback_ear(st, active_hint)
            fb == 0 && _throw_no_ear(active_count)
            push!(ear_queue, fb)
            in_queue_ver[fb] = vert_ver[fb]
        end

        ear = ear_queue[q_pos]
        q_pos += 1
        qv = in_queue_ver[ear]
        if qv == 0 || !st.is_active[ear]
            in_queue_ver[ear] = 0
            continue
        end
        if qv != vert_ver[ear] && !_is_ear_vertex(st, ear)
            in_queue_ver[ear] = 0
            continue
        end
        in_queue_ver[ear] = 0

        pr = st.prev_v[ear]
        nx = st.next_v[ear]

        push!(triangles, (pr, ear, nx))
        d = abs(pr - nx)
        if d != 1 && d != n - 1
            push!(diagonals, pr < nx ? (pr, nx) : (nx, pr))
        end

        _remove_reflex!(st.rindex, ear)
        st.next_v[pr] = nx
        st.prev_v[nx] = pr
        st.is_active[ear] = false
        active_count -= 1
        active_hint = nx

        vert_ver[pr] += 1
        vert_ver[nx] += 1
        _update_reflex!(st, pr)
        _update_reflex!(st, nx)

        for nb in (pr, nx)
            if _is_ear_vertex(st, nb)
                if in_queue_ver[nb] == 0
                    push!(ear_queue, nb)
                end
                in_queue_ver[nb] = vert_ver[nb]
            else
                in_queue_ver[nb] = 0
            end
        end
    end

    if active_count == 3
        v1 = st.is_active[active_hint] ? active_hint : findfirst(st.is_active)
        if v1 !== nothing
            v2 = st.next_v[v1]
            v3 = st.next_v[v2]
            push!(triangles, (v1, v2, v3))
        end
    end

    return PolygonTriangulation(diagonals, triangles)
end
