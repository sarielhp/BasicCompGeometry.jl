using BasicCompGeometry, Random
orient(a, b, c) = sign((b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x))
onseg(a, b, p) =
    orient(a, b, p)==0 &&
    min(a.x, b.x)<=p.x<=max(a.x, b.x) &&
    min(a.y, b.y)<=p.y<=max(a.y, b.y)
function seg_intersect(a, b, c, d)
    o1, o2, o3, o4 = orient(a, b, c), orient(a, b, d), orient(c, d, a), orient(c, d, b)
    (o1*o2<0 && o3*o4<0) && return true
    return onseg(a, b, c)||onseg(a, b, d)||onseg(c, d, a)||onseg(c, d, b)
end
function is_simple(P)
    n=length(P)
    length(unique(P))==n || return false
    for i in 1:n, j in (i+1):n
        a, b = P[i], P[mod1(i+1, n)]
        c, d = P[j], P[mod1(j+1, n)]
        if j==i+1
            (onseg(c, d, a) || onseg(a, b, d)) && return false
            continue
        elseif i==1 && j==n
            (onseg(c, d, b) || onseg(a, b, c)) && return false
            continue
        end
        seg_intersect(a, b, c, d) && return false
    end
    true
end
area2(P) =
    (n = length(P); sum(P[i].x*P[mod1(i+1, n)].y-P[mod1(i+1, n)].x*P[i].y for i in 1:n))
function pip(P, q)
    n=length(P)
    c=false
    for i in 1:n
        a, b=P[i], P[mod1(i+1, n)]
        if (a.y>q.y)!=(b.y>q.y)
            x = a.x+(q.y-a.y)*(b.x-a.x)/(b.y-a.y)
            x>q.x && (c=!c)
        end
    end
    c
end
function verify(P, res)
    n=length(P)
    length(res.triangles)==n-2 || return "ntri=$(length(res.triangles))"
    length(res.diagonals)==n-3 || return "ndiag=$(length(res.diagonals))"
    length(unique(res.diagonals))==n-3 || return "dup diag"
    for t in res.triangles
        orient(P[t[1]], P[t[2]], P[t[3]])>0 || return "non-ccw tri $t"
    end
    sum(area2([P[t[1]], P[t[2]], P[t[3]]]) for t in res.triangles)==area2(P) ||
        return "area"
    for (u, v) in res.diagonals
        a, b=P[u], P[v]
        for w in 1:n
            (w==u||w==v) && continue
            onseg(a, b, P[w]) && return "diag thru vertex"
        end
        for i in 1:n
            j=mod1(i+1, n)
            (i in (u, v) || j in (u, v)) && continue
            seg_intersect(a, b, P[i], P[j]) && return "diag crosses edge"
        end
        pip(P, (a+b)/2) || return "diag outside"
    end
    nothing
end
