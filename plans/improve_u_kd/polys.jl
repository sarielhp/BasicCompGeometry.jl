using BasicCompGeometry, Random
function star(n; rng=Xoshiro(1))
    θ = sort(rand(rng, n)) .* 2π
    [(r = 0.5+rand(rng); point(r*cos(t), r*sin(t))) for t in θ]
end
function altstar(n)  # alternating radii (as in tests)
    [point((isodd(i) ? 1.0 : 0.6)*cos(2π*i/n), (isodd(i) ? 1.0 : 0.6)*sin(2π*i/n)) for i in 1:n]
end
function comb(k)
    P = Point{2,Float64}[]
    for i in 0:k-1
        push!(P, point(2.0i, 0.0)); push!(P, point(2.0i+1, 0.0))
    end
    # teeth upward
    Q = Point{2,Float64}[point(0.0,0.0), point(2.0k,0.0)]
    for i in k-1:-1:0
        push!(Q, point(2.0i+1.5, 10.0)); push!(Q, point(2.0i+1.0, 1.0)); push!(Q, point(2.0i+0.5, 10.0)); push!(Q, point(2.0i, 1.0))
    end
    Q[1:end-1]
end
function spiral(n)  # thick spiral polygon
    m = n ÷ 2; out = Point{2,Float64}[]; inn = Point{2,Float64}[]
    for i in 0:m-1
        t = 6π*i/m; r = 1 + t
        push!(out, point((r+0.8)*cos(t), (r+0.8)*sin(t))); push!(inn, point(r*cos(t), r*sin(t)))
    end
    P = vcat(out, reverse(inn))
    twice_area(P) < 0 && reverse!(P)
    P
end
# random walk "coastline": star with noisy radius random walk
function coast(n; rng=Xoshiro(2))
    r = 1.0; P = Point{2,Float64}[]
    for i in 1:n
        r = clamp(r + 0.05*randn(rng), 0.3, 3.0)
        t = 2π*i/n; push!(P, point(r*cos(t), r*sin(t)))
    end
    P
end
