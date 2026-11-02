#!/usr/bin/env julia

# Discrete hull of a large disk: ConvexHull(D ∩ Z²), computed by the
# output-sensitive algorithm of Har-Peled (1998) in O(r^(2/3) log r) time,
# checked against the naive O(r) column-extremes hull, plus a small figure.

include(joinpath(@__DIR__, "support/runtime.jl"))

using QuickEnv
using BasicCompGeometry
using BasicCompGeometry.IpeDraw
using Printf
using Random

Random.seed!(1998)
lcross(u, v) = Int128(u[1]) * v[2] - Int128(u[2]) * v[1]

# Lattice points on the boundary of the hull polygon (Table 4.3 of the paper).
boundary_points(H) = sum(gcd(H[mod1(i + 1, end)] - H[i]...) for i in eachindex(H))

is_strictly_convex(H) =
    all(lcross(H[mod1(i + 1, end)] - H[i], H[mod1(i + 2, end)] - H[mod1(i + 1, end)]) > 0 for
        i in eachindex(H))

# 1. One large disk, verified against the naive algorithm.
r = example_size(10_000_000, 100_000)
D = LatticeDisk((rand(), rand()), r)
discrete_hull(LatticeDisk((0.5, 0.5), 10.0))               # compile
t_fast = @elapsed H = discrete_hull(D)
t_naive = @elapsed H_naive = discrete_hull_columns(D)
ok = H == H_naive && all(p -> p ∈ D, H) && is_strictly_convex(H)
@printf("Disk: center (%.4f, %.4f), radius %d\n", BasicCompGeometry.center(D)..., r)
@printf("  hull vertices        %d  (%.3f r^(2/3))\n", length(H), length(H) / r^(2 / 3))
@printf("  boundary lattice pts %d  (%.3f r^(2/3))\n", boundary_points(H),
        boundary_points(H) / r^(2 / 3))
@printf("  discrete_hull        %.4f s\n", t_fast)
@printf("  naive column hull    %.4f s\n", t_naive)
println(ok ? "  Verified: identical to the naive hull, all vertices in D, strictly convex." :
        "  VERIFICATION FAILED")
ok || exit(1)

# 2. Average size over random centers (cf. Tables 4.1-4.2: ratio ≈ 3.45).
println("\n        r   avg |DH|   |DH|/r^(2/3)   avg time (s)")
trials = example_size(20, 5)
for k = 1:example_size(7, 4)
    rk = 10^k
    sizes, time = Int[], 0.0
    for _ = 1:trials
        Dk = LatticeDisk((rand(), rand()), rk)
        time += @elapsed Hk = discrete_hull(Dk)
        push!(sizes, length(Hk))
    end
    avg = sum(sizes) / trials
    @printf("%9d %10.1f %14.3f %14.6f\n", rk, avg, avg / rk^(2 / 3), time / trials)
end

# 3. Figure: a small disk, its lattice points, and its discrete hull.
output_dir = example_output()
mkpath(output_dir)
pdf_path = joinpath(output_dir, "discrete_hull_disk.pdf")
small = LatticeDisk((0.37, 0.21), 7.3)
c, rs = BasicCompGeometry.center(small), BasicCompGeometry.radius(small)
Hs = discrete_hull(small)
grid = [Point2I(x, y) for x = floor(Int, c.x - rs)-1:ceil(Int, c.x + rs)+1
        for y = floor(Int, c.y - rs)-1:ceil(Int, c.y + rs)+1]
inside = [Point2F(p) for p in grid if p ∈ small]
outside = [Point2F(p) for p in grid if p ∉ small]
disk = Circle(Point2F(c), rs)
hull_poly = Polygon([Point2F(p) for p in Hs])

theme = Theme(
    disk = Style(stroke = :blue, pen = :heavier),
    hull = Style(fill = :lightred, stroke = :red, fill_opacity = 0.3, pen = :heavier),
)
view = BBox(Point2F(c.x - rs - 1, c.y - rs - 1), Point2F(c.x + rs + 1, c.y + rs + 4))
figure(pdf_path; fit = view, margin = 20, preview = "--preview" in ARGS) do fig
    draw!(fig, hull_poly; style = theme.hull)
    draw!(fig, disk; style = theme.disk)
    mark!(fig, outside; stroke = :gray, fill = :gray, shape = :disk, size = :tiny)
    mark!(fig, inside; stroke = :black, fill = :black, shape = :disk, size = :small)
    mark!(fig, [Point2F(p) for p in Hs]; stroke = :red, fill = :red, shape = :disk,
          size = :normal)
    legend!(fig, [raw"disk $D$" => theme.disk,
                  raw"discrete hull $\mathcal{CH}(D \cap \mathbb{Z}^2)$" => theme.hull];
            position = :northwest, width = 170)
end
println("\nGenerated: ", pdf_path)
