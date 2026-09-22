###############################################
### Halfplane type and operations

"""
    Halfplane{T}

A closed 2D halfplane bounded by a directed infinite `Line{2, T}`.
By convention, the halfplane is the closed right halfplane of its boundary line `boundary`:
the set of all points `q` such that `p`, `p + u`, and `q` form a collinear set or a right turn
(i.e., `turn_sign(p, p + u, q) <= 0`).

The complement halfplane is obtained by reversing the line direction (`Line(p, -u)`).
"""
struct Halfplane{T}
    boundary::Line{2,T}
end

"""
    Halfplane2F

A 2-dimensional halfplane with `Float64` coordinates.
"""
const Halfplane2F = Halfplane{Float64}

# Constructors


function Halfplane(p::Point{2,T}, u::Point{2,T}) where {T}
    return Halfplane{T}(Line(p, u))
end

function Halfplane(p::Point{2,T1}, u::Point{2,T2}) where {T1,T2}
    T = promote_type(T1, T2)
    return Halfplane{T}(Line(Point{2,T}(p), Point{2,T}(u)))
end

"""
    Halfplane(line::Line{2, T}, z::Point{2, S}; atol::Real = 1e-12)

Construct the closed halfplane bounded by `line` that contains the point `z`.
The point `z` must not lie on `line` (i.e. distance from `z` to `line` must exceed `atol`);
otherwise an `ArgumentError` is thrown.
"""
function Halfplane(line::Line{2,T}, z::Point{2,S}; atol::Real = 1e-12) where {T,S}
    u_len = norm(line.u)
    u_len == 0 && error("Line direction vector u cannot be zero.")
    s = turn_sign(line.p, line.p + line.u, z)
    # Distance from point z to the line is |s| / ||u||
    d = abs(s) / u_len
    if d <= atol
        throw(ArgumentError("Point z lies on the boundary line (distance = $d <= atol = $atol); cannot determine halfplane."))
    elseif s < 0
        # z is to the right of (p -> p+u), so line's right halfplane already contains z
        return Halfplane{T}(line)
    else
        # z is to the left of (p -> p+u); reversing direction puts z into the right halfplane
        return Halfplane{T}(Line(line.p, -line.u))
    end
end

"""
    Halfplane(p::Point{2}, u::Point{2}, z::Point{2}; atol::Real = 1e-12)

Construct the closed halfplane bounded by the line passing through `p` with direction `u`
that contains the point `z`.
"""
function Halfplane(p::Point{2}, u::Point{2}, z::Point{2}; atol::Real = 1e-12)
    return Halfplane(Line(p, u), z; atol = atol)
end

# Basic accessors and properties

"""
    boundary(h::Halfplane)
    boundary_line(h::Halfplane)

Return the boundary `Line{2, T}` of the halfplane `h`.
"""
@inline boundary(h::Halfplane) = h.boundary
@inline boundary_line(h::Halfplane) = h.boundary

"""
    complement(h::Halfplane)
    -h

Return the opposite (complementary) closed halfplane by reversing the boundary line direction.
"""
@inline complement(h::Halfplane{T}) where {T} = Halfplane{T}(Line(h.boundary.p, -h.boundary.u))
@inline Base.:-(h::Halfplane) = complement(h)

function Base.show(io::IO, h::Halfplane{T}) where {T}
    print(io, "Halfplane(Line(p=", h.boundary.p, ", u=", h.boundary.u, "))")
end

# Point membership & predicates

"""
    is_inside(q::Point{2}, h::Halfplane; atol::Real = 0.0)
    is_inside(h::Halfplane, q::Point{2}; atol::Real = 0.0)

Return `true` if query point `q` lies inside the closed halfplane `h`.
Under the right-halfplane convention, this checks if `p`, `p + u`, and `q` form
a right turn or are collinear (`turn_sign(p, p + u, q) <= atol`).
"""
@inline function is_inside(q::Point{2}, h::Halfplane; atol::Real = 0.0)
    p = h.boundary.p
    return turn_sign(p, p + h.boundary.u, q) <= atol
end

@inline is_inside(h::Halfplane, q::Point{2}; atol::Real = 0.0) = is_inside(q, h; atol = atol)

@inline Base.in(q::Point{2}, h::Halfplane) = is_inside(q, h)

"""
    in_interior(q::Point{2}, h::Halfplane; atol::Real = 0.0)

Return `true` if query point `q` lies strictly in the interior of the halfplane (not on boundary).
"""
@inline function in_interior(q::Point{2}, h::Halfplane; atol::Real = 0.0)
    p = h.boundary.p
    return turn_sign(p, p + h.boundary.u, q) < -atol
end

"""
    on_boundary(q::Point{2}, h::Halfplane; atol::Real = 1e-12)

Return `true` if query point `q` lies on the boundary line of `h` within tolerance `atol`.
"""
@inline function on_boundary(q::Point{2}, h::Halfplane; atol::Real = 1e-12)
    u_len = norm(h.boundary.u)
    u_len == 0 && return false
    s = turn_sign(h.boundary.p, h.boundary.p + h.boundary.u, q)
    return abs(s) / u_len <= atol
end

"""
    distance(q::Point{2}, h::Halfplane)
    distance(h::Halfplane, q::Point{2})

Return the Euclidean distance from query point `q` to the halfplane `h`.
If `q` is inside `h`, the distance is `0.0`. Otherwise, it is the distance to the boundary line.
"""
function distance(q::Point{2,S}, h::Halfplane{T}) where {S,T}
    if is_inside(q, h)
        return zero(promote_type(S, T, Float64))
    end
    return distance(q, h.boundary)
end

@inline distance(h::Halfplane, q::Point{2}) = distance(q, h)

# Depth operation

"""
    depth(halfplanes, q::Point{2})
    depth(q::Point{2}, halfplanes)

Count how many halfplanes in `halfplanes` contain the query point `q`.
"""
function depth(
    halfplanes::Union{AbstractVector{<:Halfplane}, Tuple{Vararg{Halfplane}}},
    q::Point{2},
)
    return count(h -> is_inside(q, h), halfplanes)
end

@inline function depth(
    q::Point{2},
    halfplanes::Union{AbstractVector{<:Halfplane}, Tuple{Vararg{Halfplane}}},
)
    return depth(halfplanes, q)
end

# Random halfplane generation

"""
    rand_halfplane([rng=Random.default_rng()])

Generate a random 2D halfplane by:
1. Picking a point `p` uniformly at random in the unit square `[0, 1]^2`.
2. Picking a random direction `u` from the standard 2D normal distribution and normalizing it to unit length.
3. Choosing with equal probability (50%) one of the two halfplanes bounded by the resulting line (by randomly negating `u`).
"""
function rand_halfplane(rng::AbstractRNG = Random.default_rng())
    p = Point{2,Float64}(rand(rng, Float64), rand(rng, Float64))
    dx = randn(rng, Float64)
    dy = randn(rng, Float64)
    while hypot(dx, dy) < 1e-14
        dx = randn(rng, Float64)
        dy = randn(rng, Float64)
    end
    l = hypot(dx, dy)
    u = Point{2,Float64}(dx / l, dy / l)
    if rand(rng, Bool)
        u = -u
    end
    return Halfplane(Line(p, u))
end

rand_halfplane(n::Integer; rng::AbstractRNG = Random.default_rng()) = [rand_halfplane(rng) for _ = 1:n]

"""
    random_halfplane([rng=Random.default_rng()])
    random_halfplane(n::Integer; [rng=Random.default_rng()])

Alias for `rand_halfplane`.
"""
const random_halfplane = rand_halfplane

Base.rand(rng::AbstractRNG, ::Type{Halfplane}) = rand_halfplane(rng)
Base.rand(::Type{Halfplane}) = rand_halfplane(Random.default_rng())
Base.rand(rng::AbstractRNG, ::Type{Halfplane2F}) = rand_halfplane(rng)
Base.rand(::Type{Halfplane2F}) = rand_halfplane(Random.default_rng())

Base.rand(rng::AbstractRNG, ::Type{Halfplane}, n::Integer) = [rand_halfplane(rng) for _ = 1:n]
Base.rand(::Type{Halfplane}, n::Integer) = rand(Random.default_rng(), Halfplane, n)
Base.rand(rng::AbstractRNG, ::Type{Halfplane2F}, n::Integer) = [rand_halfplane(rng) for _ = 1:n]
Base.rand(::Type{Halfplane2F}, n::Integer) = rand(Random.default_rng(), Halfplane2F, n)

# Sutherland-Hodgman polygon clipping by a halfplane

"""
    clip(poly::Vector{Point{2, T}}, h::Halfplane; eps::Real = 1e-11)

Clip a 2D polygon `poly` (given as a vertex sequence) by the halfplane `h` using the
Sutherland-Hodgman algorithm. Returns a new `Vector{Point{2, T}}`.
"""
function clip(poly::Vector{Point{2,T}}, h::Halfplane; eps::Real = 1e-11) where {T}
    n = length(poly)
    n == 0 && return Point{2,T}[]
    output = Point{2,T}[]
    p_line = h.boundary.p
    u_line = h.boundary.u
    u_len = norm(u_line)
    u_len == 0 && error("Boundary line direction is zero.")

    for i = 1:n
        cur = poly[i]
        prev = poly[i == 1 ? n : i - 1]

        cur_dist_signed = turn_sign(p_line, p_line + u_line, cur) / u_len
        prev_dist_signed = turn_sign(p_line, p_line + u_line, prev) / u_len

        cur_in = cur_dist_signed <= eps
        prev_in = prev_dist_signed <= eps

        if cur_in
            if !prev_in
                pt = _intersect_segment_line(prev, cur, h.boundary)
                pt !== nothing && push!(output, pt)
            end
            push!(output, cur)
        elseif prev_in
            pt = _intersect_segment_line(prev, cur, h.boundary)
            pt !== nothing && push!(output, pt)
        end
    end
    return output
end

function _intersect_segment_line(
    a::Point{2,T},
    b::Point{2,T},
    line::Line{2},
)::Union{Nothing,Point{2,T}} where {T}
    dp = b - a
    det_m = line.u[1] * (-dp[2]) - line.u[2] * (-dp[1])
    if abs(det_m) < 1e-14
        return nothing
    end
    diff = a - line.p
    t = (line.u[1] * diff[2] - line.u[2] * diff[1]) / det_m
    t_clamped = clamp(t, 0.0, 1.0)
    return Point{2,T}((1 - t_clamped) * a + t_clamped * b)
end

function clip(poly::PntSeq{2,T}, h::Halfplane; eps::Real = 1e-11) where {T}
    return PntSeq(clip(poly.pnts, h; eps = eps))
end

"""
    write_halfplanes(filename::String, hps::AbstractVector{<:Halfplane})

Save a collection of halfplanes to a text file. Each line contains four coordinates:
`p_x p_y u_x u_y` representing `Line(Point(p_x, p_y), Point(u_x, u_y))`.
"""
function write_halfplanes(filename::String, hps::AbstractVector{<:Halfplane})
    mkpath(dirname(abspath(filename)))
    open(filename, "w") do io
        println(io, "# Halfplanes count: ", length(hps))
        println(io, "# Format: p_x p_y u_x u_y")
        for h in hps
            p = h.boundary.p
            u = h.boundary.u
            @printf(io, "%.17g %.17g %.17g %.17g\n", p[1], p[2], u[1], u[2])
        end
    end
    return filename
end

"""
    read_halfplanes(filename::String)::Vector{Halfplane2F}

Read a collection of halfplanes from a text file formatted with `p_x p_y u_x u_y` per line.
Lines starting with `#` and empty lines are ignored.
"""
function read_halfplanes(filename::String)::Vector{Halfplane2F}
    hps = Halfplane2F[]
    open(filename, "r") do io
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            startswith(s, "#") && continue
            tokens = split(s)
            if length(tokens) >= 4
                px = parse(Float64, tokens[1])
                py = parse(Float64, tokens[2])
                ux = parse(Float64, tokens[3])
                uy = parse(Float64, tokens[4])
                push!(hps, Halfplane(Line(Point2F(px, py), Point2F(ux, uy))))
            end
        end
    end
    return hps
end

export Halfplane, Halfplane2F
export boundary, boundary_line, complement, depth
export rand_halfplane, random_halfplane
export in_interior, on_boundary, clip
export write_halfplanes, read_halfplanes
