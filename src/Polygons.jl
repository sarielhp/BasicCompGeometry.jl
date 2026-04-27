###############################################
### Polygon type

"""
    Polygon{D, T}

A wrapper around `Vector{Point{D, T}}` representing a polygon or a polygonal curve in D dimensions.
The points are stored sequentially in the `pnts` field.
"""
struct Polygon{D, T}
    pnts::Vector{Point{D, T}}
end

# Constructors
"""
    Polygon{D, T}()

Construct an empty polygon for dimension `D` and type `T`.
"""
Polygon{D, T}() where {D, T} = Polygon{D, T}(Point{D, T}[])

# Note: Struct definition already provides: Polygon(pnts::Vector{Point{D, T}})

# Accessors
"""
    Points(P::Polygon)

Return the underlying vector of points of the polygon.
"""
@inline Points(P::Polygon) = P.pnts

"""
    cardin(P)

Return the number of vertices (cardinality) in the polygon `P`.
"""
@inline cardin(P::Polygon) = length(P.pnts)

"""
    translate!(P, v)

Translate the polygon `P` in-place by the vector `v`.
"""
function translate!(P::Polygon{D, T}, v::Point{D, T}) where {D, T}
    P.pnts .+= Ref(v)
    return P
end

"""
    scale!(P, s)

Scale the polygon `P` in-place by factor `s` relative to the origin.
"""
function scale!(P::Polygon{D, T}, s::Real) where {D, T}
    P.pnts .*= s
    return P
end

"""
    geom_length(poly)

Calculate the total Euclidean length of the polygonal curve (sum of edge lengths).
"""
function geom_length(poly::Polygon{D, T}) where {D, T}
    len = 0.0
    pnts = poly.pnts
    for i in 1:length(pnts)-1
        len += dist(pnts[i], pnts[i+1])
    end
    return len
end

"""
    prefix_lengths(poly)

Return a vector containing the cumulative arc-length distances at each vertex of the polygon.
The first element is always 0.0, and the last is the total geometric length.
"""
function prefix_lengths(poly::Polygon{D, T}) where {D, T}
    n = length(poly)
    v = Vector{Float64}(undef, n)
    n == 0 && return v
    
    v[1] = 0.0
    pnts = poly.pnts
    for i in 2:n
        v[i] = v[i-1] + dist(pnts[i-1], pnts[i])
    end
    return v
end

"""
    at(poly, t)

Return a point on the polygonal curve at the normalized parameter `t` in `[0, 1]`.
`t=0` is the first vertex, `t=1` is the last vertex.
The parameterization is based on arc-length (linear interpolation along edges).
"""
function at(poly::Polygon{D, T}, t::Real) where {D, T}
    n = length(poly)
    n == 0 && throw(ArgumentError("Empty polygon"))
    n == 1 && return poly.pnts[1]
    
    plens = prefix_lengths(poly)
    total_len = plens[end]
    total_len == 0.0 && return poly.pnts[1]
    
    # Target absolute arc-length
    target = t * total_len
    if target <= 0
        return poly.pnts[1]
    elseif target >= total_len
        return poly.pnts[n]
    end
    
    # Find the segment containing the target distance
    i = searchsortedfirst(plens, target)
    i <= 1 && return poly.pnts[1]
    
    segment_len = plens[i] - plens[i-1]
    local_t = (target - plens[i-1]) / segment_len
    return convex_comb(poly.pnts[i-1], poly.pnts[i], local_t)
end

# Base Extensions for collection-like behavior
Base.length(P::Polygon) = length(P.pnts)
Base.getindex(P::Polygon, i::Int) = P.pnts[i]
Base.setindex!(P::Polygon, v, i::Int) = (P.pnts[i] = v)
Base.firstindex(P::Polygon) = 1
Base.lastindex(P::Polygon) = length(P.pnts)
Base.iterate(P::Polygon) = iterate(P.pnts)
Base.iterate(P::Polygon, state) = iterate(P.pnts, state)
Base.push!(P::Polygon, p::Point) = push!(P.pnts, p)
Base.pop!(P::Polygon) = pop!(P.pnts)
Base.empty!(P::Polygon) = empty!(P.pnts)
Base.first(P::Polygon) = first(P.pnts)
Base.last(P::Polygon) = last(P.pnts)
Base.size(P::Polygon) = size(P.pnts)

function Base.show(io::IO, poly::Polygon{D, T}) where {D, T}
    print(io, "Polygon{$D,$T} with $(length(poly)) points")
end

"""
    Matrix(P::Polygon)

Convert the polygon's vertices into a `D x N` matrix where `N` is the number of points.
"""
function Base.Matrix(P::Polygon{D, T}) where {D, T}
    m = Matrix{T}(undef, D, length(P))
    for (i, p) in enumerate(P.pnts)
        m[:, i] = p
    end
    return m
end

"""
    write_plt(filename, P)

Write the polygon coordinates to a PLT-compatible delimited file (N rows, D columns).
"""
function write_plt(filename::String, P::Polygon)
    open(filename, "w") do io
        writedlm(io, Matrix(P)', ',')
    end
end

"""
    simplify(P, r)

Perform distance-based simplification of the polygon `P`.
Vertices are skipped if they are within distance `r` of the last kept vertex.
Returns `(simplified_polygon, kept_indices)`.
"""
function simplify(P::Polygon{D, T}, r::Real) where {D, T}
    n = length(P)
    n <= 2 && return P, collect(1:n)
    
    pout = Polygon{D, T}()
    pindices = Int[]
    
    push!(pout, P[1])
    push!(pindices, 1)
    
    curr = P[1]
    for i in 2:n-1
        if dist(P[i], curr) > r
            push!(pout, P[i])
            push!(pindices, i)
            curr = P[i]
        end
    end
    
    push!(pout, P[n])
    push!(pindices, n)
    return pout, pindices
end

"""
    sample_uniformly(P, n)

Resample the polygon `P` into `n` vertices that are evenly spaced by arc-length.
"""
function sample_uniformly(P::Polygon{D, T}, n::Int) where {D, T}
    n <= 1 && throw(ArgumentError("n must be > 1"))
    new_P = Polygon{D, T}()
    for i in 1:n
        push!(new_P, at(P, (i-1)/(n-1)))
    end
    return new_P
end

"""
    rand_polygon(D, T, n)

Generate a random polygon in D dimensions with `n` vertices, each coordinate in `[0, 1)`.
"""
function rand_polygon(D::Int, T::Type{<:Real}, n::Int)
    return Polygon([rand_point(D, T) for _ in 1:n])
end

"""
    Polygon_random_sphere(D, T, n)

Generate a random polygon in D dimensions with `n` vertices sampled uniformly from the unit sphere.
"""
function Polygon_random_sphere(D::Int, T::Type{<:Real}, n::Int)
    pnts = Vector{Point{D, T}}(undef, n)
    for i in 1:n
        p = randn(T, D)
        mag = norm(p)
        if mag > 0
            p ./= mag
        end
        pnts[i] = Point{D, T}(p)
    end
    return Polygon(pnts)
end

"""
    Polygon2I

Alias for a 2D polygon with `Int64` coordinates.
"""
const Polygon2I = Polygon{2, Int64}

"""
    Polygon2F

Alias for a 2D polygon with `Float64` coordinates.
"""
const Polygon2F = Polygon{2, Float64}

"""
    Polygon3F

Alias for a 3D polygon with `Float64` coordinates.
"""
const Polygon3F = Polygon{3, Float64}

export Polygon, Polygon2I, Polygon2F, Polygon3F
export Points, cardin, geom_length, prefix_lengths
export translate!, scale!, simplify, sample_uniformly, at
export rand_polygon, Polygon_random_sphere, write_plt
