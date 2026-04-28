"""
    translation(a, b)

Return a 3x3 homogeneous translation matrix for 2D coordinates.
Translates by `a` along the X-axis and `b` along the Y-axis.
"""
function translation(a::T, b::T) where {T}
    return [
        one(T) zero(T) a;
        zero(T) one(T) b;
        zero(T) zero(T) one(T)
    ]
end

"""
    translation(p::Point{2})

Return a 3x3 translation matrix using the coordinates of 2D point `p`.
"""
function translation(p::Point{2,T}) where {T}
    return translation(p[1], p[2])
end

"""
    rotation(x)

Return a 3x3 homogeneous rotation matrix for 2D coordinates.
Rotates counter-clockwise by angle `x` (in radians).
"""
function rotation(x::T) where {T}
    return [
        cos(x) -sin(x) zero(T);
        sin(x) cos(x) zero(T);
        zero(T) zero(T) one(T)
    ]
end

"""
    apply_transform(m, p::Point{2})

Apply a 3x3 homogeneous transformation matrix `m` to a 2D point `p`.
Return the resulting `Point{2}`.
"""
function apply_transform(m::Matrix{T}, p::Point{2,T})::Point{2,T} where {T}
    # Augmented vector multiplication
    res = m * [p[1], p[2], one(T)]
    return Point{2,T}(res[1], res[2])
end

export translation, rotation, apply_transform
