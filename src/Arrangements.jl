###############################################
### Arrangement Vertical Decomposition (Trapezoidal Map)

"""
    VerticalTrapezoid{T}

A vertical trapezoid in 2D bounded by:
- `x_left::T`: left vertical boundary coordinate
- `x_right::T`: right vertical boundary coordinate
- `bottom::Union{Nothing, Line{2, T}}`: bottom boundary line (or `nothing` for bottom of view box)
- `top::Union{Nothing, Line{2, T}}`: top boundary line (or `nothing` for top of view box)
- `corners::SVector{4, Point{2, T}}`: 4 corner vertices ordered counter-clockwise:
  1. bottom-left: `(x_left, y_bottom(x_left))`
  2. bottom-right: `(x_right, y_bottom(x_right))`
  3. top-right: `(x_right, y_top(x_right))`
  4. top-left: `(x_left, y_top(x_left))`
"""
struct VerticalTrapezoid{T}
    x_left::T
    x_right::T
    bottom::Union{Nothing,Line{2,T}}
    top::Union{Nothing,Line{2,T}}
    corners::SVector{4,Point{2,T}}
end

"""
    VertTrapezoid

Alias for `VerticalTrapezoid`.
"""
const VertTrapezoid = VerticalTrapezoid

"""
    VerticalTrapezoid2F

A vertical trapezoid with `Float64` coordinates.
"""
const VerticalTrapezoid2F = VerticalTrapezoid{Float64}

"""
    area(trap::VerticalTrapezoid)

Return the area of the vertical trapezoid.
"""
function area(trap::VerticalTrapezoid)
    w = trap.x_right - trap.x_left
    h_left = trap.corners[4][2] - trap.corners[1][2]
    h_right = trap.corners[3][2] - trap.corners[2][2]
    return w * (h_left + h_right) / 2.0
end

"""
    vertices(trap::VerticalTrapezoid)

Return the corner vertices of the vertical trapezoid as an `SVector{4, Point{2, T}}`.
"""
vertices(trap::VerticalTrapezoid) = trap.corners

"""
    polygon(trap::VerticalTrapezoid)

Return the corner vertices of the trapezoid as a `PntSeq{2, T}`.
"""
polygon(trap::VerticalTrapezoid{T}) where {T} = PntSeq(Vector(trap.corners))

"""
    is_inside(p::Point{2}, trap::VerticalTrapezoid; tol::Real = 1e-11)

Return `true` if query point `p` lies inside or on the boundary of `trap`.
"""
function is_inside(p::Point{2}, trap::VerticalTrapezoid; tol::Real = 1e-11)
    if p[1] < trap.x_left - tol || p[1] > trap.x_right + tol
        return false
    end
    w = trap.x_right - trap.x_left
    if w <= tol
        return min(trap.corners[1][2], trap.corners[4][2]) - tol <=
               p[2] <=
               max(trap.corners[1][2], trap.corners[4][2]) + tol
    end
    t = (p[1] - trap.x_left) / w
    y_bot = (1 - t) * trap.corners[1][2] + t * trap.corners[2][2]
    y_top = (1 - t) * trap.corners[4][2] + t * trap.corners[3][2]
    return y_bot - tol <= p[2] <= y_top + tol
end

Base.in(p::Point{2}, trap::VerticalTrapezoid) = is_inside(p, trap)

function Base.show(io::IO, trap::VerticalTrapezoid{T}) where {T}
    print(
        io,
        "VerticalTrapezoid(x=[",
        round(trap.x_left, digits = 4),
        "..",
        round(trap.x_right, digits = 4),
        "], ",
    )
    print(
        io,
        "y_left=[",
        round(trap.corners[1][2], digits = 4),
        "..",
        round(trap.corners[4][2], digits = 4),
        "], ",
    )
    print(
        io,
        "y_right=[",
        round(trap.corners[2][2], digits = 4),
        "..",
        round(trap.corners[3][2], digits = 4),
        "])",
    )
end

###############################################
### ArrVertDecomp type

"""
    ArrVertDecomp{T}

Vertical decomposition of an arrangement of lines inside a view rectangle.

Fields:
- `lines::Vector{Line{2, T}}`: original bundle of lines (or boundary lines of halfplanes)
- `vertices::Vector{Point{2, T}}`: arrangement vertices (intersections) inside the view
- `view::BBox{2, T}`: view bounding box
- `trapezoids::Vector{VerticalTrapezoid{T}}`: computed vertical trapezoids
"""
struct ArrVertDecomp{T}
    lines::Vector{Line{2,T}}
    vertices::Vector{Point{2,T}}
    view::BBox{2,T}
    trapezoids::Vector{VerticalTrapezoid{T}}
end

"""
    ArrVertDecomp2F

Alias for `ArrVertDecomp{Float64}`.
"""
const ArrVertDecomp2F = ArrVertDecomp{Float64}

# Accessors
trapezoids(decomp::ArrVertDecomp) = decomp.trapezoids
vertices(decomp::ArrVertDecomp) = decomp.vertices
view_box(decomp::ArrVertDecomp) = decomp.view
lines(decomp::ArrVertDecomp) = decomp.lines
area(decomp::ArrVertDecomp) = sum(area, decomp.trapezoids)

Base.length(decomp::ArrVertDecomp) = length(decomp.trapezoids)
Base.getindex(decomp::ArrVertDecomp, i::Int) = decomp.trapezoids[i]
Base.iterate(decomp::ArrVertDecomp, state...) = iterate(decomp.trapezoids, state...)

"""
    locate(p::Point{2}, decomp::ArrVertDecomp; tol::Real = 1e-11)

Locate the vertical trapezoid in `decomp` containing query point `p`.
Returns `VerticalTrapezoid` or `nothing` if `p` is outside the view box.
"""
function locate(p::Point{2}, decomp::ArrVertDecomp; tol::Real = 1e-11)
    !is_inside(p, decomp.view) && return nothing
    for trap in decomp.trapezoids
        if is_inside(p, trap; tol = tol)
            return trap
        end
    end
    return nothing
end

function Base.show(io::IO, decomp::ArrVertDecomp{T}) where {T}
    print(
        io,
        "ArrVertDecomp{$T}(lines=",
        length(decomp.lines),
        ", vertices=",
        length(decomp.vertices),
        ", trapezoids=",
        length(decomp.trapezoids),
        ")",
    )
end

###############################################
### Algorithm: vertical_decomposition

@inline function _line_y_at_x(l::Line{2}, x::Real)
    return l.p[2] + (l.u[2] / l.u[1]) * (x - l.p[1])
end

struct _VertWall{T}
    y_bot::T
    y_top::T
end

mutable struct _ActiveTrap{T}
    x_left::T
    x_right::T
    bottom::Union{Nothing,Line{2,T}}
    top::Union{Nothing,Line{2,T}}
end

# Intersection vertices inside the view (clamped to it), deduplicated.
function _vd_in_view_vertices(
    lines::AbstractVector{Line{2,R}},
    x_min,
    x_max,
    y_min,
    y_max,
    tol,
) where {R}
    in_view_verts = Point{2,R}[]
    for v in vertices(lines; tol = tol)
        if (x_min - tol <= v[1] <= x_max + tol) && (y_min - tol <= v[2] <= y_max + tol)
            push!(
                in_view_verts,
                Point{2,R}(clamp(v[1], x_min, x_max), clamp(v[2], y_min, y_max)),
            )
        end
    end
    return _unique_points(in_view_verts; tol = tol)
end

# Push x-coordinates where `l` crosses the top/bottom view sides, or where it is vertical.
function _vd_push_line_events!(x_events, l, x_min, x_max, y_min, y_max)
    if abs(l.u[2]) > 1e-14
        for y_b in (y_min, y_max)
            t = (y_b - l.p[2]) / l.u[2]
            x_b = l.p[1] + t * l.u[1]
            if x_min < x_b < x_max
                push!(x_events, x_b)
            end
        end
    end
    # Vertical lines (u_x == 0)
    if abs(l.u[1]) <= 1e-14
        x_v = l.p[1]
        if x_min < x_v < x_max
            push!(x_events, x_v)
        end
    end
    return nothing
end

# Sorted, deduplicated x-coordinates of all slab boundaries; first is `x_min`, last `x_max`.
function _vd_event_xs(
    lines::AbstractVector{Line{2,R}},
    in_view_verts,
    x_min,
    x_max,
    y_min,
    y_max,
    tol,
) where {R}
    x_events = R[x_min, x_max]
    for v in in_view_verts
        push!(x_events, v[1])
    end

    # Also add x-coordinates where lines cross y = y_min or y = y_max
    for l in lines
        _vd_push_line_events!(x_events, l, x_min, x_max, y_min, y_max)
    end
    return _vd_dedup_xs(x_events, x_min, x_max, tol)
end

# Sort and deduplicate slab boundaries, forcing the extremes to `x_min` and `x_max`.
function _vd_dedup_xs(x_events::Vector{R}, x_min, x_max, tol) where {R}
    sort!(x_events)
    unique_x = R[]
    for x in x_events
        if isempty(unique_x) || x - unique_x[end] > tol
            push!(unique_x, x)
        end
    end
    if isempty(unique_x) || unique_x[1] > x_min
        pushfirst!(unique_x, x_min)
    else
        unique_x[1] = x_min
    end
    if unique_x[end] < x_max
        push!(unique_x, x_max)
    else
        unique_x[end] = x_max
    end
    return unique_x
end

# Vertical walls at each x-event: through every vertex, and along every vertical line.
function _vd_walls(
    lines::AbstractVector{Line{2,R}},
    in_view_verts,
    unique_x,
    y_min,
    y_max,
    tol,
) where {R}
    walls_at_event = [_VertWall{R}[] for _ = 1:length(unique_x)]

    for v in in_view_verts
        k_idx = findfirst(x -> abs(x - v[1]) <= tol, unique_x)
        k_idx === nothing && continue
        x_val = unique_x[k_idx]

        non_vert = filter(l -> abs(l.u[1]) > 1e-14, lines)
        y_vals = [_line_y_at_x(l, x_val) for l in non_vert]

        above_y = filter(y -> y > v[2] + tol, y_vals)
        below_y = filter(y -> y < v[2] - tol, y_vals)

        wall_top = isempty(above_y) ? y_max : min(minimum(above_y), y_max)
        wall_bot = isempty(below_y) ? y_min : max(maximum(below_y), y_min)

        push!(walls_at_event[k_idx], _VertWall{R}(wall_bot, wall_top))
    end

    # Vertical lines (u_x == 0) have vertical wall spanning [y_min, y_max]
    for l in lines
        abs(l.u[1]) <= 1e-14 || continue
        k_idx = findfirst(x -> abs(x - l.p[1]) <= tol, unique_x)
        k_idx === nothing && continue
        push!(walls_at_event[k_idx], _VertWall{R}(y_min, y_max))
    end
    return walls_at_event
end

# y of `line` at `x` clamped to the view; `default` stands for a missing (view) boundary.
_vd_clamped_y(line, x, default, y_min, y_max) =
    line === nothing ? default : clamp(_line_y_at_x(line, x), y_min, y_max)

function _vd_make_trap(at::_ActiveTrap, ::Type{R}, y_min, y_max) where {R}
    c1_y = _vd_clamped_y(at.bottom, at.x_left, y_min, y_min, y_max)
    c2_y = _vd_clamped_y(at.bottom, at.x_right, y_min, y_min, y_max)
    c3_y = _vd_clamped_y(at.top, at.x_right, y_max, y_min, y_max)
    c4_y = _vd_clamped_y(at.top, at.x_left, y_max, y_min, y_max)

    c1 = Point{2,R}(at.x_left, c1_y)
    c2 = Point{2,R}(at.x_right, c2_y)
    c3 = Point{2,R}(at.x_right, c3_y)
    c4 = Point{2,R}(at.x_left, c4_y)
    return VerticalTrapezoid{R}(
        at.x_left,
        at.x_right,
        at.bottom,
        at.top,
        SVector{4,Point{2,R}}(c1, c2, c3, c4),
    )
end

# Non-coincident lines crossing the interior of slab [x_l, x_r], sorted bottom to top.
function _vd_slab_lines(lines::AbstractVector{Line{2,R}}, x_l, x_r, y_min, y_max, tol) where {R}
    x_mid = (x_l + x_r) / 2.0
    slab_lines = Line{2,R}[]
    for l in lines
        abs(l.u[1]) <= 1e-14 && continue
        y_l = _line_y_at_x(l, x_l)
        y_r = _line_y_at_x(l, x_r)
        if max(y_l, y_r) > y_min + tol && min(y_l, y_r) < y_max - tol
            push!(slab_lines, l)
        end
    end

    sort!(slab_lines, by = l -> _line_y_at_x(l, x_mid))

    # Deduplicate coincident lines in slab
    unique_slab = Line{2,R}[]
    for l in slab_lines
        if isempty(unique_slab) ||
           abs(_line_y_at_x(l, x_mid) - _line_y_at_x(unique_slab[end], x_mid)) > tol
            push!(unique_slab, l)
        end
    end
    return unique_slab
end

# Can active trapezoid `at` (bounded by `bot_b`, `top_b`) continue across the wall at `x_l`?
function _vd_can_extend(at, bot_b, top_b, x_l, walls, y_min, y_max, tol)
    (at.bottom == bot_b && at.top == top_b) || return false
    y_bot_xl = _vd_clamped_y(bot_b, x_l, y_min, y_min, y_max)
    y_top_xl = _vd_clamped_y(top_b, x_l, y_max, y_min, y_max)
    (y_top_xl - y_bot_xl) <= tol && return false
    for w in walls
        if max(w.y_bot, y_bot_xl) < min(w.y_top, y_top_xl) - tol
            return false
        end
    end
    return true
end

function _vd_find_active(active_traps, k, bot_b, top_b, x_l, walls, y_min, y_max, tol)
    k > 1 || return nothing
    for (a_idx, at) in enumerate(active_traps)
        if _vd_can_extend(at, bot_b, top_b, x_l, walls, y_min, y_max, tol)
            return a_idx
        end
    end
    return nothing
end

"""
    vertical_decomposition(bundle, view; tol::Real = 1e-11)
    vertical_decomposition(bundle; factor::Real = 1.1, tol::Real = 1e-11)

Compute the vertical decomposition (trapezoidal map) of the arrangement of lines
inside `view` (a `BBox{2}`).
`bundle` can be an `AbstractVector{<:Line{2}}` or `AbstractVector{<:Halfplane}`.
If `view` is omitted, `default_view(bundle; factor=factor)` is used.
Returns an `ArrVertDecomp`.
"""
function vertical_decomposition(
    bundle_lines::AbstractVector{<:Line{2,T}},
    view::BBox{2,S};
    tol::Real = 1e-11,
) where {T,S}
    R_type = promote_type(T, S, Float64)
    x_min = R_type(view.mini[1])
    x_max = R_type(view.maxi[1])
    y_min = R_type(view.mini[2])
    y_max = R_type(view.maxi[2])

    if x_min >= x_max || y_min >= y_max
        throw(ArgumentError("View bounding box must have positive width and height: $view"))
    end

    lines_typed =
        Line{2,R_type}[Line{2,R_type}(Point{2,R_type}(l.p), Point{2,R_type}(l.u)) for l in bundle_lines]

    in_view_verts = _vd_in_view_vertices(lines_typed, x_min, x_max, y_min, y_max, tol)
    unique_x = _vd_event_xs(lines_typed, in_view_verts, x_min, x_max, y_min, y_max, tol)
    walls_at_event = _vd_walls(lines_typed, in_view_verts, unique_x, y_min, y_max, tol)

    # Sweep slabs and build/merge trapezoids
    completed = VerticalTrapezoid{R_type}[]
    active_traps = _ActiveTrap{R_type}[]

    for k = 1:(length(unique_x)-1)
        x_l = unique_x[k]
        x_r = unique_x[k+1]
        slab_lines = _vd_slab_lines(lines_typed, x_l, x_r, y_min, y_max, tol)

        bounds = Vector{Union{Nothing,Line{2,R_type}}}()
        push!(bounds, nothing) # bottom of view
        append!(bounds, slab_lines)
        push!(bounds, nothing) # top of view

        new_active = _ActiveTrap[]
        for p_idx = 1:(length(bounds)-1)
            bot_b = bounds[p_idx]
            top_b = bounds[p_idx+1]
            matched_idx = _vd_find_active(
                active_traps,
                k,
                bot_b,
                top_b,
                x_l,
                walls_at_event[k],
                y_min,
                y_max,
                tol,
            )
            if matched_idx !== nothing
                at = active_traps[matched_idx]
                at.x_right = x_r
                push!(new_active, at)
                deleteat!(active_traps, matched_idx)
            else
                push!(new_active, _ActiveTrap(x_l, x_r, bot_b, top_b))
            end
        end

        for at in active_traps
            push!(completed, _vd_make_trap(at, R_type, y_min, y_max))
        end
        active_traps = new_active
    end

    for at in active_traps
        push!(completed, _vd_make_trap(at, R_type, y_min, y_max))
    end

    return ArrVertDecomp{R_type}(lines_typed, in_view_verts, view, completed)
end

function vertical_decomposition(
    halfplanes::AbstractVector{<:Halfplane},
    view::BBox{2};
    tol::Real = 1e-11,
)
    lines = [h.boundary for h in halfplanes]
    return vertical_decomposition(lines, view; tol = tol)
end

function vertical_decomposition(
    bundle::AbstractVector;
    factor::Real = 1.1,
    tol::Real = 1e-11,
)
    view = default_view(bundle; factor = factor, tol = tol)
    return vertical_decomposition(bundle, view; tol = tol)
end

# Constructors for ArrVertDecomp
ArrVertDecomp(bundle::AbstractVector, view::BBox{2}; tol::Real = 1e-11) =
    vertical_decomposition(bundle, view; tol = tol)

ArrVertDecomp(bundle::AbstractVector; factor::Real = 1.1, tol::Real = 1e-11) =
    vertical_decomposition(bundle; factor = factor, tol = tol)

###############################################
### Face reconstruction from Trapezoidal Map

# Fill `right_neighbor`/`left_neighbor` for trapezoids whose vertical sides touch and
# overlap, skipping x-coordinates that are input vertical lines.
function _tc_link_neighbors!(right_neighbor, left_neighbor, traps, vert_line_x, tol)
    n = length(traps)
    for i = 1:n, j = 1:n
        i == j && continue
        ti, tj = traps[i], traps[j]
        # ti must be to the left of tj
        abs(ti.x_right - tj.x_left) <= tol || continue
        # Make sure this x is not an input vertical line
        any(xv -> abs(xv - ti.x_right) <= tol, vert_line_x) && continue
        # Check vertical boundary overlap between ti's right edge [c2.y, c3.y]
        # and tj's left edge [c1.y, c4.y]
        y_bot = max(ti.corners[2][2], tj.corners[1][2])
        y_top = min(ti.corners[3][2], tj.corners[4][2])
        if y_top - y_bot > tol
            right_neighbor[i] = j
            left_neighbor[j] = i
        end
    end
    return nothing
end

"""
    trapezoid_chains(decomp::ArrVertDecomp; tol::Real = 1e-9)

Partition the vertical trapezoids in `decomp` into chains of trapezoids that share vertical
boundaries (which are internal cuts of the arrangement faces).
Returns a `Vector{Vector{Int}}` where each inner vector is a sequence of trapezoid indices
ordered from left to right forming a single face of the arrangement.
"""
function trapezoid_chains(decomp::ArrVertDecomp; tol::Real = 1e-9)
    traps = decomp.trapezoids
    n = length(traps)
    n == 0 && return Vector{Vector{Int}}()

    right_neighbor = fill(0, n)
    left_neighbor = fill(0, n)

    # Detect vertical lines in original input lines if any (to avoid connecting across them)
    vert_line_x = Float64[l.p[1] for l in decomp.lines if abs(l.u[1]) <= 1e-12]

    _tc_link_neighbors!(right_neighbor, left_neighbor, traps, vert_line_x, tol)

    chains = Vector{Vector{Int}}()
    visited = fill(false, n)

    for i = 1:n
        if left_neighbor[i] == 0
            chain = Int[i]
            visited[i] = true
            curr = i
            while right_neighbor[curr] != 0
                curr = right_neighbor[curr]
                push!(chain, curr)
                visited[curr] = true
            end
            push!(chains, chain)
        end
    end

    # Fallback for any unvisited trapezoids (e.g. isolates or cycles)
    for i = 1:n
        if !visited[i]
            push!(chains, Int[i])
            visited[i] = true
        end
    end

    return chains
end

"""
    faces(decomp::ArrVertDecomp; tol::Real = 1e-9)
    faces(bundle, view; tol::Real = 1e-9)
    faces(bundle; factor::Real = 1.1, tol::Real = 1e-9)

Compute the convex polygonal faces of the arrangement inside the view box.
Each face corresponds to the union of a chain of vertical trapezoids sharing vertical boundaries.
Returns a `Vector{PntSeq{2, Float64}}`.
"""
function faces(decomp::ArrVertDecomp; tol::Real = 1e-9)
    chains = trapezoid_chains(decomp; tol = tol)
    traps = decomp.trapezoids
    face_polys = PntSeq{2,Float64}[]

    for c in chains
        corners = Point{2,Float64}[]
        for trap_idx in c
            append!(corners, traps[trap_idx].corners)
        end
        poly = convex_hull(corners)
        push!(face_polys, poly)
    end
    return face_polys
end

faces(bundle::AbstractVector, view::BBox{2}; tol::Real = 1e-9) =
    faces(vertical_decomposition(bundle, view; tol = tol); tol = tol)

faces(bundle::AbstractVector; factor::Real = 1.1, tol::Real = 1e-9) =
    faces(vertical_decomposition(bundle; factor = factor, tol = tol); tol = tol)

export VerticalTrapezoid, VerticalTrapezoid2F
export ArrVertDecomp, ArrVertDecomp2F
export vertical_decomposition, trapezoids, view_box, polygon
export trapezoid_chains
