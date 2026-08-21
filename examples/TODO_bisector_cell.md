TODO: Fix examples/bisector_cell.jl

The bisector of a line and a point is a parabola, not a hyperbola.
The example currently uses `Hyperbola(line, origin)` which constructs
a hyperbola approximation. It should use `Parabola(origin, line)` instead.

Changes needed:
1. Replace `Hyperbola(line, origin)` with `Parabola(origin, line)` in `cell_containing_origin`
2. Update drawing code in `draw_bisector_cell`:
   - Replace `at(h, t; branch=...)` with `at(h, t)` (parabola has no branches)
   - Replace `b_coeff`, `rotation_angle` with `p_coeff`, `axis_direction`
   - The parameter range and sampling logic needs to change
3. The `angle_bisectors` function and vertex computation logic may also need adjustment

CRITICAL: `Hyperbola(line, point)` constructor in `src/Hyperbolas.jl` is fundamentally
wrong — the bisector of a line and a point is a parabola, not a hyperbola, and the
current implementation just makes up two foci that don't satisfy the distance
difference property. This constructor should either be removed entirely or renamed
to something like `Hyperbola(line, point; tol=...)` with a clear warning that it's
an approximation. The `bisector(line, point)` function that calls it is also wrong.