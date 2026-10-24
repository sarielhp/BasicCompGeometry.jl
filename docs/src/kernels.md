# Kernels and numerical scope

## Choosing predicate arithmetic

Ordinary calls use floating-point predicates. Kernel-aware operations accept
`kernel=` so that exact decisions can be requested locally:

```@example kernels
using BasicCompGeometry

points = PntSeq([point(0, 0), point(4, 0), point(4, 3), point(0, 3), point(2, 1)])
hull = convex_hull(points; kernel=ExactKernel())
is_left_turn(hull[1], hull[2], hull[3]; kernel=ExactKernel())
```

The same kernel passes through every predicate in a kernel-aware algorithm.
No global setting, context object, or new point type is required.

| Kernel | Arithmetic |
|:-------|:-----------|
| `FloatingKernel()` | `Float64` predicate arithmetic |
| `FloatingKernel(Float32)` | `Float32` predicate arithmetic |
| `ExactKernel()` | Arbitrary-precision integers and rationals; algebraic comparison of circular radicals |

Finite floating-point inputs to `ExactKernel()` mean their stored binary values,
not the decimal values a caller may have intended. Integer predicates also use
floating-point arithmetic by default; request the exact kernel when integer
exactness matters.

## Supported operations

Kernel-aware operations currently include:

- 2D turn predicates (`turn_sign`, `orient`, `orient_sign`) and lexicographic point comparison (`compare_xy`).
- Closed segment membership (`on_segment`) and proper interior segment crossing (`proper_crossing`).
- Counterclockwise angular order of direction vectors (`compare_angle`, `angle_less`).
- Canonical undirected segment keys (`segment_key`).
- Exact point constructions: `midpoint`, `affine_interpolate`, and nonvertical segment evaluation at `x` (`segment_at_x`).
- Polygon signed double area (`twice_area`) and linear-time monotone polygon triangulation (`triangulate_monotone`).
- Distance and squared-distance threshold comparisons (`compare_distance`, `compare_squared_distance`).
- Sphere containment and right-circle-boundary comparisons (`compare_right_arcs`, `compare_right_arc_to_point`).
- The 2D convex hull (`convex_hull`).

Other algorithms and coordinate constructions have not been converted. An exact
predicate cannot repair an already rounded intersection. `BigFloat` alone is
not an exact kernel, and `FloatingKernel` is not a certified filter.

Most algorithms use ordinary floating-point decisions. Degenerate or
ill-conditioned inputs can produce incorrect combinatorial results, especially
in hull and arrangement code. The package targets research code, experiments,
and figure generation; select arithmetic and algorithms appropriate to the
inputs. The WSPD implementation is tested for separation and exact pair coverage,
which is distinct from exact arithmetic.

## Custom kernels

Custom kernels subtype `AbstractKernel` and implement kernel-first predicate
methods, such as `turn_sign(kernel, p, q, r)` and
`compare_xy(kernel, p, q)`. Missing methods fail rather than silently falling
back to another arithmetic policy. Algorithms must forward the chosen kernel
to every supported predicate.

The interface can accommodate a future certified filtered kernel; this release
does not provide one. See the [kernel reference](reference/predicates.md) and
the [turn predicates](reference/primitives.md).
