# Plan: diagonal-split triangulation (`:diagonal_split`)

Started 2026-11-13. Prototype in `src/DiagonalSplitTriangulation.jl`.

## Idea

Split the polygon recursively by diagonals between reflex vertices. A diagonal between two
reflex vertices resolves (makes convex) each endpoint in at least one of the two pieces, so
the reflex count never increases. Pieces with no reflex vertex are fanned; pieces with one
reflex vertex are fanned from it (a polygon with a single reflex vertex is star-shaped from
it). No queries are needed for these.

For a piece with `k >= 2` reflex vertices (collinear vertices count as reflex):

1. Pick a random reflex vertex `a` and the reflex vertex `b` that is `k/2` positions after it
   in boundary order.
2. `try(a, b)`: if `ab` is a valid diagonal, split there. Otherwise pick a middle vertex `m`
   of the chain `a -> b`: the middle reflex vertex strictly between them if one exists, else
   the middle vertex of the (convex) chain. Recurse on `try(a, m)`, then `try(m, b)`. Stop
   when the chain is a single edge.
3. If `try(a, b)` and `try(b, a)` both fail, fall back to u_kd ear clipping on the piece.

## Diagonal test

`ab` (both vertices of the current piece, non-adjacent) is a diagonal of the piece iff

- `ab` lies in the open interior cone at `a` and at `b`, with respect to the *piece*
  neighbors, and
- no *original* polygon edge meets the closed segment `ab` (edges incident to `a` or `b` only
  matter if their other endpoint lies on `ab`).

Earlier diagonals need not be indexed: leaving the piece through a diagonal enters a
component of `P \ piece` that is attached to the piece only through that diagonal, so
getting back to `b` would cross an original edge.

The reflex u_kd tree cannot answer this (a blocking notch can cross `ab` with its reflex
tip anywhere beyond the segment), so the edges get their own static Kd-tree: built over edge
midpoints, node boxes refit to the edges' extents (outward rounded to Float64), pruned by box
overlap and by a conservative line-side test. Leaf tests are exact under the kernel.

## Pieces

A piece's boundary is a subsequence of the original boundary order: `verts` (vertex ids)
plus `refl` (sorted positions of its reflex vertices). Splitting at positions `i < j` copies
`verts[i:j]` and `verts[j:end] ∪ verts[1:i]`; the reflex lists are slices of `refl`, with
the two endpoints re-classified locally. Copying is `O(piece size)` per split (balanced
splits: `O(n log n)` total; unbalanced splits can be worse — revisit if it shows up).

## Measurements to make

- Diagonal tests / n, failed tests / n, fallback pieces, edge-tree nodes / test.
- Wall clock vs `:sweep` and `:ear_cut_u_kd_tree` on `star`, `altstar`, `coast`, `spiral`,
  comb (`plans/improve_u_kd/polys.jl`).
- Correctness: `plans/improve_u_kd/fuzz_exact.jl` / `fuzz_float.jl` with the new algorithm.

## Open questions / later

- After `try(a, m)` succeeds, continue with `try(m, b)` in the remainder instead of
  re-picking a random pair.
- Start the edge query near `a` (climb) instead of from the root.
- Avoid copying the larger side on unbalanced splits.
