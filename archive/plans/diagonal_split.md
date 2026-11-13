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
  comb (`archive/plans/improve_u_kd/polys.jl`).
- Correctness: `archive/plans/improve_u_kd/fuzz_exact.jl` / `fuzz_float.jl` with the new algorithm.

## Open questions / later

- After `try(a, m)` succeeds, continue with `try(m, b)` in the remainder instead of
  re-picking a random pair.
- Start the edge query near `a` (climb) instead of from the root.
- Avoid copying the larger side on unbalanced splits.

## Status 2026-11-13

### Iteration log

1. **Pair + halving search** (commit 5f48f1a). Correct, but slow: star 1157 ms vs u_kd 56
   (2^16). The spiral is quadratic (5169 tests / n): its reflex vertices are all on the inner
   wall and never see each other, so the halving search burns O(k) failed (cone) tests per
   split.
2. **Gallop on convex chains + conservative Float64 predicates** (9076f28). Galloping did
   not help the spiral (the waste is in ranges that contain reflex vertices). The Float64
   fuzz found one silently wrong result (a vertex 1e-16 from a candidate diagonal); all
   orientation decisions now treat values within Shewchuk's static bound as collinear,
   which is always the conservative side.
3. **Ray resolution** (this commit). After a failed pair test, shoot a ray from the reflex
   endpoint along its interior bisector; the first hit (original edge via the edge tree, or
   one of the piece's own diagonals via a linear scan) gives a triangle (a, h, x) whose
   angularly first reflex vertex (or x) is visible from `a` — found by a branch-and-bound
   query on a static Kd-tree of the input's reflex vertices. Validated by `_diag_ok`.
   Every ray resolution so far has produced a valid diagonal; the halving search and the
   ear-clip fallback are never reached on the benchmark inputs.

   Storing diagonals in the edge tree was tried and dropped: long diagonals end up in the
   root's list and every query scans them (comb: 3900 segment tests per query).

### Measurements (2^16, ms; quiet-ish machine)

| input | sweep | u_kd | diagonal_split | tests/n | work/n |
|---|---|---|---|---|---|
| star | 29 | 54 | 1199 | 0.71 | 327 |
| altstar | 17 | 7.7 | 24 | 0.51 | 16 |
| coast | 20 | 24 | 104 | 0.59 | 94 |
| spiral | 7.6 | 11 | 31 | 0.55 | 20 |
| comb | 15 | 7.0 | 40 | 0.65 | 43 |

`work/n` = total size of the pieces that were split, per vertex (balanced: ~log2 n = 16).

### Remaining problems

- **Unbalanced splits** (coast, star): ray resolution often cuts off a small bay, and the
  large remainder is copied and rescanned (piece diagonals, `findfirst`) each time:
  ~n^1.5 work on coast. Fix: linked-list pieces with smaller-half relabeling (O(n log n)),
  or prefer the larger of the two candidate splits.
- **Edge boxes on the star**: its edges are long radial spikes, so a chord through the
  center crosses ~20% of all edge boxes (nodes / test grows almost linearly in n). Box
  pruning cannot fix this; it needs edges clipped into cells (each leaf lists the edges
  crossing its cell) so a query walks only the cells along the segment.
- Even on balanced inputs (altstar, spiral) it is ~3x slower than u_kd per vertex, despite
  0.5 tests / n vs ~1.1-2.5 ear tests / n: a segment query (~50 nodes + exact edge tests)
  plus the edge-tree build (~16% of the time) cost more than the saved ear tests.

## Outcome 2026-11-13: abandoned

Dropped. The diagonal split saves ear tests, but each segment query costs more than the
tests it saves, unbalanced splits cost ~n^1.5 on coast, and the star needs a different edge
index altogether. `:ear_cut_u_kd_tree` stays as the ear-clipping algorithm. The code
(`src/DiagonalSplitTriangulation.jl`, `plans/diagonal_split/bench.jl`) was removed; it is
recoverable from commit 8939c08.
