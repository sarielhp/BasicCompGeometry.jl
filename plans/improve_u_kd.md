# Plan: improve `triangulate_ear_cut_u_kd_tree`

Status as of 2026-11-13. Two improvements are implemented and committed (see
"Implemented"), validated for correctness, and benchmarked only on an overloaded machine
(load ~25 on 16 cores), so timings must be re-measured on a quiet machine before committing
any performance claims.

## Where the time goes (diagnosis)

Instrumented counts (`n = 2^18`, u_kd, original code):

| input | ear tests / n | DFS nodes / test | leaf points / test |
|---|---|---|---|
| random star (`star`) | 1.45 | 83 | 185 |
| alternating star (`altstar`) | 2.0 | 3 | 4.5 |
| coastline (`coast`) | 1.65 | 19 | 28 |
| thick spiral (`spiral`) | 2.5 | 9 | 4.3 |

- ~85% of time is in `_dfs_ear_blocked!` (profile on `star`).
- Virtually every DFS call *proves emptiness* (blocked ears are caught earlier by
  `blocker_hint`), so DFS cost is the cost of empty triangle range queries.
- On the random star, nodes/test grows like `sqrt(n)`, and the most expensive queries are
  the *smallest* ears (span 2: three consecutive input vertices, 41% of node visits, ~94
  nodes and ~239 leaf points each). Those ears are long radial slivers (radius 0.5 -> 1.5),
  whose bounding box covers a large part of the polygon, so the bbox prefilter passes most
  points and the Kd-tree query is line-like (`Theta(sqrt r)` nodes). This is inherent to the
  input for any Kd-tree; sweep (`O(n log n)`) remains faster there.

## Tried and rejected

- **Priority queue by ear size** (smallest bbox side / area / perimeter / new-diagonal length
  first instead of FIFO): slower on every input (heap overhead, no reduction in query cost).
- **Leaf size** 4 / 8 / 16 / 32 for the reflex Kd-tree: differences within noise; keep 16.
- **Fully lazy queue seeded with all vertices**: fewer tests but builds fans of long
  triangles from one corner on the comb (`_hcomb`), blowing up Z-order orient counts
  (104/n vs < 10) and failing the comb performance tests.
- **Lazy queue seeded with initial ears, queued candidates keep their old queue position**:
  still fans on the comb for Z-order (32/n).

## Implemented (committed as "Speed up ear clipping: lazy re-testing queue and leaf pre-scan")

1. **Branch-free leaf pre-scan** (`src/EarCutKdTreeTriangulation.jl`,
   `_leaf_may_block`): for the Float64 query, each leaf's live points are first scanned with
   a SIMD-friendly loop that only ORs "inside closed triangle and not equal to `pr`/`nx` by
   coordinates" (valid because input vertices are now pairwise distinct and the ear tip is
   convex, hence never live). Only on a hit does the exact scalar scan run (to find the
   blocker for `blocker_hint`). Same predicate formula as `_pt_in_tri`, so a leaf is skipped
   only when the scalar scan would find nothing. In-process A/B: star 0.77x, coast 0.83x,
   altstar/spiral neutral.
2. **Ear queue: lazy re-testing with re-queue-to-back** (`src/MonotoneTriangulation.jl`,
   `_ear_clip_run!`; shared by all ear cutters): the queue is seeded with the initial ears
   (verified eagerly, as before). When an ear is clipped, each neighbor is appended to the
   back as an unverified candidate, superseding any earlier entry (`entry_ver` / `pushed_ver`
   stamps), and is tested only when that entry is dequeued. Initial ears whose neighbors
   never changed are not re-tested (`verified` stamp). In-process A/B vs the old loop (same
   index code):

   | index | star | altstar | coast | spiral | comb |
   |---|---|---|---|---|---|
   | u_kd | 0.83 | 0.82 | 0.79 | 0.57 | 0.98 |
   | kd | 0.82 | 0.93 | 0.78 | 0.57 | 1.00 |
   | z-order | – | 0.86 | 0.74 | 0.13 | 0.89 |

   This changes the ear order, hence the output diagonals of `:ear_cut` and all indexed ear
   cutters (they still agree with each other exactly).

Validation done: `test/test_monotone_triangulation.jl` passes (81,697 tests, including the
comb performance thresholds); exact-arithmetic fuzz (`fuzz_exact.jl`, seeds 11, 12) clean for
all six algorithms. Pending at time of writing: full `Pkg.test()`, exact fuzz seeds 21/22, and
the Float64 fuzz for `ear_cut`, `ear_cut_u_kd_tree`, `ear_cut_z_order`.

## Re-measured on a quiet machine (2026-11-13): confirmed, kept

Correctness: full `Pkg.test()` passes (156,337 tests); exact fuzz seeds 11, 12, 21, 22, 31
(400 polygons each) clean for all six algorithms; Float64 fuzz (4,500 cases each) for
`ear_cut`, `ear_cut_u_kd_tree`, `ear_cut_z_order`: 0 wrong results, 9 `ArgumentError`s each
(previously 6 for `ear_cut`).

New vs old queue (`ab.jl`, new/old time):

| index | star | altstar | coast | spiral | comb |
|---|---|---|---|---|---|
| u_kd | 0.81 | 1.04 | 0.77 | 0.54 | 0.99 |
| kd | 0.80 | 0.88 | 0.80 | 0.61 | 1.02 |
| z-order | – | 0.83 | 0.72 | 0.13 | 0.91 |

Leaf pre-scan (`ab2.jl`, with/without): star 0.70, altstar 1.12, coast 0.78, spiral 1.03.
It costs ~10% on altstar, where a query touches only ~4 points per leaf.

Comb: ~2.4 ear tests per vertex for all three ear cutters, no fallbacks.

Absolute times (ms) vs `:sweep`:

| input | 2^16 sweep | 2^16 u_kd | 2^18 sweep | 2^18 u_kd |
|---|---|---|---|---|
| star | 29.2 | 54.2 | 138 | 356 |
| altstar | 16.9 | 8.0 | 70 | 41 |
| coast | 21.9 | 24.8 | 104 | 135 |
| spiral | 7.9 | 10.9 | 41 | 55 |

u_kd beats sweep only on altstar; the star stays ~2.5x slower (line-like queries on radial
slivers). Accepted as good enough. Follow-ups tried and dropped: a sliding window of recently
visited nodes (high reuse but little work saved) and a diagonal-split algorithm
(`plans/diagonal_split.md`). A u_kd → sweep hand-off on a work budget was discussed but not
built. CHANGELOG entry added under [Unreleased].

## Next steps (on a quiet machine) — done, see above

1. Re-run correctness:
   - `julia --project=. -e 'using Pkg; Pkg.test()'`
   - `cd plans/improve_u_kd && julia --project=../.. fuzz_exact.jl <seed> 400` for a few seeds
     (expect `fails=Dict()`).
   - `julia --project=../.. fuzz_float.jl 5 1500 <alg>` for `ear_cut`, `ear_cut_u_kd_tree`,
     `ear_cut_z_order`. Before these changes: `ear_cut` 0 wrong / 6 `ArgumentError` of ~4500.
     An earlier fully-lazy variant produced 1 silently wrong result; confirm the final
     variant gives 0 (any "FAIL ... non-ccw tri" line is a silently wrong triangulation).
2. Re-measure performance:
   - `julia --project=../.. ab.jl`: old vs new queue in one process (`old_run.jl` is the
     pre-change `_ear_clip_run!`, loaded as `_ear_clip_run_old!`).
   - `julia --project=../.. ab2.jl`: with vs without the leaf pre-scan.
   - Absolute timings vs `:sweep` on the inputs in `polys.jl` (`star`, `altstar`, `coast`,
     `spiral`) at `n = 2^16, 2^18`, plus the comb (`comb.jl` also prints ear-test counts).
3. If confirmed: add a CHANGELOG entry under
   [Unreleased] (u_kd / kd / z-order speedups; ear order and therefore `:ear_cut` output
   diagonals changed), and update any doc text that describes the ear queue. If not, revert the commit (or just the part that does not pay off).
4. Further ideas, not yet tried:
   - **Cheaper child tests in the DFS**: `active_count` lives in a separate array (extra
     cache miss per child); pack it next to the node or test it after the bbox check.
   - **Dynamic box shrinking**: node boxes are built over the *initial* reflex set and never
     shrink as reflex vertices deactivate; refitting boxes of nodes whose active count drops
     sharply could prune more on long-running inputs.
   - **Start node for the climb**: for slivers, the unshrunk-cell climb often goes high;
     starting from the leaf of the ear's nearest reflex vertex, or descending from the root
     with the triangle's bbox, might visit fewer nodes.
   - **Sliver-specific test**: for very thin ears (bbox much larger than the triangle area
     suggests), a segment/stabbing-style query along the long edge may beat the box DFS.

## Scripts (`plans/improve_u_kd/`)

- `polys.jl`: input generators (`star`, `altstar`, `coast`, `spiral`).
- `ab.jl`, `old_run.jl`: in-process A/B of old vs new ear queue for u_kd / kd / z-order.
- `ab2.jl`: in-process A/B of the leaf pre-scan.
- `comb.jl`: ear-test counters on the comb for `ear_cut`, u_kd, z-order.
- `fuzz_exact.jl`, `lib.jl`: exact-arithmetic random simple polygons, all six algorithms.
- `fuzz_float.jl`: Float64 polygons checked exactly; `ArgumentError` counted separately from
  wrong results. Writes its last case to `output/tmp/`.
