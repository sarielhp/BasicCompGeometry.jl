# Archived code

Historical source removed from the active package, and finished plans with their
experiment scripts, are kept here for reference.
These files are not included by the package, documented as supported APIs, or
run by the test suite. Review and test archived code before reusing it.

| File | Original location | Last active revision | Reason |
| --- | --- | --- | --- |
| [MetricSpace_vanity.jl](MetricSpace_vanity.jl) | `src/MetricSpace.jl` | `f7063bef584d6dfa46ca41cbf58ba55b021814e0` | Unused vanity permutation and its private helper; removed before v0.5.0. |
| [plans/improve_u_kd.md](plans/improve_u_kd.md), [plans/improve_u_kd/](plans/improve_u_kd/) | `plans/` | `75b63d2` | Finished plan and benchmark/fuzz scripts for the u_kd ear-clipping speedups. Run scripts from their directory with `--project=../../..`. |
| [plans/diagonal_split.md](plans/diagonal_split.md) | `plans/` | `75b63d2` | Abandoned diagonal-split triangulation; code recoverable from `8939c08`. |
