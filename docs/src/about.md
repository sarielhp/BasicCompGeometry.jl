# About the package

```@docs
BasicCompGeometry
```

## Origins

Sariel Har-Peled originally developed this library as part of `FrechetDist`,
with substantial portions based on computational geometry algorithms and C++
code written over 25 years ago. BasicCompGeometry does not implement Fréchet
distance.

The package has since been reorganized around a flat module hierarchy,
type-generic implementations, `StaticArrays.jl`, and methods for standard
`Base` operations through multiple dispatch.

## AI contributions and code ownership

Recent modernizations, extensions, tests, documentation, and predicate kernels
were written and generated with assistance from AI coding tools, including
Google Gemini CLI and OpenAI Codex.

As author and maintainer, Sariel Har-Peled understands the codebase, curates its
mathematical and algorithmic design, and takes responsibility for all code and
its maintenance, bug fixes, and development.

## Maintenance and license

Maintenance is best-effort, with correctness bugs prioritized. Bug reports with
a small reproducing example are welcome through
[GitHub issues](https://github.com/sarielhp/BasicCompGeometry.jl/issues).

The package is distributed under the
[MIT License](https://github.com/sarielhp/BasicCompGeometry.jl/blob/main/LICENSE).
