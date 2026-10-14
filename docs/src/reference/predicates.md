# Kernels and predicates

See [Kernels and numerical scope](../kernels.md) for arithmetic policies and
their limits. Turn predicates such as `turn_sign` and `is_left_turn` are in the
[primitive reference](primitives.md). Circle-boundary comparisons belong to
[Extras](../imports.md).

```@autodocs
Modules = [BasicCompGeometry]
Pages = ["Kernels.jl", "KernelPredicates.jl"]
Filter = obj -> obj in Main.documented_api
```
