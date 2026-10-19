# Imports and Extras

## The default API

Starting with v0.5.0, `using BasicCompGeometry` brings 130 core names into scope:
geometric types, useful coordinate aliases, common operations, main algorithms,
and feature namespaces such as `BBT` and `IpeDraw`.

Another 91 supported names are available through qualified calls, explicit
imports, or the opt-in `Extras` submodule. These include generic helpers such as
`area` and `width`, specialized accessors, compatibility aliases, auxiliary
namespaces, and Cairo/LaTeX utilities. Their definitions, signatures, and behavior
are unchanged.

## Qualified calls and explicit imports

Use the alias `CG` when combining packages that export overlapping names:

```@example imports
import BasicCompGeometry as CG
bb = CG.BBox(CG.point(0.0, 0.0), CG.point(2.0, 3.0))
CG.width(bb)
```

An explicit import works for any supported helper, even if it is not exported:

```@example imports
using BasicCompGeometry: width
width(bb)
```

Feature submodules have their own exports. For example, `CG.BBT.depth` remains
exported by `BBT`, while `CG.depth` is available through qualification or `Extras`.
Qualify core types such as `CG.Point` when another package exports the same name.

## Loading both sets of exports

```@example extras
import BasicCompGeometry as CG
using .CG, .CG.Extras

bb = BBox(point(0.0, 0.0), point(2.0, 3.0))
width(bb)
```

The equivalent without an alias is:

```julia
using BasicCompGeometry
using BasicCompGeometry.Extras
```

`Extras` exports an explicit list of supported public names and shares their
original bindings. It adds no wrappers and does not export internal
implementation utilities. Julia permits access to unexported bindings, but
that does not make implementation helpers supported API.

Loading `Extras` does not activate optional extensions. Load `Cairo`, `Colors`,
or `LaTeXStrings` separately as described in [Cairo and LaTeX](cairo.md).

## Migrating existing code

If a previously unqualified helper is undefined after updating, qualify it with
`CG`, import it explicitly, or opt into `Extras`. See the
[complete v0.5.0 migration list](https://github.com/sarielhp/BasicCompGeometry.jl/blob/main/CHANGELOG.md#050---unreleased)
for all 91 moved names. The split changes imports rather than implementations.

```@docs
BasicCompGeometry.Extras
```
