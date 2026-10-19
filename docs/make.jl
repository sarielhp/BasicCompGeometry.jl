using Documenter
using BasicCompGeometry

# Extras reexports supported bindings whose parent names are intentionally
# unexported. Filter by the public objects so their docstrings appear too.
const api_modules = [
    BasicCompGeometry, BasicCompGeometry.Extras, BasicCompGeometry.BBT,
    BasicCompGeometry.WSPD, BasicCompGeometry.MVBB, BasicCompGeometry.MetricSpace,
    BasicCompGeometry.LongestConvexSubset, BasicCompGeometry.ConvexHull3D,
    BasicCompGeometry.IpeDraw, BasicCompGeometry.ReadWrite, BasicCompGeometry.VirtArray,
    BasicCompGeometry.PntSeqHierarchy,
]
const documented_api = Set(
    getfield(mod, name) for mod in api_modules for name in names(mod)
)
const ipe_composition = Set(
    getfield(BasicCompGeometry.IpeDraw, name) for name in (
        :IpeDraw, :IpeCanvas, :Viewport, :PageBox, :Style, :Theme, :publication_theme,
        :figure, :open_ipe, :edit_ipe, :add_preamble!, :draw!, :mark!, :label!,
        :layer, :with_style, :fit!, :page_space, :inset, :clip_to, :legend!, :scale_bar!,
    )
)

makedocs(
    sitename = "BasicCompGeometry.jl",
    modules = [BasicCompGeometry],
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting-started.md",
        "Guides" => [
            "Imports and Extras" => "imports.md",
            "Geometry and point sequences" => "geometry.md",
            "Kernels and numerical scope" => "kernels.md",
            "Algorithms and spatial structures" => "algorithms.md",
            "Finite metric spaces" => "metrics.md",
            "Vector figures" => "figures.md",
            "Cairo and LaTeX" => "cairo.md",
        ],
        "API reference" => [
            "Primitives and transformations" => "reference/primitives.md",
            "Point sequences" => "reference/sequences.md",
            "Circles and ellipses" => "reference/round-shapes.md",
            "Curves and conics" => "reference/curves.md",
            "Kernels and predicates" => "reference/predicates.md",
            "Halfplanes and arrangements" => "reference/arrangements.md",
            "Hull and diameter algorithms" => "reference/algorithms.md",
            "Spatial structures" => "reference/spatial.md",
            "Metric spaces" => "reference/metrics.md",
            "Simplification hierarchies" => "reference/hierarchy.md",
            "Ipe canvas and composition" => "reference/ipe.md",
            "Ipe drawing and export" => "reference/ipe-drawing.md",
            "Cairo and LaTeX" => "reference/cairo.md",
            "Data files" => "reference/io.md",
        ],
        "About" => "about.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link = "main",
    ),
)

deploydocs(repo = "github.com/sarielhp/BasicCompGeometry.jl.git", devbranch = "main")
