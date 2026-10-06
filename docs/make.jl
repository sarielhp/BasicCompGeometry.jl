using Documenter
using BasicCompGeometry

makedocs(
    sitename = "BasicCompGeometry.jl",
    modules = [BasicCompGeometry],
    pages = ["Home" => "index.md"],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link = "main",
    ),
    warnonly = [:missing_docs, :cross_references, :autodocs_block],
)

deploydocs(repo = "github.com/sarielhp/BasicCompGeometry.jl.git", devbranch = "main")
