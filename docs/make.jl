using Documenter
using BasicCompGeometry

# Ensure we include the module path
push!(LOAD_PATH, "../src/")

makedocs(
    sitename = "BasicCompGeometry.jl",
    modules = [BasicCompGeometry],
    pages = ["Home" => "index.md"],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    warnonly = [:missing_docs, :cross_references, :autodocs_block],
)

deploydocs(repo = "github.com/sarielhp/BasicCompGeometry.jl.git", devbranch = "main")
