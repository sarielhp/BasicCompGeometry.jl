# Shared controls for reproducible, bounded example runs. Normal runs keep their defaults.
const EXAMPLE_SMOKE = get(ENV, "BCG_EXAMPLE_SMOKE", "false") == "true"
example_size(normal, small) = EXAMPLE_SMOKE ? small : normal
example_output(parts...) = joinpath(
    get(ENV, "BCG_EXAMPLE_OUTPUT_DIR", normpath(joinpath(@__DIR__, "..", "..", "output"))),
    parts...,
)
if EXAMPLE_SMOKE
    using Random
    Random.seed!(42)
end
