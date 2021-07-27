using Documenter, OrdinalMultinomialModels

makedocs(
    format = Documenter.HTML(),
    sitename = "OrdinalMultinomialModels.jl",
    clean = true,
    debug = true,
    pages = [
        "index.md"
    ]
)

deploydocs(
    repo   = "github.com/OpenMendel/OrdinalMultinomialModels.jl.git",
    target = "build",
    deps = nothing,
    make = nothing
)

