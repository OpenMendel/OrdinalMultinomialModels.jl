using Documenter, OrdinalMultinomialModels

ENV["DOCUMENTER_DEBUG"] = "true"

makedocs(
    format = :html,
    sitename = "OrdinalMultinomialModels",
    modules = [PolrModels]
)

deploydocs(
    repo   = "github.com/OpenMendel/OrdinalMultinomialModels.jl.git",
    target = "build"
)
