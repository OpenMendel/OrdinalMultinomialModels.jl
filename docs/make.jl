using Documenter, PolrModels

ENV["DOCUMENTER_DEBUG"] = "true"

makedocs(
    format = :html,
    sitename = "PolrModels",
    modules = [PolrModels]
)

deploydocs(
    repo   = "github.com/OpenMendel/PolrModels.jl.git",
    target = "build"
)
