module PolrModels

using Distributions, Reexport, StatsBase
@reexport using GLM
@reexport using Ipopt
@reexport using NLopt
using MathProgBase
import Base.LinAlg: BlasReal

export
    # types
    PolrModel,
    PolrScoreTest,
    # functions
    coeftable,
    deviance,
    dof_residual,
    loglikelihood,
    nobs,
    polrtest,
    polrfun!,
    polrmle,
    polrtest,
    predict,
    rpolr,
    stderr,
    vcov

struct PolrModel{TY<:Integer, T<:BlasReal, TL<:GLM.Link} <: MathProgBase.AbstractNLPEvaluator
    # dims
    n::Int
    p::Int
    J::Int
    npar::Int
    # data
    Y::Vector{TY}
    X::Matrix{T}
    # parameters
    θ::Vector{T}
    α::Vector{T}
    β::Vector{T}
    link::TL
    # working parameters
    η::Vector{T}
    dθdα::Matrix{T}
    "Jacobian: dθdα[i, j] = dθj/dαi"
    ∇::Vector{T}
    "gradient wrt (θ, β)"
    H::Matrix{T}
    "Hessian wrt (θ, β)"
    vcov::Matrix{T}
    "vcov = inv(-H)"
    wt∂β::Vector{T}
    "weight vector, n x 1"
    wt∂θ∂β::Matrix{T}
    "weight matrix, n x (J - 1)"
    wt∂β∂β::Vector{T}
    "weight vector, n x 1"
    scratchm1::Matrix{T}
    "scratch matrix of same size as X"
end

# Constructor
function PolrModel(
    y::Vector{TY},
    X::Matrix{T},
    link::GLM.Link) where TY <: Integer where T <: BlasReal
    # check y has observations in each category
    yct = counts(y)
    J   = length(yct)
    J < 2 && throw(ArgumentError("Response must have 3 or more levels"))
    for j in 1:J
        yct[j] < 1 && throw(ArgumentError("No observations in category $j"))
    end
    n, p   = size(X)
    npar   = J - 1 + p
    η      = zeros(T, n)
    θ      = zeros(T, J - 1)
    α      = zeros(T, J - 1)
    β      = zeros(T, p)
    dθdα   = zeros(T, J - 1, J - 1)
    ∇      = zeros(T, J - 1 + p)
    H      = zeros(T, J - 1 + p, J - 1 + p)
    vcov   = zeros(T, J - 1 + p, J - 1 + p)
    wt∂β   = zeros(T, n)
    wt∂θ∂β = zeros(T, n, J - 1)
    wt∂β∂β = zeros(T, n)
    scratchm1 = zero(X)
    PolrModel{eltype(y), eltype(X), typeof(link)}(n, p, J, npar, y, X, θ, α, β,
        link, η, dθdα, ∇, H, vcov, wt∂β, wt∂θ∂β, wt∂β∂β, scratchm1)
end

coef(m::PolrModel) = [m.θ; m.β]
deviance(m::PolrModel) = -2polyrfun!(m, false, false)
dof_residual(m::PolrModel) = m.n - m.npar
fitted(m::PolrModel) = nothing # TODO
loglikelihood(m::PolrModel) = polyrfun!(m, false, false)
nobs(m::PolrModel) = m.n
predict(m::PolrModel) = nothing # TODO
stderr(m::PolrModel) = sqrt.(diag(m.vcov))
vcov(m::PolrModel) = m.vcov

function coeftable(m::PolrModel)
    cc = coef(m)
    se = stderr(m)
    tt = cc ./ se
    CoefTable(hcat(cc,se,tt,ccdf.(FDist(1, dof_residual(m)), abs2.(tt))),
              ["Estimate","Std.Error","t value", "Pr(>|t|)"],
              ["x$i" for i = 1:m.npar], 4)
end

include("polrrand.jl")
include("polrfit.jl")
include("polrtest.jl")

end # module
