__precompile__()

module PolrModels

using Distributions, Reexport, StatsBase
@reexport using GLM
@reexport using Ipopt
@reexport using NLopt
using MathProgBase
import Base.LinAlg: BlasReal
import StatsBase: coef, coeftable, confint, deviance, nulldeviance, dof, dof_residual,
                    loglikelihood, nullloglikelihood, nobs, stderr, vcov, residuals,
                    predict, fit, model_response, r2, r², adjr2, adjr², PValue

export
    # types
    AbstractPolrModel,
    PolrModel,
    PolrScoreTest,
    # functions
    coef,
    coeftable,
    confint,
    cor,
    deviance,
    dof_residual,
    loglikelihood,
    nobs,
    polr,
    polrtest,
    polrfun!,
    polrfit,
    polrtest,
    predict,
    rpolr,
    stderr,
    vcov

abstract type AbstractPolrModel <: RegressionModel end

"""
    PolrModel

The data, parameters, and various derived variables for the proportional oddss
logistic regression model.
"""
struct PolrModel{TY<:Integer, T<:BlasReal, TL<:GLM.Link} <: MathProgBase.AbstractNLPEvaluator
    # dimensions
    "`n`: number of observations"
    n::Int
    "`n`: number of covariates, excluding intercept"
    p::Int
    "`J`: number of categories"
    J::Int
    "`npar`: number of parameters"
    npar::Int
    # data
    "`Y`: response vector"
    Y::Vector{TY}
    "`X`: covariate matrix"
    X::Matrix{T}
    "`wts`: prior observation weights, can be empty"
    wts::Vector{T}
    # parameters
    "`θ`: intecept parameters, satisfying `θ[1]≤...≤θ[J-1]`"
    θ::Vector{T}
    "`α`: unconstrained parameterization of `θ`"
    α::Vector{T}
    "`β`: regression coefficients"
    β::Vector{T}
    "`link`: link function"
    link::TL
    # working parameters
    "`η`: linear systematic component `Xβ`"
    η::Vector{T}
    "`dθdα`: Jacobian dθdα[i, j] = dθj/dαi"
    dθdα::Matrix{T}
    "`∇`: gradient wrt (θ, β)"
    ∇::Vector{T}
    "`H`: Hessian wrt (θ, β)"
    H::Matrix{T}
    "`vcov`:  `inv(-H)`"
    vcov::Matrix{T}
    "`wtwk`: working weights"
    wtwk::Vector{T}
    "`wt∂β``: weight vector for computing derivative, `n x 1`"
    wt∂β::Vector{T}
    "`wt∂θ∂β`: weight matrix for computing Hessian, n x (J - 1)"
    wt∂θ∂β::Matrix{T}
    "`wt∂β∂β`: weight vector for computing Hessian, n x 1"
    wt∂β∂β::Vector{T}
    "`scratchm1`: scratch matrix of same size as `X`"
    scratchm1::Matrix{T}
end

# Constructor
function PolrModel(
    X::Matrix{T},
    y::Vector{TY},
    wts::Vector{T} = similar(X, 0),
    link::GLM.Link = LogitLink()) where TY <: Integer where T <: BlasReal
    # check y has observations in each category
    yct = counts(y)
    J   = length(yct)
    J < 2 && throw(ArgumentError("Response must have 3 or more levels"))
    for j in 1:J
        yct[j] < 1 && throw(ArgumentError("No observations in category $j"))
    end
    if !isempty(wts)
        lw, ly = length(wts), length(y)
        lw ≠ ly && throw(ArgumentError("wts has length $lw, should be 0 or $ly"))
        minw = minimum(wts)
        minw < 0 && throw(ArgumentError("wts should be nonnegative, found entry "))
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
    wtwk   = zeros(T, n)
    wt∂β   = zeros(T, n)
    wt∂θ∂β = zeros(T, n, J - 1)
    wt∂β∂β = zeros(T, n)
    scratchm1 = zero(X)
    PolrModel{eltype(y), eltype(X), typeof(link)}(n, p, J, npar, y, X, wts,
        θ, α, β, link, η, dθdα, ∇, H, vcov, wtwk, wt∂β, wt∂θ∂β, wt∂β∂β, scratchm1)
end
PolrModel(X, y, link) = PolrModel(X, y, similar(X, 0), link)

coef(m::PolrModel) = [m.θ; m.β]
deviance(m::PolrModel) = -2polrfun!(m, false, false)
dof_residual(m::PolrModel) = m.n - m.npar
fitted(m::PolrModel) = nothing # TODO
loglikelihood(m::PolrModel) = polrfun!(m, false, false)
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
              [["θ$i" for i = 1:m.J-1]; ["β$i" for i = 1:m.p]], 4)
end

confint(m::PolrModel, level::Real) = hcat(coef(m), coef(m)) +
    stderr(m) * quantile(Normal(), (1. - level) / 2.) * [1. -1.]
confint(m::PolrModel) = confint(m, 0.95)

function cor(m::PolrModel)
    Σ = vcov(m)
    invstd = similar(Σ, size(Σ, 1))
    for i in eachindex(invstd)
        invstd[i] = 1 / sqrt(Σ[i, i])
    end
    scale!(invstd, scale!(Σ, invstd))
end

include("polrrand.jl")
include("polrfit.jl")
include("polrtest.jl")

end # module
