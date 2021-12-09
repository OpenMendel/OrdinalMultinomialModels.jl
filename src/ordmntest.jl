function polrtest(nm::OrdinalMultinomialModel, Z::AbstractVecOrMat; test=:score)
    if test == :score
        polrtest(OrdinalMultinomialScoreTest(nm, reshape(Z, size(Z, 1), size(Z, 2))))
    elseif test == :LRT
        polrtest(OrdinalMultinomialLrtTest(nm, reshape(Z, size(Z, 1), size(Z, 2))))
    else
        throw(ArgumentError("unrecognized test $test"))
    end
end

polrtest(nm::StatsModels.TableRegressionModel{<:OrdinalMultinomialModel}, Z::AbstractVecOrMat; kwargs...) = 
polrtest(nm.model, Z; kwargs...)

###########################
# Score test
###########################

struct OrdinalMultinomialScoreTest{TY<:Integer, T<:BlasReal, TL<:GLM.Link}
    "`q`: number of covariates to test significance"
    q::Int
    "`Z`: n-by-q covariate matrix to test significance"
    Z::Matrix{T}
    "`nm`: fitted null model"
    nm::OrdinalMultinomialModel{TY,T,TL}
    "`scoreγ`: score vector of γ"
    scoreγ::Vector{T}
    "`∂γ∂θβ`: information of γ vs (θ,β)"
    ∂γ∂θβ::Matrix{T}
    "`∂γ∂γ`: information of γ vs γ"
    ∂γ∂γ::Matrix{T}
    "`scratchv1`: working vector, same size as `scoreγ`"
    scratchv1::Vector{T}
    "`scratchm1`: working matrix, same size as `Z`"
    scratchm1::Matrix{T}
    "`scratchm2`: working matrix, same size as `∂γ∂θβ`"
    scratchm2::Matrix{T}
    "`scratchm3`: working matrix, same size as `∂γ∂γ`"
    scratchm3::Matrix{T}
end

function OrdinalMultinomialScoreTest(nm::OrdinalMultinomialModel, Z::Matrix)
    TY, T, TL = eltype(nm.Y), eltype(nm.X), typeof(nm.link)
    q = size(Z, 2)
    scoreγ = zeros(T, q)
    ∂γ∂θβ  = zeros(T, q, nm.npar)
    ∂γ∂γ   = zeros(T, q, q)
    scratchv1 = similar(scoreγ)
    scratchm1 = similar(Z)
    scratchm2 = similar(∂γ∂θβ)
    scratchm3 = similar(∂γ∂γ)
    OrdinalMultinomialScoreTest{TY, T, TL}(q, Z, nm, scoreγ, ∂γ∂θβ, ∂γ∂γ, scratchv1, scratchm1, scratchm2, scratchm3)
end


function polrtest(d::OrdinalMultinomialScoreTest)
    #
    # partition the informatrix matrix as
    # [P  W]
    # [W' Q],
    # test statistics is R^T (Q - W' P^{-1} W)^{-1} * R, where R is the score
    #
    # compute R = d.scoreγ, score of γ
    if isempty(d.nm.wts)
        mul!(d.scoreγ, transpose(d.Z), d.nm.wt∂β)
    else
        d.nm.wtwk .= d.nm.wts .* d.nm.wt∂β
        mul!(d.scoreγ, transpose(d.Z), d.nm.wtwk)
    end
    # compute W' = d.∂γ∂θβ
    if isempty(d.nm.wts)
        @views mul!(d.∂γ∂θβ[:, 1:d.nm.J-1], transpose(d.Z), d.nm.wt∂θ∂β) # wt∂θ∂β is n-by-(J-1)
    else
        copyto!(d.scratchm1, d.Z)
        lmul!(Diagonal(d.nm.wts), d.scratchm1)
        @views mul!(d.∂γ∂θβ[:, 1:d.nm.J-1], transpose(d.scratchm1), d.nm.wt∂θ∂β)
    end
    if isempty(d.nm.wts)
        copyto!(d.nm.wtwk, d.nm.wt∂β∂β)
    else
        d.nm.wtwk .= d.nm.wts .* d.nm.wt∂β∂β
    end
    copyto!(d.nm.scratchm1, d.nm.X)
    lmul!(Diagonal(d.nm.wtwk), d.nm.scratchm1)
    @views mul!(d.∂γ∂θβ[:, d.nm.J:end], transpose(d.Z), d.nm.scratchm1)
    # compute Q = ∂γ∂γ
    copyto!(d.scratchm1, d.Z)
    lmul!(Diagonal(d.nm.wtwk), d.scratchm1)
    mul!(d.∂γ∂γ, transpose(d.Z), d.scratchm1)
    # compute scratchm3 = Q - W' P^{-1} W
    mul!(d.scratchm2, d.∂γ∂θβ, d.nm.vcov)
    mul!(d.scratchm3, d.scratchm2, transpose(d.∂γ∂θβ))
    d.scratchm3 .= d.∂γ∂γ .- d.scratchm3
    # compute ts = R^T (Q - W' P^{-1} W)^{-1} R
    eigfact = eigen!(Symmetric(d.scratchm3))
    mul!(d.scratchv1, transpose(eigfact.vectors), d.scoreγ)
    T = eltype(d.scoreγ)
    atol = 1e-8 # tolerance for determining rank
    rk = 0 # rank
    ts = zero(T)
    for j in 1:d.q
        if eigfact.values[j] > atol
            ts += abs2(d.scratchv1[j]) / eigfact.values[j]
            rk += 1
        end
    end
    pval = ts ≤ 0 ? 1.0 : ccdf(Chisq(rk), ts)
end


###########################
# LRT test
###########################

struct OrdinalMultinomialLrtTest{TY<:Integer, T<:BlasReal, TL<:GLM.Link}
    "`q`: number of covariates to test significance"
    q::Int
    "`Z`: n-by-q covariate matrix to test significance"
    Z::Matrix{T}
    "`nm`: fitted null model"
    nm::OrdinalMultinomialModel{TY,T,TL}
    "`Xaug`: augmented X matrix for full model"
    Xaug::Matrix{T}
end

function OrdinalMultinomialLrtTest(nm::OrdinalMultinomialModel, Z::AbstractVecOrMat)
    TY, T, TL = eltype(nm.Y), eltype(nm.X), typeof(nm.link)
    q = size(Z, 2)
    Xaug = Matrix{T}(undef, nm.n, nm.p + q)
    @views copyto!(Xaug[:, 1:nm.p], nm.X)
    @views copyto!(Xaug[:, nm.p+1:end], Z)
    OrdinalMultinomialLrtTest{TY, T, TL}(q, Z, nm, Xaug)
end

function polrtest(
    d::OrdinalMultinomialLrtTest; 
    solver::TSOLVER = NLopt.Optimizer()
    ) where {TSOLVER <: MOI.AbstractOptimizer}
    if typeof(solver) == NLopt.Optimizer
        algo = MOI.get(solver, MOI.RawOptimizerAttribute("algorithm"))
        if algo == :none
            MOI.set(solver, MOI.RawOptimizerAttribute("algorithm"), :LD_SLSQP)
            MOI.set(solver, MOI.RawOptimizerAttribute("max_iter"), 4000)
        end
    end
    am = polr(d.Xaug, d.nm.Y, d.nm.link, solver; wts = d.nm.wts)
    ccdf(Chisq(d.q), deviance(d.nm) - deviance(am))
end
