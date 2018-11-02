function polrtest(nm::PolrModel, Z::AbstractVecOrMat; test=:score)
    if test == :score
        polrtest(PolrScoreTest(nm, reshape(Z, size(Z, 1), size(Z, 2))))
    elseif test == :LRT
        polrtest(PolrLrtTest(nm, reshape(Z, size(Z, 1), size(Z, 2))))
    else
        throw(ArgumentError("unrecognized test $test"))
    end
end

polrtest(nm::StatsModels.DataFrameRegressionModel{<:PolrModel}, Z::AbstractVecOrMat; kwargs...) = 
polrtest(nm.model, Z; kwargs...)

###########################
# Score test
###########################

struct PolrScoreTest{TY<:Integer, T<:BlasReal, TL<:GLM.Link}
    "`q`: number of covariates to test significance"
    q::Int
    "`Z`: n-by-q covariate matrix to test significance"
    Z::Matrix{T}
    "`nm`: fitted null model"
    nm::PolrModel{TY,T,TL}
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

function PolrScoreTest(nm::PolrModel, Z::Matrix)
    TY, T, TL = eltype(nm.Y), eltype(nm.X), typeof(nm.link)
    q = size(Z, 2)
    scoreγ = zeros(T, q)
    ∂γ∂θβ  = zeros(T, q, nm.npar)
    ∂γ∂γ   = zeros(T, q, q)
    scratchv1 = similar(scoreγ)
    scratchm1 = similar(Z)
    scratchm2 = similar(∂γ∂θβ)
    scratchm3 = similar(∂γ∂γ)
    PolrScoreTest{TY, T, TL}(q, Z, nm, scoreγ, ∂γ∂θβ, ∂γ∂γ, scratchv1, scratchm1, scratchm2, scratchm3)
end


function polrtest(d::PolrScoreTest)
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
    d.∂γ∂θβ .*= -1
    # compute Q = ∂γ∂γ
    copyto!(d.scratchm1, d.Z)
    lmul!(Diagonal(d.nm.wtwk), d.scratchm1)
    mul!(d.∂γ∂γ, transpose(d.Z), d.scratchm1)
    d.∂γ∂γ .*= -1
    # compute scratchm3 = Q - W' P^{-1} W
    mul!(d.scratchm2, d.∂γ∂θβ, d.nm.vcov)
    mul!(d.scratchm3, d.scratchm2, transpose(d.∂γ∂θβ))
    d.scratchm3 .= d.∂γ∂γ .- d.scratchm3
    # compute ts = R^T (Q - W' P^{-1} W)^{-1} R
    ldiv!(d.scratchv1, bunchkaufman!(Symmetric(d.scratchm3)), d.scoreγ)
    ts = dot(d.scoreγ, d.scratchv1)
    # p-value
    ts ≤ 0 ? 1.0 : ccdf(Chisq(d.q), ts)
end


###########################
# LRT test
###########################

struct PolrLrtTest{TY<:Integer, T<:BlasReal, TL<:GLM.Link}
    "`q`: number of covariates to test significance"
    q::Int
    "`Z`: n-by-q covariate matrix to test significance"
    Z::Matrix{T}
    "`nm`: fitted null model"
    nm::PolrModel{TY,T,TL}
    "`Xaug`: augmented X matrix for full model"
    Xaug::Matrix{T}
end

function PolrLrtTest(nm::PolrModel, Z::AbstractVecOrMat)
    TY, T, TL = eltype(nm.Y), eltype(nm.X), typeof(nm.link)
    q = size(Z, 2)
    Xaug = Matrix{T}(undef, nm.n, nm.p + q)
    @views copyto!(Xaug[:, 1:nm.p], nm.X)
    @views copyto!(Xaug[:, nm.p+1:end], Z)
    PolrLrtTest{TY, T, TL}(q, Z, nm, Xaug)
end

function polrtest(
    d::PolrLrtTest; 
    solver = NLoptSolver(algorithm=:LD_SLSQP, maxeval=4000))
    am = polr(d.Xaug, d.nm.Y, d.nm.link, solver; wts = d.nm.wts)
    ccdf(Chisq(d.q), deviance(d.nm) - deviance(am))
end
