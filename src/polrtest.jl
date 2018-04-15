struct PolrScoreTest{TY<:Integer, T<:BlasReal, TL<:GLM.Link}
    "`q`: number of covariates to test significance"
    q::Int
    "`Z`: covariates to test significance"
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

polrtest(nm::PolrModel, Z::AbstractMatrix) = polrtest(PolrScoreTest(nm, Z))

polrtest(nm::PolrModel, Z::AbstractVector) = polrtest(nm, reshape(Z, length(Z), 1))

function polrtest(d::PolrScoreTest)
    At_mul_B!(d.scoreγ, d.Z, d.nm.wt∂β)
    @views At_mul_B!(d.∂γ∂θβ[:, 1:d.nm.J-1], d.Z, d.nm.wt∂θ∂β)
    @views At_mul_B!(d.∂γ∂θβ[:, d.nm.J:end], d.Z, d.nm.scratchm1)
    d.∂γ∂θβ .*= -1
    copy!(d.scratchm1, d.Z)
    scale!(d.nm.wt∂β∂β, d.scratchm1)
    At_mul_B!(d.∂γ∂γ, d.Z, d.scratchm1)
    d.∂γ∂γ .*= -1
    A_mul_B!(d.scratchm2, d.∂γ∂θβ, d.nm.vcov)
    A_mul_Bt!(d.scratchm3, d.scratchm2, d.∂γ∂θβ)
    d.scratchm3 .= d.∂γ∂γ .- d.scratchm3
    cf = cholfact!(Symmetric(d.scratchm3), Val{true})
    if rank(cf) < d.q; return 1.0; end
    A_ldiv_B!(d.scratchv1, cf, d.scoreγ)
    ts = dot(d.scoreγ, d.scratchv1)
    ccdf(Chisq(d.q), ts)
end
