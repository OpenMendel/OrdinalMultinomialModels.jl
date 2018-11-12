"""
polrfun!(m, [needgrad=false], [needhess=false])

Evaluate the log-likelihood and optionally gradient and Hessian of
ordered multinomial model.

# Input

* `m::PolrModel`: a `PolrModel` type. Log-likelihood is evaluated based on
fields `Y`, `X`, `α`, `β`, and `link` of `m`.

# Output

* `logl`: log-likelihood. If `needgrad=true`, field `m.∇` is overwritten by the
gradient (score) with respect to `(θ, β)`. If `needhess=true`, field `m.H` is
overwritten by the Hessian with respect to `(θ, β)`.
"""
function polrfun!(
    m::PolrModel,
    needgrad = false,
    needhess = false)
    
    T = eltype(m.X)
    needgrad && fill!(m.∇, 0)
    needhess && fill!(m.H, 0)
    # update η according to X and β
    mul!(m.η, m.X, m.β)
    # update θ according to α
    # α[2] to α[J-2] overwritten by their exponential
    fill!(m.dθdα, 0)
    m.θ[1] = m.α[1]
    m.dθdα[1, :] .= 1
    for j in 2:m.J-1
        m.θ[j] = m.θ[j-1] + exp(m.α[j])
        m.dθdα[j, j:end] .= exp(m.α[j])
    end
    minζ, maxζ = convert(T, -100), convert(T, 100)
    logl = zero(T)
    for obs in 1:m.n
        wtobs = isempty(m.wts) ? one(T) : m.wts[obs]
        wtobs == 0 && continue
        ### log-likelihood
        yobs = m.Y[obs]
        ζ1 = yobs == m.J ? maxζ : m.θ[yobs]   - m.η[obs]
        ζ2 = yobs == 1   ? minζ : m.θ[yobs-1] - m.η[obs]
        ζ1 = clamp(ζ1, minζ, maxζ) # for numerical stability
        ζ2 = clamp(ζ2, minζ, maxζ)
        cumprob1 = GLM.linkinv(m.link, ζ1)
        cumprob2 = GLM.linkinv(m.link, ζ2)
        py = cumprob1 - cumprob2
        logl += py ≤ 0 ? typemin(T) : wtobs * log(py)
        ### gradient
        needgrad || continue
        if py ≤ 0
            # derivative of θ
            yobs < m.J && (m.∇[yobs  ] += typemax(T))
            yobs > 1   && (m.∇[yobs-1] += typemin(T))
            # derivative of β = wt[obs] * X[:. obs]
            m.wt∂β[obs] = zero(T)
        else
            # derivative of θ
            pyinv  = inv(py)
            deriv1 = pyinv * GLM.mueta(m.link, ζ1)
            deriv2 = pyinv * GLM.mueta(m.link, ζ2)
            derivΔ = deriv1 - deriv2
            yobs < m.J && (m.∇[yobs  ] += wtobs * deriv1)
            yobs > 1   && (m.∇[yobs-1] -= wtobs * deriv2)
            # derivative of β = wt[obs] * X[:. obs]
            m.wt∂β[obs] = - derivΔ
        end
        ### hessian
        needhess || continue
        if py > 0
            h1 = pyinv * muetad2(m.link, ζ1)
            h2 = pyinv * muetad2(m.link, ζ2)
            # hessian for ∂θ^2
            if     yobs > 1  ; m.H[yobs-1, yobs-1] -= wtobs * (h2 + deriv2 * deriv2); end
            if 1 < yobs < m.J; m.H[yobs  , yobs-1] += wtobs * (deriv1 * deriv2); end
            if     yobs < m.J; m.H[yobs  , yobs  ] += wtobs * (h1 - deriv1 * deriv1); end
            # hessian for ∂θ∂β
            yobs > 1   && (m.wt∂θ∂β[obs, yobs - 1] += h2 - deriv2 * derivΔ)
            yobs < m.J && (m.wt∂θ∂β[obs, yobs    ] -= h1 - deriv1 * derivΔ)
            # hessian for ∂β∂β
            m.wt∂β∂β[obs] += h1 - h2 - derivΔ * derivΔ
        end
    end
    if needgrad
        if isempty(m.wts)
            @views mul!(m.∇[m.J:end], transpose(m.X), m.wt∂β)
        else
            m.wtwk .= m.wts .* m.wt∂β
            @views mul!(m.∇[m.J:end], transpose(m.X), m.wtwk)
        end
    end
    if needhess
        if isempty(m.wts)
            @views mul!(m.H[m.J:end, 1:m.J-1], transpose(m.X), m.wt∂θ∂β)
            copyto!(m.scratchm1, m.X)
            lmul!(Diagonal(m.wt∂β∂β), m.scratchm1)
            @views mul!(m.H[m.J:end, m.J:end], transpose(m.X), m.scratchm1)
        else
            copyto!(m.scratchm1, m.X)
            lmul!(Diagonal(m.wts), m.scratchm1)
            @views mul!(m.H[m.J:end, 1:m.J-1], transpose(m.scratchm1), m.wt∂θ∂β)
            m.wtwk .= m.wts .* m.wt∂β∂β
            copyto!(m.scratchm1, m.X)
            lmul!(Diagonal(m.wtwk), m.scratchm1)
            @views mul!(m.H[m.J:end, m.J:end], transpose(m.X), m.scratchm1)
        end
        LinearAlgebra.copytri!(m.H, 'L')
    end
    logl
end

polr(X, y, args...; kwargs...) = fit(AbstractPolrModel, X, y, args...; kwargs...)

"""
fit(AbstractPolrModel, X, y, link, solver)

Fit ordered multinomial model by maximum likelihood estimation.

# Input
- `y::Vector`: integer vector taking values in `1,...,J`.
- `X::Matrix`: `n x p` covariate matrix excluding intercept.
- `link::GLM.Link`: `LogitLink()`, `ProbitLink()`, `CauchitLink()`, or `CloglogLink()`.
- `solver`: `IpoptSolver()` or `NLoptSolver()`.

# Output
- `dd:PolrModel`: a `PolrModel` type.
"""
function fit(
    ::Type{M},
    X::AbstractMatrix,
    y::AbstractVector,
    link::GLM.Link = LogitLink(),
    solver = NLoptSolver(algorithm=:LD_SLSQP, maxeval=4000);
    wts::AbstractVector = similar(X, 0)
    ) where M <: AbstractPolrModel
    
    ydata = denserank(y)
    dd = PolrModel(X, ydata, convert(Vector{eltype(X)}, wts), link)
    m = MathProgBase.NonlinearModel(solver)
    lb = fill(-Inf, dd.npar)
    ub = fill( Inf, dd.npar)
    MathProgBase.loadproblem!(m, dd.npar, 0, lb, ub, Float64[], Float64[], :Max, dd)
    # initialize from LS solution
    β0 = [ones(length(y)) X] \ ydata
    par0 = [β0[1] - dd.J / 2 + 1; zeros(dd.J - 2); β0[2:end]]
    MathProgBase.setwarmstart!(m, par0)
    MathProgBase.optimize!(m)
    
    # ouput
    stat = MathProgBase.status(m)
    stat == :Optimal || @warn("Optimization unsuccesful; got $stat")
    xsol = MathProgBase.getsolution(m)
    copyto!(dd.α, 1, xsol, 1, dd.J - 1)
    copyto!(dd.β, 1, xsol, dd.J, dd.p)
    polrfun!(dd, true, true)
    Hbk = bunchkaufman!(Symmetric(dd.H), check=false)
    if issuccess(Hbk)
        dd.vcov[:] = -inv(Hbk)
    else # Hessian is singular
        dd.vcov .= Inf
    end
    return dd
end

function MathProgBase.initialize(m::PolrModel,
    requested_features::Vector{Symbol})
    for feat in requested_features
        if !(feat in [:Grad])
            error("Unsupported feature $feat")
        end
    end
end

MathProgBase.features_available(m::PolrModel) = [:Grad]

function MathProgBase.eval_f(m::PolrModel, par::Vector)
    copyto!(m.α, 1, par, 1, m.J - 1)
    copyto!(m.β, 1, par, m.J, m.p)
    polrfun!(m, false, false)
end

function MathProgBase.eval_grad_f(m::PolrModel, grad::Vector, par::Vector)
    copyto!(m.α, 1, par, 1, m.J - 1)
    copyto!(m.β, 1, par, m.J, m.p)
    polrfun!(m, true, false)
    @views mul!(grad[1:m.J-1], m.dθdα, m.∇[1:m.J-1])
    copyto!(grad, m.J, m.∇, m.J, m.p)
end

"""
muetad2(::CuachitLink, η)

`d^2μ/dη^2` for CauchitLink.
"""
muetad2(::CauchitLink, η) = - 2η / (pi * (one(η) + abs2(η)))^2

"""
muetad2(::CloglogLink, η)

`d^2μ/dη^2` for CloglogLink.
"""
function muetad2(::CloglogLink, η::T) where T<:Real
    expη = exp(η)
    expη * exp(-expη) * (one(η) - expη)
end

"""
muetad2(::LogitLink, η)

`d^2μ/dη^2` for LogitLink.
"""
function muetad2(::LogitLink, η)
    expabs = exp(-abs(η))
    denominv = inv(1 + expabs)
    term1 = expabs * denominv * denominv
    η ≥ 0 ? term1 * (1 - 2denominv) : term1 * (1 - 2expabs * denominv)
end

"""
muetad2(::ProbitLink, η)

`d^2μ/dη^2` for ProbitLink.
"""
muetad2(::ProbitLink, η) = - (exp(-abs2(η) / 2) / GLM.sqrt2π) * η
