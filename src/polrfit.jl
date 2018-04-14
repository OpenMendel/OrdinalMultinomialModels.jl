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
    A_mul_B!(m.η, m.X, m.β)
    # update θ according to α
    # α[2] to α[J-2] overwritten by their exponential
    fill!(m.dθdα, 0)
    m.θ[1] = m.α[1]
    m.dθdα[1, :] = 1
    for j in 2:m.J-1
        m.θ[j] = m.θ[j-1] + exp(m.α[j])
        m.dθdα[j, j:end] = exp(m.α[j])
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
        ### gradent
        needgrad || continue
        # derivative of θ
        pyinv  = inv(py)
        deriv1 = pyinv * mueta(m.link, ζ1)
        deriv2 = pyinv * mueta(m.link, ζ2)
        derivΔ = deriv1 - deriv2
        yobs < m.J && (m.∇[yobs  ] += wtobs * deriv1)
        yobs > 1   && (m.∇[yobs-1] -= wtobs * deriv2)
        # derivative of β = wt[obs] * X[:. obs]
        m.wt∂β[obs] = - derivΔ
        ### hessian
        needhess || continue
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
    if needgrad
        if isempty(m.wts)
            @views At_mul_B!(m.∇[m.J:end], m.X, m.wt∂β)

        else
            m.wtwk .= m.wts .* m.wt∂β
            @views At_mul_B!(m.∇[m.J:end], m.X, m.wtwk)
        end
    end
    if needhess
        if isempty(m.wts)
            @views At_mul_B!(m.H[m.J:end, 1:m.J-1], m.X, m.wt∂θ∂β)
            copy!(m.scratchm1, m.X)
            scale!(m.wt∂β∂β, m.scratchm1)
            @views At_mul_B!(m.H[m.J:end, m.J:end], m.X, m.scratchm1)
        else
            copy!(m.scratchm1, m.X)
            scale!(m.wts, m.scratchm1)
            @views At_mul_B!(m.H[m.J:end, 1:m.J-1], m.scratchm1, m.wt∂θ∂β)
            m.wtwk .= m.wts .* m.wt∂β∂β
            copy!(m.scratchm1, m.X)
            scale!(m.wtwk, m.scratchm1)
            @views At_mul_B!(m.H[m.J:end, m.J:end], m.X, m.scratchm1)
        end
    LinAlg.copytri!(m.H, 'L')
    end
    logl
end

polr(X, y, args...; kwargs...) = fit(AbstractPolrModel, X, y, args...; kwargs...)

"""
    fit(AbstractPolrModel, X, y, link, solver)

Fit ordered multinomial model by maximum likelihood estimation.

# Input

* `y::Vector`: integer vector taking values in `1,...,J`.
* `X::Matrix`: `p x n` covariate matrix excluding intercept.
* `link::GLM.Link`: `LogitLink()`, `ProbitLink()`, `CauchitLink()`, or `CloglogLink()`.
* `solver`: `IpoptSolver()` or `NLoptSolver()`

# Output

* `dd:PolrModel`: a `PolrModel` type.
"""
function fit(
    ::Type{M},
    X::AbstractMatrix,
    y::AbstractVector,
    link::GLM.Link = LogitLink(),
    solver = NLoptSolver(algorithm=:LD_SLSQP);
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
    stat == :Optimal || warn("Optimization unsuccesful; got $stat")
    xsol = MathProgBase.getsolution(m)
    copy!(dd.α, 1, xsol, 1, dd.J - 1)
    copy!(dd.β, 1, xsol, dd.J, dd.p)
    polrfun!(dd, true, true)
    dd.vcov[:] = inv(-dd.H)
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
    copy!(m.α, 1, par, 1, m.J - 1)
    copy!(m.β, 1, par, m.J, m.p)
    polrfun!(m, false, false)
end

function MathProgBase.eval_grad_f(m::PolrModel, grad::Vector, par::Vector)
    copy!(m.α, 1, par, 1, m.J - 1)
    copy!(m.β, 1, par, m.J, m.p)
    polrfun!(m, true, false)
    @views A_mul_B!(grad[1:m.J-1], m.dθdα, m.∇[1:m.J-1])
    copy!(grad, m.J, m.∇, m.J, m.p)
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
