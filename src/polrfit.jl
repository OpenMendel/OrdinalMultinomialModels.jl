"""
    polrfun!(m, [needgrad=false], [needhess=false])

Evaluate the log-likelihood and optionally gradient and Hessian of
ordered multinomial model.

# Positional arguments

- `m::PolrModel`: a `PolrModel` type. Log-likelihood is evaluated based on
fields `Y`, `X`, `α`, `β`, and `link` of `m`.  
- `needgrad::Bool=false`: evaluate gradient or not.  
- `needhess::Bool=false`: evaluate Hessian or not.

# Output

- `logl`: log-likelihood. If `needgrad=true`, field `m.∇` is overwritten by the
gradient (score) with respect to `(θ, β)`. If `needhess=true`, field `m.H` is
overwritten by the Hessian with respect to `(θ, β)`.
"""
function polrfun!(
    m::PolrModel,
    needgrad = false,
    needhess = false)
    
    T = eltype(m.X)
    needgrad && fill!(m.∇, 0)
    needhess && fill!(m.FIM, 0)
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
    minη, maxη = convert(T, -100), convert(T, 100)
    logl = zero(T)
    for obs in 1:m.n
        wtobs = isempty(m.wts) ? one(T) : m.wts[obs]
        wtobs == 0 && continue
        ### log-likelihood
        yobs = m.Y[obs]
        ηy   = yobs == m.J ? maxη : m.θ[yobs]   - m.η[obs]
        ηym1 = yobs == 1   ? minη : m.θ[yobs-1] - m.η[obs]
        ηy   = clamp(ηy  , minη, maxη) # for numerical stability
        ηym1 = clamp(ηym1, minη, maxη) # ym1 reads y-minus-one
        cumproby   = GLM.linkinv(m.link, ηy)
        cumprobym1 = GLM.linkinv(m.link, ηym1)
        py = cumproby - cumprobym1
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
            pyinv    = inv(py)
            derivy   = GLM.mueta(m.link, ηy)
            derivym1 = GLM.mueta(m.link, ηym1)
            derivΔ = derivy - derivym1
            yobs < m.J && (m.∇[yobs  ] += wtobs * pyinv * derivy)
            yobs > 1   && (m.∇[yobs-1] -= wtobs * pyinv * derivym1)
            # derivative of β = wt[obs] * X[:. obs]
            m.wt∂β[obs] = - pyinv * derivΔ
        end
        ### hessian
        needhess || continue
        ηjm1 = ηj = ηjp1 = zero(T)
        cumprobjm1 = cumprobj = cumprobjp1 = zero(T)
        derivjm1 = derivj = derivjp1 = zero(T)
        for j in 1:m.J-1
            # systematic components at j-1, j, j+1
            ηjm1       = j == 1 ? minη : ηj
            ηj         = j == 1 ? m.θ[1] - m.η[obs] : ηjp1
            ηjp1       = j == m.J-1 ? maxη : m.θ[j+1] - m.η[obs]
            # cumulative probabilities at j-1, j, j+1
            cumprobjm1 = j == 1 ? zero(T) : cumprobj
            cumprobj   = j == 1 ? GLM.linkinv(m.link, ηj) : cumprobjp1
            cumprobjp1 = j == m.J-1 ? one(T) : GLM.linkinv(m.link, ηjp1)
            # point mass at j-1, j, j+1
            pj         = cumprobj - cumprobjm1
            pjp1       = cumprobjp1 - cumprobj
            pjinv      =   pj ≤ 0 ? zero(T) : inv(pj)
            pjp1inv    = pjp1 ≤ 0 ? zero(T) : inv(pjp1)
            # mueta at j-1, j, j+1
            derivjm1   = j == 1 ? zero(T) : derivj
            derivj     = j == 1 ? GLM.mueta(m.link, ηj) : derivjp1
            derivjp1   = j == m.J-1 ? zero(T) : GLM.mueta(m.link, ηjp1)
            # FIM for ∂θ^2
            m.FIM[j, j] += wtobs * (pjinv + pjp1inv) * derivj * derivj
            if j > 1
                m.FIM[j, j-1] -= wtobs * pjinv * derivj * derivjm1
            end
            # FIM for ∂θ∂β, weights only
            m.wt∂θ∂β[obs, j] += (pjinv * derivjm1 - (pjinv + pjp1inv) * derivj +
                pjp1inv * derivjp1) * derivj
            # FIM for ∂β∂β, weights only
            m.wt∂β∂β[obs] += pjinv * (derivj - derivjm1)^2
            j == m.J - 1 &&  (m.wt∂β∂β[obs] += pjp1inv * derivj * derivj) # p_iJ
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
            # ∂θ∂β
            @views mul!(m.FIM[m.J:end, 1:m.J-1], transpose(m.X), m.wt∂θ∂β)
            # ∂β∂β
            copyto!(m.scratchm1, m.X)
            lmul!(Diagonal(m.wt∂β∂β), m.scratchm1)
            @views mul!(m.FIM[m.J:end, m.J:end], transpose(m.X), m.scratchm1)
        else
            # ∂θ∂β
            copyto!(m.scratchm1, m.X)
            lmul!(Diagonal(m.wts), m.scratchm1)
            @views mul!(m.FIM[m.J:end, 1:m.J-1], transpose(m.scratchm1), m.wt∂θ∂β)
            # ∂β∂β
            m.wtwk .= m.wts .* m.wt∂β∂β
            copyto!(m.scratchm1, m.X)
            lmul!(Diagonal(m.wtwk), m.scratchm1)
            @views mul!(m.FIM[m.J:end, m.J:end], transpose(m.X), m.scratchm1)
        end
        LinearAlgebra.copytri!(m.FIM, 'L') # copy lower triangular to upper triangular
    end
    logl
end

"""
    polr(formula, df, link, solver=NLoptSolver(algorithm=:LD_SLSQP, maxeval=4000))
    polr(X, y, link, solver=NLoptSolver(algorithm=:LD_SLSQP, maxeval=4000))

Fit ordered multinomial model by maximum likelihood estimation.

# Positional arguments

- `formula::Formula`: a model formula specifying responses and regressors.
- `df::DataFrame`: a dataframe. Response variable should take integer values 
    starting from 1.
- `y::Vector`: integer vector taking values in `1,...,J`.
- `X::Matrix`: `n x p` covariate matrix excluding intercept.
- `link::GLM.Link`: `LogitLink()` (default), `ProbitLink()`, `CauchitLink()`, or `CloglogLink()`.
- `solver`: `NLoptSolver()` (default) or `IpoptSolver()`.

# Keyword arguments

- `wts::AbstractVector=similar(X, 0)`: observation weights.

# Output
- `dd:PolrModel`: a `PolrModel` type.
"""
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
    FIMbk = bunchkaufman(Symmetric(dd.FIM), check=false)
    if issuccess(FIMbk)
        dd.vcov[:] = inv(FIMbk)
    else # FIM is singular
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
