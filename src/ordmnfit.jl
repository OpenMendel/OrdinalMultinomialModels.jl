"""
    loglikelihood!(m, [needgrad=false], [needhess=false])

Evaluate the log-likelihood and optionally gradient and Hessian of
ordered multinomial model.

# Positional arguments
- `m::PolrModel`: a `PolrModel` type. Log-likelihood is evaluated based on
fields `Y`, `X`, `α`, `β`, and `link` of `m`.  
- `needgrad::Bool=false`: evaluate gradient or not.  
- `needhess::Bool=false`: evaluate Fisher information matrix (FIM) (negative Hessian) or not.

# Output
- `logl`: log-likelihood. If `needgrad=true`, field `m.∇` is overwritten by the
gradient (score) with respect to `(θ, β)`. If `needhess=true`, field `m.FIM` is
overwritten by the Fisher information matrix (FIM), equivalent to negative Hessian, 
with respect to `(θ, β)`.
"""
function loglikelihood!(
    m::OrdinalMultinomialModel,
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
            j > 1 && (m.FIM[j, j-1] -= wtobs * pjinv * derivj * derivjm1)
            # FIM for ∂θ∂β, weights only
            m.wt∂θ∂β[obs, j] += (pjinv * derivjm1 - (pjinv + pjp1inv) * derivj +
                pjp1inv * derivjp1) * derivj
            # FIM for ∂β∂β, weights only
            m.wt∂β∂β[obs] += pjinv * (derivj - derivjm1)^2
            j == m.J - 1 &&  (m.wt∂β∂β[obs] += pjp1inv * derivj * derivj) # contribution from p_iJ
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
- `solver`: An instance of `Ipopt.Optimizer()` or `NLopt.Optimizer()` (default).

# Keyword arguments
- `wts::AbstractVector=similar(X, 0)`: observation weights.

# Output
- `dd:PolrModel`: a `PolrModel` type.
"""
polr(X::AbstractMatrix, y::AbstractVector, args...; kwargs...) = 
    fit(AbstractOrdinalMultinomialModel, X, y, args...; kwargs...)
polr(f::FormulaTerm, df, args...; kwargs...) = 
    fit(AbstractOrdinalMultinomialModel, f, df, args...; kwargs...)

"""
    fit(AbstractOrdinalMultinomialModel, X, y, link, solver)

Fit ordered multinomial model by maximum likelihood estimation.

# Input
- `y::Vector`: integer vector taking values in `1,...,J`.
- `X::Matrix`: `n x p` covariate matrix excluding intercept.
- `link::GLM.Link`: `LogitLink()`, `ProbitLink()`, `CauchitLink()`, or `CloglogLink()`.
- `solver`: An instance of `Ipopt.Optimizer()` or `NLopt.Optimizer()` (default).

# Output
- `dd:OrdinalMultinomialModel`: an `OrdinalMultinomialModel` type.
"""
function fit(
    ::Type{M},
    X::AbstractMatrix,
    y::AbstractVecOrMat,
    link::GLM.Link = LogitLink(),
    solver::TSOLVER = NLopt.Optimizer();
    wts::AbstractVector = similar(X, 0)
    ) where {M<:AbstractOrdinalMultinomialModel, TSOLVER<:MOI.AbstractOptimizer}
    T = eltype(X)
    ydata = Vector{Int}(undef, size(y, 1))
    # set up optimization
    # This is a hack since we can't initialize a solver with parameters already set
    if typeof(solver) == NLopt.Optimizer
        algo = MOI.get(solver, MOI.RawOptimizerAttribute("algorithm"))
        if algo == :none
            MOI.set(solver, MOI.RawOptimizerAttribute("algorithm"), :LD_SLSQP)
            MOI.set(solver, MOI.RawOptimizerAttribute("max_iter"), 4000)
        end
    end
    
    if size(y, 2) == 1
        ydata = denserank(y)
    else #y is encoded via dummy-encoding 
        for i in 1:size(y, 1)
            idx = findfirst(view(y, i, :) .== 1)
            ydata[i] = idx == nothing ? 1 : idx + 1
        end
    end
    #ydata = denserank(y) # dense ranking of y, http://juliastats.github.io/StatsBase.jl/stable/ranking.html#StatsBase.denserank
    dd = OrdinalMultinomialModel(X, ydata, convert(Vector{eltype(X)}, wts), link)
    lb = T[]
    ub = T[]
    NLPBlock = MOI.NLPBlockData(
        MOI.NLPBoundsPair.(lb, ub), dd, true
    )
    # initialize from LS solution
    β0 = [ones(length(ydata)) X] \ ydata
    par0 = [β0[1] - dd.J / 2 + 1; zeros(dd.J - 2); β0[2:end]]
    
    params = MOI.add_variables(solver, dd.npar)
    for i in 1:dd.npar
        MOI.set(solver, MOI.VariablePrimalStart(), params[i], par0[i])
    end

    MOI.set(solver, MOI.NLPBlock(), NLPBlock)
    MOI.set(solver, MOI.ObjectiveSense(), MOI.MAX_SENSE)

    MOI.optimize!(solver)
    # output
    stat = MOI.get(solver, MOI.TerminationStatus())
    stat == MOI.LOCALLY_SOLVED || @warn("Optimization unsuccessful; got $stat")
    xsol = similar(par0)
    for i in eachindex(xsol)
        xsol[i] = MOI.get(solver, MOI.VariablePrimal(), MOI.VariableIndex(i))
    end
    copyto!(dd.α, 1, xsol, 1, dd.J - 1)
    copyto!(dd.β, 1, xsol, dd.J, dd.p)
    loglikelihood!(dd, true, true)
    FIMbk = bunchkaufman(Symmetric(dd.FIM), check=false)
    if issuccess(FIMbk)
        dd.vcov[:] = inv(FIMbk)
    else # FIM is singular
        dd.vcov .= Inf
    end
    return dd

end

function MathOptInterface.initialize(
    m::OrdinalMultinomialModel,
    requested_features::Vector{Symbol})
    for feat in requested_features
        if !(feat in MOI.features_available(m))
            error("Unsupported feature $feat")
        end
    end
end

MathOptInterface.features_available(m::OrdinalMultinomialModel) = [:Grad]

function MathOptInterface.eval_objective(
    m::OrdinalMultinomialModel,
    par::Vector)
    copyto!(m.α, 1, par, 1, m.J - 1)
    copyto!(m.β, 1, par, m.J, m.p)
    loglikelihood!(m, false, false)
end

function MathOptInterface.eval_objective_gradient(
    m::OrdinalMultinomialModel,
    grad::Vector,
    par::Vector)
    copyto!(m.α, 1, par, 1, m.J - 1)
    copyto!(m.β, 1, par, m.J, m.p)
    loglikelihood!(m, true, false)
    @views mul!(grad[1:m.J-1], m.dθdα, m.∇[1:m.J-1])
    copyto!(grad, m.J, m.∇, m.J, m.p)
end

function StatsModels.coeftable(mod::StatsModels.TableRegressionModel{T, S} 
    where {T <: OrdinalMultinomialModel, S <: Matrix} )
    ct = coeftable(mod.model)
    cfnames = [["intercept$i|$(i+1)" for i in 1:(mod.model.J - 1)]; coefnames(mod)]
    if length(ct.rownms) == length(cfnames)
        ct.rownms = cfnames
    end
    ct
end