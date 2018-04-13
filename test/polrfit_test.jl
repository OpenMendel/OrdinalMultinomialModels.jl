module PolrfitTest

using Base.Test, BenchmarkTools, PolrModels

srand(123)

n, p, J = 500, 5, 7
X = randn(n, p)
β = ones(p)
θ = collect(1.0:J-1)
α = [θ[1]; [log(θ[i] - θ[i-1]) for i in 2:J-1]]
link = LogitLink() # LogitLink(), ProbitLink(), CauchitLink(), CloglogLink()

Y = rpolr(X, β, θ, link)
# @code_warntype rpolr(X[1, :], β, θ, :logit)
# @benchmark rpolyr(X, β, θ, :logit)
# Profile.clear_malloc_data()
# Profile.clear()
# @profile rpolyr(X, β, θ, :logit)
# Profile.print(format=:flat)

m = PolrModel(Y, X, link)
# polrfun!(m, true, true)
# @code_warntype polrfun!(m, true, true)
# @benchmark polrfun!(m, true, true)
# Profile.clear_malloc_data()
# Profile.clear()
# @profile polrfun!(m, true, true)
# Profile.print(format=:flat)

# solver = IpoptSolver()
# Gradient based: LD_LBFGS, :LD_MMA, :LD_SLSQP, :LD_CCSAQ, :LD_TNEWTON_PRECOND_RESTART, :LD_TNEWTON_PRECOND, :LD_TNEWTON_RESTART, :LD_VAR2, :LD_VAR1
# Gradient free: :LN_COBYLA
# solver = NLoptSolver(algorithm=:LD_SLSQP)
# solver = IpoptSolver() # more stable but take a lot more iterations
solver = IpoptSolver(mehrotra_algorithm="yes")
@time dd = polrmle(Y, X, link, solver)
# [[dd.θ; dd.β] stderr(dd)]
coeftable(dd)
end
