using BenchmarkTools, PolrModels

srand(123)
n, p, J, q = 1000000, 5, 7, 3
Xtrue = randn(n, p + q)
X = Xtrue[:, 1:p] # null model matrix
Z = Xtrue[:, p+1:p+q] # covariates to be tested
βtrue = [ones(p); 0.1ones(q)]
θ = collect(1.0:J-1)
α = [θ[1]; [log(θ[i] - θ[i-1]) for i in 2:J-1]]
link = LogitLink() # LogitLink(), ProbitLink(), CauchitLink(), CloglogLink()

Y = rpolr(Xtrue, βtrue, θ, link)
# @code_warntype rpolyr(X[1, :], β, θ, :logit)
# @benchmark rpolyr(X, β, θ, :logit)
# Profile.clear_malloc_data()
# Profile.clear()
# @profile rpolyr(X, β, θ, :logit)
# Profile.print(format=:flat)

# m = PolrModel(X, Y, link)
# polrfun!(m, false, false)
# @code_warntype polrfun!(m, true, true)
# @benchmark polrfun!(m, true, true)
# Profile.clear_malloc_data()
# Profile.clear()
# @profile polyrfun!(m, true, true)
# Profile.print(format=:flat)

# solver = IpoptSolver()
# Gradient based: LD_LBFGS, :LD_MMA, :LD_SLSQP, :LD_CCSAQ, :LD_TNEWTON_PRECOND_RESTART, :LD_TNEWTON_PRECOND, :LD_TNEWTON_RESTART, :LD_VAR2, :LD_VAR1
# Gradient free: :LN_COBYLA
# solver = NLoptSolver(algorithm=:LD_SLSQP)
# solver = IpoptSolver() # more stable but take a lot more iterations
solver = IpoptSolver(mehrotra_algorithm="yes")
@time dd = polr(X, Y, link, solver)
coeftable(dd)

# testing
# polrtest(dd, Z)
ts = PolrScoreTest(dd, Z)
# @code_warntype polrtest(PolrScoreTest(dd, Z))
@benchmark polrtest(ts)
