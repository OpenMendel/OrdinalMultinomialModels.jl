module PolrtestTest

using Test, PolrModels, Random

Random.seed!(123)

n, p, J, q = 1000, 5, 7, 5
Xtrue = randn(n, p + q)
X = Xtrue[:, 1:p] # null model matrix
Z = Xtrue[:, p+1:p+q] # covariates to be tested
βtrue = [ones(p); 0.1ones(q)]
θ = collect(1.0:J-1)
α = [θ[1]; [log(θ[i] - θ[i-1]) for i in 2:J-1]]
link = LogitLink() # LogitLink(), ProbitLink(), CauchitLink(), CloglogLink()

Y = rpolr(Xtrue, βtrue, θ, link)

# Gradient based: LD_LBFGS, :LD_MMA, :LD_SLSQP, :LD_CCSAQ, :LD_TNEWTON_PRECOND_RESTART, :LD_TNEWTON_PRECOND, :LD_TNEWTON_RESTART, :LD_VAR2, :LD_VAR1
# Gradient free: :LN_COBYLA
# solver = NLoptSolver(algorithm=:LD_LBFGS)
solver = IpoptSolver() # more stable but take a lot more iterations
@time dd = polr(X, Y, link, solver)
@show dd.θ, dd.β
@show stderror(dd)

# testing
@show polrtest(dd, Z; test = :score)
@show polrtest(dd, Z; test = :LRT)

end
