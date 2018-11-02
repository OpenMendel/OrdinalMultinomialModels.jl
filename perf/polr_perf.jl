module PerfTuning

using InteractiveUtils, Libdl, Random
if VERSION ≥ v"0.7.0"
    using Profile
end
using BenchmarkTools, PolrModels

Random.seed!(123)
n, p, J, q = 1_000_000, 5, 7, 5
Xtrue = randn(n, p + q)
X = Xtrue[:, 1:p] # null model matrix
Z = Xtrue[:, p+1:p+q] # covariates to be tested
βtrue = [ones(p); 0.01ones(q)]
θ = collect(1.0:J-1)
α = [θ[1]; [log(θ[i] - θ[i-1]) for i in 2:J-1]]
link = LogitLink() # LogitLink(), ProbitLink(), CauchitLink(), CloglogLink()

@info "generate responses by rpolr"
Y = rpolr(Xtrue, βtrue, θ, link)
# @code_warntype rpolyr(X[1, :], β, θ, :logit)
# @benchmark rpolyr(X, β, θ, :logit)
# Profile.clear_malloc_data()
# Profile.clear()
# @profile rpolyr(X, β, θ, :logit)
# Profile.print(format=:flat)

m = PolrModel(X, Y, link)
println()
@info "loglikelihood calculated by `polrfun!`"
@show polrfun!(m, false, false)
# println()
# @info "type inference on `polrfun`"
# @code_warntype polrfun!(m, true, true)
println()
@info "benchmark `polrfun`"
display(@benchmark polrfun!($m, true, false))
println(); println()
@info "profiling `polrfun`"
Profile.clear()
@profile polrfun!(m, true, true)
Profile.print(format=:flat)

@info "optimization"
# Gradient based: LD_LBFGS, :LD_MMA, :LD_SLSQP, :LD_CCSAQ, :LD_TNEWTON_PRECOND_RESTART, :LD_TNEWTON_PRECOND, :LD_TNEWTON_RESTART, :LD_VAR2, :LD_VAR1
# Gradient free: :LN_COBYLA
# solver = NLoptSolver(algorithm=:LD_SLSQP)
solver = IpoptSolver() # more stable but take a lot more iterations
# solver = IpoptSolver(mehrotra_algorithm="yes")
dd = polr(X, Y, link, solver)
display(coeftable(dd))

@info "score test"
ts = PolrScoreTest(dd, Z)
@show polrtest(ts)
@info "type inference"
@code_warntype polrtest(ts)
println()
@info "benchmark `polrtest`"
display(@benchmark polrtest($ts))
println(); println()
@info "profiling `polrtest`"
Profile.clear()
@profile polrtest(ts)
Profile.print(format=:flat)

end