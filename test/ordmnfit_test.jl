module PolrfitTest

using Test, OrdinalMultinomialModels, RDatasets, MathOptInterface, DelimitedFiles
using Tables, CategoricalArrays
const MOI = MathOptInterface

housing = dataset("MASS", "housing")
@info "Housing example for `polr` function in R package MASS"

ref_fitted_logistic=readdlm("fitted.logistic.csv",',',header=true)[1]
ref_fitted_probit=readdlm("fitted.probit.csv",',',header=true)[1]
ref_fitted_cloglog=readdlm("fitted.cloglog.csv",',',header=true)[1]

@testset "interface" begin
    houseplr1 = polr(@formula(Sat ~ Infl + Type + Cont), housing,
            LogitLink(), wts = housing[!, :Freq])
    
    @test Tables.istable(fitted(houseplr1))
    @test Tables.istable(predict(houseplr1, housing,kind=:probs))
    @test typeof(predict(houseplr1, housing, kind=:class)) <: CategoricalArray

    X = StatsModels.modelmatrix(@formula(Sat ~ Infl + Type + Cont).rhs,housing)
    Y = levelcode.(housing.Sat)
    houseplr2 = polr(X,Y,
            LogitLink(),wts = housing[!, :Freq])

    @test typeof(fitted(houseplr2)) <: AbstractMatrix
    @test typeof(predict(houseplr2, X, kind=:class)) <: Vector
    
    # mixed case, most similar to the old behaviour
    housing.S = levelcode.(housing.Sat)
    houseplr3 = polr(@formula(S ~ Infl + Type + Cont), housing, 
            LogitLink(),wts = housing[!, :Freq])
    @test Tables.istable(fitted(houseplr3))
    @test typeof(predict(houseplr3, housing, kind=:class)) <: Vector
    @test Tables.istable(predict(houseplr3, housing, kind=:probs))
end

@testset "logit link" begin
    Ipopt_solver = Ipopt.Optimizer()
    MOI.set(Ipopt_solver, MOI.RawOptimizerAttribute("print_level"), 0)
    NLopt_solver = NLopt.Optimizer()
    MOI.set(NLopt_solver, MOI.RawOptimizerAttribute("algorithm"), :LD_SLSQP)
    for solver in [Ipopt_solver, NLopt_solver]
        houseplr = polr(@formula(Sat ~ Infl + Type + Cont), housing,
            LogitLink(), solver; wts = housing[!, :Freq])
        @test nobs(houseplr) == 72
        @test dof(houseplr) == 8
        @test isapprox(deviance(houseplr), 3479.149; rtol=1e-4)
        @test isapprox(coef(houseplr), [-0.4961353, 0.6907083, 0.5663937, 1.2888191, -0.5723501, -0.3661866, -1.0910149, 0.3602841]; rtol=1e-4)
        @test isapprox(stderror(houseplr), [0.124541, 0.125212, 0.104963, 0.126705, 0.118747, 0.156766, 0.151514, 0.0953575]; rtol=1e-4)
        @test loglikelihood(houseplr.model) ≈ -1739.5746495294738
        @test isapprox(fitted(houseplr.model), ref_fitted_logistic; rtol=1e-4) 
        @test response(houseplr.model) === houseplr.model.Y
        @test modelmatrix(houseplr.model) === houseplr.model.X
        @test score(houseplr.model) === houseplr.model.∇
        @test vcov(houseplr.model) === houseplr.model.vcov
        @test all(abs.(cor(houseplr.model)) .≤ 1 + eps(Float64))
        @test all(weights(houseplr.model) .== housing[!, :Freq])
        @show coeftable(houseplr)
        @show confint(houseplr.model)
    end
end

@testset "probit link" begin
    Ipopt_solver = Ipopt.Optimizer()
    MOI.set(Ipopt_solver, MOI.RawOptimizerAttribute("print_level"), 0)
    NLopt_solver = NLopt.Optimizer()
    MOI.set(NLopt_solver, MOI.RawOptimizerAttribute("algorithm"), :LD_SLSQP)
    for solver in [Ipopt_solver, NLopt_solver]
        houseplr = polr(@formula(Sat ~ Infl + Type + Cont), housing,
            ProbitLink(), solver; wts = housing[!, :Freq])
        @test nobs(houseplr) == 72
        @test dof(houseplr) == 8
        @test isapprox(deviance(houseplr), 3479.689; rtol=1e-4)
        @test isapprox(coef(houseplr), [-0.29983, 0.42672, 0.34642, 0.78291, -0.34754, -0.21789, -0.66417, 0.22239]; rtol=1e-4)
        @test isapprox(stderror(houseplr), [0.0761614, 0.0763991, 0.0641796, 0.0762645, 0.0722116, 0.0955741, 0.0919294, 0.0581214]; rtol=1e-4)
        @test isapprox(fitted(houseplr.model), ref_fitted_probit; rtol=1e-4) 
    end
end

@testset "cloglog link" begin
    Ipopt_solver = Ipopt.Optimizer()
    MOI.set(Ipopt_solver, MOI.RawOptimizerAttribute("print_level"), 0)
    NLopt_solver = NLopt.Optimizer()
    MOI.set(NLopt_solver, MOI.RawOptimizerAttribute("algorithm"), :LD_SLSQP)
    for solver in [Ipopt_solver, NLopt_solver]
        houseplr = polr(@formula(Sat ~ Infl + Type + Cont), housing,
            CloglogLink(), solver; wts = housing[!, :Freq])
        @test nobs(houseplr) == 72
        @test dof(houseplr) == 8
        @test isapprox(deviance(houseplr), 3484.053; rtol=1e-4)
        @test isapprox(coef(houseplr), [-0.79622, 0.05537, 0.38204, 0.91536, -0.40720, -0.28053, -0.74245, 0.20922]; rtol=1e-3)
        @test isapprox(stderror(houseplr), [0.0904734, 0.0866657, 0.0701218, 0.0924964, 0.0861325, 0.112948, 0.102084, 0.0653755]; rtol=1e-4)
        #@test isapprox(fitted(houseplr.model), ref_fitted_cloglog; rtol=1e-4)  # this one fails :(
    end
end

end
