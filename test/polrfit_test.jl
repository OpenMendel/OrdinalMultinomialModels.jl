module PolrfitTest

using Test, PolrModels, RDatasets

housing = dataset("MASS", "housing")
@info "Housing example for `polr` function in R package MASS"

@testset "logit link" begin
    for solver in [IpoptSolver(print_level=0), NLoptSolver(algorithm=:LD_SLSQP)]
        houseplr = polr(@formula(Sat ~ 0 + Infl + Type + Cont), housing,
            LogitLink(), solver; wts = housing[:Freq])
        @test nobs(houseplr) == 72
        @test isapprox(deviance(houseplr), 3479.149; rtol=1e-4)
        @test isapprox(coef(houseplr), [-0.4961353, 0.6907083, 0.5663937, 1.2888191, -0.5723501, -0.3661866, -1.0910149, 0.3602841]; rtol=1e-4)
        @test isapprox(stderror(houseplr), [0.12485, 0.12547, 0.104653, 0.127156, 0.119238, 0.155173, 0.151486, 0.095536]; rtol=1e-4)
    end
end

@testset "probit link" begin
    for solver in [IpoptSolver(print_level=0), NLoptSolver(algorithm=:LD_SLSQP)]
        houseplr = polr(@formula(Sat ~ 0 + Infl + Type + Cont), housing,
            ProbitLink(), solver; wts = housing[:Freq])
        @test nobs(houseplr) == 72
        @test isapprox(deviance(houseplr), 3479.689; rtol=1e-4)
        @test isapprox(coef(houseplr), [-0.29983, 0.42672, 0.34642, 0.78291, -0.34754, -0.21789, -0.66417, 0.22239]; rtol=1e-4)
        @test isapprox(stderror(houseplr), [0.07615, 0.07640, 0.064137, 0.076426, 0.072291, 0.094766, 0.091800, 0.058123]; rtol=1e-4)
    end
end

@testset "cloglog link" begin
    for solver in [IpoptSolver(print_level=0), NLoptSolver(algorithm=:LD_SLSQP)]
        houseplr = polr(@formula(Sat ~ 0 + Infl + Type + Cont), housing,
            CloglogLink(), solver; wts = housing[:Freq])
        @test nobs(houseplr) == 72
        @test isapprox(deviance(houseplr), 3484.053; rtol=1e-4)
        @test isapprox(coef(houseplr), [-0.79622, 0.05537, 0.38204, 0.91536, -0.40720, -0.28053, -0.74245, 0.20922]; rtol=1e-3)
        @test isapprox(stderror(houseplr), [0.08965, 0.08560, 0.070260, 0.092560, 0.086071, 0.111149, 0.101331, 0.065106]; rtol=1e-4)
    end
end


end
