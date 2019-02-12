module PolrfitTest

using Test, OrdinalMultinomialModels, RDatasets

housing = dataset("MASS", "housing")
@info "Housing example for `polr` function in R package MASS"

@testset "logit link" begin
    for solver in [IpoptSolver(print_level=0), NLoptSolver(algorithm=:LD_SLSQP)]
        houseplr = polr(@formula(Sat ~ Infl + Type + Cont), housing,
            LogitLink(), solver; wts = housing[:Freq])
        @test nobs(houseplr) == 72
        @test isapprox(deviance(houseplr), 3479.149; rtol=1e-4)
        @test isapprox(coef(houseplr), [-0.4961353, 0.6907083, 0.5663937, 1.2888191, -0.5723501, -0.3661866, -1.0910149, 0.3602841]; rtol=1e-4)
        @test isapprox(stderror(houseplr), [0.124541, 0.125212, 0.104963, 0.126705, 0.118747, 0.156766, 0.151514, 0.0953575]; rtol=1e-4)
    end
end

@testset "probit link" begin
    for solver in [IpoptSolver(print_level=0), NLoptSolver(algorithm=:LD_SLSQP)]
        houseplr = polr(@formula(Sat ~ Infl + Type + Cont), housing,
            ProbitLink(), solver; wts = housing[:Freq])
        @test nobs(houseplr) == 72
        @test isapprox(deviance(houseplr), 3479.689; rtol=1e-4)
        @test isapprox(coef(houseplr), [-0.29983, 0.42672, 0.34642, 0.78291, -0.34754, -0.21789, -0.66417, 0.22239]; rtol=1e-4)
        @test isapprox(stderror(houseplr), [0.0761614, 0.0763991, 0.0641796, 0.0762645, 0.0722116, 0.0955741, 0.0919294, 0.0581214]; rtol=1e-4)
    end
end

@testset "cloglog link" begin
    for solver in [IpoptSolver(print_level=0), NLoptSolver(algorithm=:LD_SLSQP)]
        houseplr = polr(@formula(Sat ~ Infl + Type + Cont), housing,
            CloglogLink(), solver; wts = housing[:Freq])
        @test nobs(houseplr) == 72
        @test isapprox(deviance(houseplr), 3484.053; rtol=1e-4)
        @test isapprox(coef(houseplr), [-0.79622, 0.05537, 0.38204, 0.91536, -0.40720, -0.28053, -0.74245, 0.20922]; rtol=1e-3)
        @test isapprox(stderror(houseplr), [0.0904734, 0.0866657, 0.0701218, 0.0924964, 0.0861325, 0.112948, 0.102084, 0.0653755]; rtol=1e-4)
    end
end

end
