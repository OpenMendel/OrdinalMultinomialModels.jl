var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "PolrModels.jl",
    "title": "PolrModels.jl",
    "category": "page",
    "text": ""
},

{
    "location": "#PolrModels.jl-1",
    "page": "PolrModels.jl",
    "title": "PolrModels.jl",
    "category": "section",
    "text": "PolrModels.jl provides Julia utilities to fit ordered multinomial models, including proportional odds model and ordered Probit model as special cases. "
},

{
    "location": "#Installation-1",
    "page": "PolrModels.jl",
    "title": "Installation",
    "category": "section",
    "text": "This package requires Julia v0.7 or later. The package has not yet been registered and must be installed using the repository location. Start julia and use the ] key to switch to the package manager REPL(v0.7) pkg> add https://github.com/OpenMendel/PolrModels.git# Machine info for results in this tutorial\nversioninfo()Julia Version 0.7.0\nCommit a4cb80f3ed (2018-08-08 06:46 UTC)\nPlatform Info:\n  OS: macOS (x86_64-apple-darwin14.5.0)\n  CPU: Intel(R) Core(TM) i7-6920HQ CPU @ 2.90GHz\n  WORD_SIZE: 64\n  LIBM: libopenlibm\n  LLVM: libLLVM-6.0.0 (ORCJIT, skylake)\nEnvironment:\n  JULIA_EDITOR = code# for use in this tutorial\nusing PolrModels, BenchmarkTools┌ Info: Recompiling stale cache file /Users/huazhou/.julia/compiled/v0.7/PolrModels/KD4ui.ji for PolrModels [b51826c6-2bc1-4eb0-a3af-c080fd29f1d0]\n└ @ Base loading.jl:1185"
},

{
    "location": "#Example-data-1",
    "page": "PolrModels.jl",
    "title": "Example data",
    "category": "section",
    "text": "housing is a data set from R package MASS. The outcome of interest is Sat (satisfication) that takes values Low, Medium, or High. Predictors include Infl (inflation, categorical), Type (housing type, categorical), and Cont (categorical). Freq codes number of observation for each combination of levels.using RDatasets\n\nhousing = dataset(\"MASS\", \"housing\")<table class=\"data-frame\"><thead><tr><th></th><th>Sat</th><th>Infl</th><th>Type</th><th>Cont</th><th>Freq</th></tr><tr><th></th><th>Categorical…</th><th>Categorical…</th><th>Categorical…</th><th>Categorical…</th><th>Int32</th></tr></thead><tbody><tr><th>1</th><td>Low</td><td>Low</td><td>Tower</td><td>Low</td><td>21</td></tr><tr><th>2</th><td>Medium</td><td>Low</td><td>Tower</td><td>Low</td><td>21</td></tr><tr><th>3</th><td>High</td><td>Low</td><td>Tower</td><td>Low</td><td>28</td></tr><tr><th>4</th><td>Low</td><td>Medium</td><td>Tower</td><td>Low</td><td>34</td></tr><tr><th>5</th><td>Medium</td><td>Medium</td><td>Tower</td><td>Low</td><td>22</td></tr><tr><th>6</th><td>High</td><td>Medium</td><td>Tower</td><td>Low</td><td>36</td></tr><tr><th>7</th><td>Low</td><td>High</td><td>Tower</td><td>Low</td><td>10</td></tr><tr><th>8</th><td>Medium</td><td>High</td><td>Tower</td><td>Low</td><td>11</td></tr><tr><th>9</th><td>High</td><td>High</td><td>Tower</td><td>Low</td><td>36</td></tr><tr><th>10</th><td>Low</td><td>Low</td><td>Apartment</td><td>Low</td><td>61</td></tr><tr><th>11</th><td>Medium</td><td>Low</td><td>Apartment</td><td>Low</td><td>23</td></tr><tr><th>12</th><td>High</td><td>Low</td><td>Apartment</td><td>Low</td><td>17</td></tr><tr><th>13</th><td>Low</td><td>Medium</td><td>Apartment</td><td>Low</td><td>43</td></tr><tr><th>14</th><td>Medium</td><td>Medium</td><td>Apartment</td><td>Low</td><td>35</td></tr><tr><th>15</th><td>High</td><td>Medium</td><td>Apartment</td><td>Low</td><td>40</td></tr><tr><th>16</th><td>Low</td><td>High</td><td>Apartment</td><td>Low</td><td>26</td></tr><tr><th>17</th><td>Medium</td><td>High</td><td>Apartment</td><td>Low</td><td>18</td></tr><tr><th>18</th><td>High</td><td>High</td><td>Apartment</td><td>Low</td><td>54</td></tr><tr><th>19</th><td>Low</td><td>Low</td><td>Atrium</td><td>Low</td><td>13</td></tr><tr><th>20</th><td>Medium</td><td>Low</td><td>Atrium</td><td>Low</td><td>9</td></tr><tr><th>21</th><td>High</td><td>Low</td><td>Atrium</td><td>Low</td><td>10</td></tr><tr><th>22</th><td>Low</td><td>Medium</td><td>Atrium</td><td>Low</td><td>8</td></tr><tr><th>23</th><td>Medium</td><td>Medium</td><td>Atrium</td><td>Low</td><td>8</td></tr><tr><th>24</th><td>High</td><td>Medium</td><td>Atrium</td><td>Low</td><td>12</td></tr><tr><th>25</th><td>Low</td><td>High</td><td>Atrium</td><td>Low</td><td>6</td></tr><tr><th>26</th><td>Medium</td><td>High</td><td>Atrium</td><td>Low</td><td>7</td></tr><tr><th>27</th><td>High</td><td>High</td><td>Atrium</td><td>Low</td><td>9</td></tr><tr><th>28</th><td>Low</td><td>Low</td><td>Terrace</td><td>Low</td><td>18</td></tr><tr><th>29</th><td>Medium</td><td>Low</td><td>Terrace</td><td>Low</td><td>6</td></tr><tr><th>30</th><td>High</td><td>Low</td><td>Terrace</td><td>Low</td><td>7</td></tr><tr><th>&vellip;</th><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td></tr></tbody></table>There are 72 unique combination of levels and the total number of observations is 1,681.size(housing, 1), sum(housing[:Freq])(72, 1681)"
},

{
    "location": "#PolrModels.polr",
    "page": "PolrModels.jl",
    "title": "PolrModels.polr",
    "category": "function",
    "text": "polr(formula, df, link, solver=NLoptSolver(algorithm=:LD_SLSQP, maxeval=4000))\npolr(X, y, link, solver=NLoptSolver(algorithm=:LD_SLSQP, maxeval=4000))\n\nFit ordered multinomial model by maximum likelihood estimation.\n\nPositional arguments\n\nformula::Formula: a model formula specifying responses and regressors.\ndf::DataFrame: a dataframe. Response variable should take integer values    starting from 1.\ny::Vector: integer vector taking values in 1,...,J.\nX::Matrix: n x p covariate matrix excluding intercept.\nlink::GLM.Link: LogitLink() (default), ProbitLink(), CauchitLink(), or CloglogLink().\nsolver: NLoptSolver() (default) or IpoptSolver().\n\nKeyword arguments\n\nwts::AbstractVector=similar(X, 0): observation weights.\n\nOutput\n\ndd:PolrModel: a PolrModel type.\n\n\n\n\n\n"
},

{
    "location": "#Syntax-1",
    "page": "PolrModels.jl",
    "title": "Syntax",
    "category": "section",
    "text": "polr is the main function of fitting ordered multinomial model. For documentation, type ?polr at Julia REPL.polr"
},

{
    "location": "#Fit-ordered-multinomial-models-1",
    "page": "PolrModels.jl",
    "title": "Fit ordered multinomial models",
    "category": "section",
    "text": ""
},

{
    "location": "#Proportional-odds-model-1",
    "page": "PolrModels.jl",
    "title": "Proportional odds model",
    "category": "section",
    "text": "To fit an ordered multinomial model using default link link=LogitLink(), i.e., proportional odds modelhouse_po = polr(@formula(Sat ~ 0 + Infl + Type + Cont), housing, \n    wts = housing[:Freq])StatsModels.DataFrameRegressionModel{PolrModel{Int64,Float64,LogitLink},Array{Float64,2}}\n\nFormula: Sat ~ Infl + Type + Cont\n\nCoefficients:\n      Estimate Std.Error  t value Pr(>|t|)\nθ1   -0.496141  0.124847 -3.97398   0.0002\nθ2    0.690706  0.125472  5.50486    <1e-6\nβ1    0.566392  0.104653   5.4121    <1e-6\nβ2     1.28881  0.127156  10.1357   <1e-14\nβ3   -0.572352  0.119238 -4.80008    <1e-5\nβ4   -0.366182  0.155173 -2.35983   0.0213\nβ5    -1.09101  0.151486 -7.20206    <1e-9\nβ6    0.360284 0.0955358   3.7712   0.0004note: Note\nIt is necessary to exclude intercept in formula because ordered multinomial model automatically includes intercept for proper modeling.Since there are J=3 categories in Sat, the fitted model has 2 intercept parameters  that satisfy theta_1 le theta_2. beta_1 beta_2 are  regression coefficients for Infl (3 levels), beta_3 beta_4 beta_5 for Type (4 levels), and beta_6 for Cont (2 levels). Deviance (-2 loglikelihood) of the fitted model isdeviance(house_po)3479.149299072586Estimated regression coefficients arecoef(house_po)8-element Array{Float64,1}:\n -0.49614061341454607\n  0.6907056439400201 \n  0.5663915864515744 \n  1.2888108250034271 \n -0.5723518911944651 \n -0.3661822718560424 \n -1.091011574737562  \n  0.3602844417528419with standard errorsstderror(house_po)8-element Array{Float64,1}:\n 0.12484725578571425\n 0.12547194322484329\n 0.10465278622775007\n 0.12715610075413825\n 0.1192380112188433 \n 0.15517334323365534\n 0.15148598827529147\n 0.09553578572175248"
},

{
    "location": "#Ordered-probit-model-1",
    "page": "PolrModels.jl",
    "title": "Ordered probit model",
    "category": "section",
    "text": "To fit an ordered probit model, we use link ProbitLink()house_op = polr(@formula(Sat ~ 0 + Infl + Type + Cont), housing, ProbitLink(), \n    wts = housing[:Freq])StatsModels.DataFrameRegressionModel{PolrModel{Int64,Float64,ProbitLink},Array{Float64,2}}\n\nFormula: Sat ~ Infl + Type + Cont\n\nCoefficients:\n      Estimate Std.Error  t value Pr(>|t|)\nθ1   -0.299829 0.0761537 -3.93716   0.0002\nθ2     0.42672 0.0764043  5.58502    <1e-6\nβ1    0.346423 0.0641371  5.40129    <1e-5\nβ2    0.782914 0.0764262  10.2441   <1e-14\nβ3   -0.347538 0.0722909 -4.80749    <1e-5\nβ4   -0.217889 0.0947661 -2.29923   0.0248\nβ5   -0.664175    0.0918 -7.23502    <1e-9\nβ6    0.222386 0.0581227  3.82616   0.0003deviance(house_op)3479.6888425652414"
},

{
    "location": "#Proportional-hazards-model-1",
    "page": "PolrModels.jl",
    "title": "Proportional hazards model",
    "category": "section",
    "text": "To fit a proportional hazards model, we use CloglogLink()house_ph = polr(@formula(Sat ~ 0 + Infl + Type + Cont), housing, CloglogLink(), \n    wts = housing[:Freq])StatsModels.DataFrameRegressionModel{PolrModel{Int64,Float64,CloglogLink},Array{Float64,2}}\n\nFormula: Sat ~ Infl + Type + Cont\n\nCoefficients:\n      Estimate Std.Error  t value Pr(>|t|)\nθ1   -0.796225 0.0896499 -8.88149   <1e-12\nθ2   0.0553535 0.0855972 0.646674   0.5202\nβ1    0.382045 0.0702599  5.43759    <1e-6\nβ2    0.915384 0.0925608  9.88954   <1e-13\nβ3   -0.407228 0.0860718 -4.73125    <1e-4\nβ4    -0.28056   0.11115 -2.52416   0.0141\nβ5   -0.742481  0.101331 -7.32728    <1e-9\nβ6    0.209235 0.0651057  3.21377   0.0021deviance(house_ph)3484.0531705626913From the deviances, we see the proportional odds model (logit link) has the best fit among all three models.deviance(house_po), deviance(house_op), deviance(house_ph)(3479.149299072586, 3479.6888425652414, 3484.0531705626913)"
},

{
    "location": "#Alternative-syntax-without-using-DataFrame-1",
    "page": "PolrModels.jl",
    "title": "Alternative syntax without using DataFrame",
    "category": "section",
    "text": "An alternative syntax is useful when it is inconvenient to use DataFramepolr(X, y, link, solver; wts)where y is the response vector and X is the n x p predictor matrix excluding intercept."
},

{
    "location": "#Optimization-algorithms-1",
    "page": "PolrModels.jl",
    "title": "Optimization algorithms",
    "category": "section",
    "text": "PolrModels.jl relies on nonlinear programming (NLP) optimization algorithms to find the maximum likelihood estimate (MLE). User can input any solver supported by the MathProgBase.jl package (see http://www.juliaopt.org) as the 4th argument of polr function. Common choices are:  Ipopt solver: IpoptSolver(print_level=0). See Ipopt.jl for numerous arugments to IpoptSolver. For example, setting print_level=5 is useful for diagnosis purpose.   \nNLopt package: NLoptSolver(algorithm=:LD_SLSQP), NLoptSolver(algorithm=:LD_LBFGS). See NLopt algorithms for all algorithms in NLopt.jl.When optimization fails, user can always try another algorithm.Use Ipopt (interior-point) solverpolr(@formula(Sat ~ 0 + Infl + Type + Cont), housing, LogitLink(), IpoptSolver(print_level=3); wts = housing[:Freq])******************************************************************************\nThis program contains Ipopt, a library for large-scale nonlinear optimization.\n Ipopt is released as open source code under the Eclipse Public License (EPL).\n         For more information visit http://projects.coin-or.org/Ipopt\n******************************************************************************\n\nTotal number of variables............................:        8\n                     variables with only lower bounds:        0\n                variables with lower and upper bounds:        0\n                     variables with only upper bounds:        0\nTotal number of equality constraints.................:        0\nTotal number of inequality constraints...............:        0\n        inequality constraints with only lower bounds:        0\n   inequality constraints with lower and upper bounds:        0\n        inequality constraints with only upper bounds:        0\n\n\nNumber of Iterations....: 38\n\n                                   (scaled)                 (unscaled)\nObjective...............:   2.0594068522767617e+02    1.7395746495294738e+03\nDual infeasibility......:   8.1478478823959181e-09    6.8824621073629131e-08\nConstraint violation....:   0.0000000000000000e+00    0.0000000000000000e+00\nComplementarity.........:   0.0000000000000000e+00    0.0000000000000000e+00\nOverall NLP error.......:   8.1478478823959181e-09    6.8824621073629131e-08\n\n\nNumber of objective function evaluations             = 91\nNumber of objective gradient evaluations             = 39\nNumber of equality constraint evaluations            = 0\nNumber of inequality constraint evaluations          = 0\nNumber of equality constraint Jacobian evaluations   = 0\nNumber of inequality constraint Jacobian evaluations = 0\nNumber of Lagrangian Hessian evaluations             = 0\nTotal CPU secs in IPOPT (w/o function evaluations)   =      0.057\nTotal CPU secs in NLP function evaluations           =      0.005\n\nEXIT: Optimal Solution Found.\n\n\n\n\n\nStatsModels.DataFrameRegressionModel{PolrModel{Int64,Float64,LogitLink},Array{Float64,2}}\n\nFormula: Sat ~ Infl + Type + Cont\n\nCoefficients:\n      Estimate Std.Error  t value Pr(>|t|)\nθ1   -0.496135  0.124847 -3.97394   0.0002\nθ2    0.690708  0.125472  5.50488    <1e-6\nβ1    0.566394  0.104653  5.41212    <1e-6\nβ2     1.28882  0.127156  10.1357   <1e-14\nβ3    -0.57235  0.119238 -4.80006    <1e-5\nβ4   -0.366186  0.155173 -2.35985   0.0213\nβ5    -1.09101  0.151486 -7.20208    <1e-9\nβ6    0.360284 0.0955358  3.77119   0.0004Use SLSQP (sequential quadratic programming) in NLopt.jl packagepolr(@formula(Sat ~ 0 + Infl + Type + Cont), housing, LogitLink(), NLoptSolver(algorithm=:LD_SLSQP); wts = housing[:Freq])StatsModels.DataFrameRegressionModel{PolrModel{Int64,Float64,LogitLink},Array{Float64,2}}\n\nFormula: Sat ~ Infl + Type + Cont\n\nCoefficients:\n      Estimate Std.Error  t value Pr(>|t|)\nθ1   -0.496141  0.124847 -3.97398   0.0002\nθ2    0.690706  0.125472  5.50486    <1e-6\nβ1    0.566392  0.104653   5.4121    <1e-6\nβ2     1.28881  0.127156  10.1357   <1e-14\nβ3   -0.572352  0.119238 -4.80008    <1e-5\nβ4   -0.366182  0.155173 -2.35983   0.0213\nβ5    -1.09101  0.151486 -7.20206    <1e-9\nβ6    0.360284 0.0955358   3.7712   0.0004Use LBFGS (quasi-Newton algorithm) in NLopt.jl packagepolr(@formula(Sat ~ 0 + Infl + Type + Cont), housing, LogitLink(), NLoptSolver(algorithm=:LD_LBFGS); wts = housing[:Freq])StatsModels.DataFrameRegressionModel{PolrModel{Int64,Float64,LogitLink},Array{Float64,2}}\n\nFormula: Sat ~ Infl + Type + Cont\n\nCoefficients:\n      Estimate Std.Error  t value Pr(>|t|)\nθ1   -0.496111  0.124847 -3.97375   0.0002\nθ2    0.690732  0.125472  5.50507    <1e-6\nβ1    0.566394  0.104653  5.41212    <1e-6\nβ2     1.28882  0.127156  10.1357   <1e-14\nβ3   -0.572352  0.119238 -4.80008    <1e-5\nβ4    -0.36616  0.155173 -2.35968   0.0214\nβ5    -1.09102  0.151486 -7.20212    <1e-9\nβ6    0.360319 0.0955359  3.77156   0.0004"
},

{
    "location": "#Likelihood-ratio-test-(LRT)-1",
    "page": "PolrModels.jl",
    "title": "Likelihood ratio test (LRT)",
    "category": "section",
    "text": "polr function calculates the Wald test (or t-test) p-value for each predictor in the model. To carry out the potentially more powerful likelihood ratio test (LRT), we need to fill the null and alternative models separately.Step 1: Fit the null model with only Infl and Type factors.house_null = polr(@formula(Sat ~ 0 + Infl + Type), housing; wts = housing[:Freq])StatsModels.DataFrameRegressionModel{PolrModel{Int64,Float64,LogitLink},Array{Float64,2}}\n\nFormula: Sat ~ Infl + Type\n\nCoefficients:\n      Estimate Std.Error  t value Pr(>|t|)\nθ1   -0.672949  0.115689  -5.8169    <1e-6\nθ2    0.505629  0.115203  4.38903    <1e-4\nβ1    0.548392  0.104308  5.25743    <1e-5\nβ2      1.2373  0.125978  9.82155   <1e-13\nβ3   -0.521441  0.118215 -4.41094    <1e-4\nβ4   -0.289347  0.153387 -1.88639   0.0637\nβ5    -1.01404  0.149562 -6.78008    <1e-8Step 2: To test significance of the Cont variable, we use polrtest function. The first argument is the fitted null model, the second argument is the predictor vector to be tested# last column of model matrix is coding for Cont (2-level factor)\ncont = modelmatrix(house_po.model)[:, end]\n# calculate p-value\npolrtest(house_null, cont; test=:LRT)0.000155351855453278"
},

{
    "location": "#Score-test-1",
    "page": "PolrModels.jl",
    "title": "Score test",
    "category": "section",
    "text": "User can perform score test using the polrtest function too. Score test has the advantage that, when testing a huge number of predictors such as in genomewide association studies (GWAS), one only needs to fit the null model once and then testing each predictor is cheap. Both Wald and likelihood ratio test (LRT) need to fit a separate alternative model for each predictor being tested.Step 1: Fit the null model with only Infl and Type factors.house_null = polr(@formula(Sat ~ 0 + Infl + Type), housing; wts = housing[:Freq])StatsModels.DataFrameRegressionModel{PolrModel{Int64,Float64,LogitLink},Array{Float64,2}}\n\nFormula: Sat ~ Infl + Type\n\nCoefficients:\n      Estimate Std.Error  t value Pr(>|t|)\nθ1   -0.672949  0.115689  -5.8169    <1e-6\nθ2    0.505629  0.115203  4.38903    <1e-4\nβ1    0.548392  0.104308  5.25743    <1e-5\nβ2      1.2373  0.125978  9.82155   <1e-13\nβ3   -0.521441  0.118215 -4.41094    <1e-4\nβ4   -0.289347  0.153387 -1.88639   0.0637\nβ5    -1.01404  0.149562 -6.78008    <1e-8Step 2: To test significance of the Cont variable, we use polrtest function. The first argument is the fitted null model, the second argument is the predictor vector to be tested# last column of model matrix is coding for Cont (2-level factor)\ncont = modelmatrix(house_po.model)[:, end]\n# calculate p-value\npolrtest(house_null, cont; test=:score)0.00015872510523744706Step 3: Now suppose we want to test significance of another predictor, z1. We just need to call polrtest with z1 and the same fiited null model. No model fitting is needed.For demonstration purpose, we generate z1 randomly. The score test p-value of z1 is, not suprisingly, large.z1 = randn(nobs(house_null))\npolrtest(house_null, z1)0.33703682134036256Step 4: We can also test a set of precitors or a factor.z3 = randn(nobs(house_null), 3)\npolrtest(house_null, z3)1.585206970323263e-14"
},

]}
