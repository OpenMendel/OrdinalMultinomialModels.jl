"""
    rpolyr(x, β, θ)

Generate a random integer `Y` such that `P(Y≤j) = g^{-1}(θ[j]-x'β)`.
`θ` has to be monotone `θ[1] ≤ ... ≤ θ[J-1]`.
"""
function rpolr(
    x::AbstractVector,
    β::AbstractVector,
    θ::AbstractVector,
    link::GLM.Link
    )
    # number of categories
    J = length(θ) + 1
    # check monotonicity of θj
    issorted(θ) || throw(ArgumentError("θ[j] should be nondecreasing."))
    # generate category according to cumulative probability
    iprod = dot(x, β)
    q = rand()
    for j in 1:J-1
        ηj = θ[j] - iprod
        cumprobj = GLM.linkinv(link, ηj)
        q ≤ cumprobj && (return j)
    end
    return J
end

"""
    rpolr(X, β, θ)

Generate a vector of random integers `Y` such that
`P(Y[i]≤j) = g^{-1}(θ[j]-X[:,j]'β)`.
`θ` has to be monotone `θ[1] ≤ ... ≤ θ[J-1]`.
"""
function rpolr(
    X::AbstractMatrix,
    β::AbstractVector,
    θ::AbstractVector,
    link::GLM.Link
    )
    @views Y = [rpolr(X[i, :], β, θ, link) for i in 1:size(X, 1)]
end
