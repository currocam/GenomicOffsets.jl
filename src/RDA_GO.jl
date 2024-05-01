using LinearAlgebra
using StatsBase
"""
    RDA_GO(Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}, Xstar::AbstractMatrix{T2}) where {T1<:Real,T2<:Real}

Compute the redundancy Analysis (RDA) genomic offset. 

# Arguments
- `Y::AbstractMatrix{T1}`: NxL matrix with individuals genotypes or allele frequencies.
- `X::AbstractMatrix{T2}`: NxP matrix with environmental data.
- `Xstar::Matrix{T2}`: NxP matrix with the predicted / altered environmental matrix.
# Returns
  - A Vector{Float64} of length N with the RDA genomic offset values.

"""
function RDA_GO(Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}, Xpred::AbstractMatrix{T2}) where {T1<:Real,T2<:Real}
    # Add intercept to X and Xstar
    X = hcat(ones(size(X, 1)), X)
    Xpred = hcat(ones(size(Xpred, 1)), Xpred)
    # Fit linear model
    B = X \ Y
    # Project predicted frequencies
    loadings = eigen(cov(X * B)).vectors
    projX = (X * B) * loadings
    projXstar = (Xpred * B) * loadings
    return sum((projX - projXstar) .^ 2, dims=2) / size(Y, 2)
  end
  