using LinearAlgebra
using StatsBase

"""
Redundancy Analysis (RDA) genomic offset. TODO: Add reference
"""
struct RDAGO{T<:Real}
  W::Matrix{T}
end

"""
  fit(::Type{RDAGO}, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}) where {T1<:Real,T2<:Real}

Fit the RDA model (that is, apply Redundancy Analysis) using the data `Y` and `X`.

# Arguments
- `Y::AbstractMatrix{T1}`: NxL matrix with individuals genotypes or allele frequencies.
- `X::AbstractMatrix{T}`: NxP matrix with environmental data.
# Returns
  - A RDAGO model.
"""
function fit(::Type{RDAGO}, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}) where {T1<:Real,T2<:Real}
  X = hcat(ones(size(X, 1)), X)
  # Fit linear model
  B = X \ Y
  # Project predicted frequencies
  loadings = eigen(cov(X * B)).vectors
  W = B * loadings
  return RDAGO(W)
end

"""
  genomic_offset(model::RDAGO, X::AbstractMatrix{T}, Xpred::AbstractMatrix{T}) where T<:Real
  
  Compute the RDA genomic offset.

# Arguments
- `model::RDAGO`: A RDAGO model.
- `X::AbstractMatrix{T}`: NxP matrix with environmental data.
- `Xpred::Matrix{T}`: NxP matrix with the predicted / altered environmental matrix.

# Returns
  - A Vector{Float64} of length N with the RDAGO genomic offset values.
"""
function genomic_offset(model::RDAGO, X::AbstractMatrix{T}, Xpred::AbstractMatrix{T}) where T<:Real
  projX =  hcat(ones(size(X, 1)), X) * model.W
  projXpred = hcat(ones(size(Xpred, 1)), Xpred) * model.W
  return sum((projX - projXpred) .^ 2, dims=2) / size(model.W, 2)
end