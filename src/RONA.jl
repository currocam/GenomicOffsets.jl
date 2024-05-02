using Random
"""
Risk of non-adaptnesss (RONA) genomic offset. First proposed by [Rellstab et al. (2016)]( https://doi.org/10.1111/mec.13889).
"""
struct RONA{T<:Real} <: AbstractGO
  B::Matrix{T}
end


"""
    fit(::Type{RONA}, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}) where {T1<:Real,T2<:Real}

Fit the RONA model (that is, solve the least squares problem) using the data `Y` and `X`.

# Arguments
- `Y::AbstractMatrix{T1}`: NxL matrix with individuals genotypes or allele frequencies.
- `X::AbstractMatrix{T}`: NxP matrix with environmental data.
# Returns
  - A RONA model.
"""
function fit(::Type{RONA}, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}) where {T1<:Real,T2<:Real}
  X = hcat(ones(size(X, 1)), X)
  B = X \ Y
  return RONA(B)
end

"""
    genomic_offset(model::RONA, X::AbstractMatrix{T}, Xpred::AbstractMatrix{T}) where T<:Real
  
  Compute the genomic offset for the RONA model.

# Arguments
- `model::RONA`: A RONA model.
- `X::AbstractMatrix{T}`: NxP matrix with environmental data.
- `Xpred::Matrix{T}`: NxP matrix with the predicted / altered environmental matrix.

# Returns
  - A Vector{Float64} of length N with the RONA genomic offset values.
"""
function genomic_offset(model::RONA, X::AbstractMatrix{T}, Xpred::AbstractMatrix{T}) where T<:Real
  X = hcat(ones(size(X, 1)), X)
  Xpred = hcat(ones(size(Xpred, 1)), Xpred)
  return sum(abs.(X * model.B - Xpred * model.B), dims=2) / size(model.B, 2)
end