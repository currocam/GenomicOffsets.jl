"""
    RONA(Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}, Xpred::AbstractMatrix{T2}) where {T1<:Real,T2<:Real}

Compute the RONA (Risk of non-adaptnesss) genomic offset. First proposed by [Rellstab et al. (2016)]( https://doi.org/10.1111/mec.13889). 

# Arguments
- `Y::AbstractMatrix{T1}`: NxL matrix with individuals genotypes or allele frequencies.
- `X::AbstractMatrix{T2}`: NxP matrix with environmental data.
- `Xpred::Matrix{T2}`: NxP matrix with the predicted / altered environmental matrix.
# Returns
  - A Vector{Float64} of length N with the RONA genomic offset values.

"""
function RONA(Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}, Xpred::AbstractMatrix{T2}) where {T1<:Real,T2<:Real}
    X = hcat(ones(size(X, 1)), X)
    Xpred = hcat(ones(size(Xpred, 1)), Xpred)
    # Solve LSL
    B = X \ Y
    return sum(abs.(X * B - Xpred * B), dims=2) / size(Y, 2)
end