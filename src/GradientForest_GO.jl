"""
    GradientForest_GO(Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}, Xpred::AbstractMatrix{T2}, ntrees::Int=100; center=true,scale=true) where {T1<:Real,T2<:Real}

Compute the Gradient Forest genomic offset. TODO: Add citation.

# Arguments
- `Y::AbstractMatrix{T1}`: NxL matrix with individuals genotypes or allele frequencies.
- `X::AbstractMatrix{T2}`: NxP matrix with environmental data.
- `Xpred::Matrix{T2}`: NxP matrix with the predicted / altered environmental matrix.
- `ntrees`: Number of trees each forest (one per allele) will have.
# Returns
  - A Vector{Float64} of length N with the RDA genomic offset values.

"""
function GradientForest_GO(Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}, Xpred::AbstractMatrix{T2}; ntrees::Int=100) where {T1<:Real,T2<:Real}
    gf = GenomicOffsets.GradientForest.gradient_forest(Y, X; ntrees = ntrees)
    N, P = size(X)
    distance = zeros(N)
    for p in 1:P
        F = gf.F[p]
        distance .+= abs2.(F.(X[:, p]) .- F.(Xpred[:, p]))
    end
    sqrt.(distance)
end
  