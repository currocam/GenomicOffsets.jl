"""
Geometric_GO(Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}, Xpred::AbstractMatrix{T2}) where {T1<:Real,T2<:Real}

Compute the geometric genomic offset as described by Gain et al. TODO: Add citation.

# Arguments
- `Y::AbstractMatrix{T1}`: NxL matrix with individuals genotypes or allele frequencies.
- `X::AbstractMatrix{T2}`: NxP matrix with environmental data.
- `Xpred::Matrix{T2}`: NxP matrix with the predicted / altered environmental matrix.
- `K`: Number of latent factors to use when fitting the LFMM model. It can be determined using the elbow rule of thumb. If nothing is provided, the Tracy-Widom statistic to determine the number of latent factors (which might be overconservative).
- Optional arguments:
    - `center::Bool=true`: Center both the genotype and environmental data. Predicted environmental data is centered using the same mean as X. If false, data is assumed to be centered.
    - `scale::Bool=true`: Scale the environmental data. Predicted environmental data is scaled using the same standard deviation as X. If false, data is assumed to be scaled.
    - tw_threshold::Real=1e-3: Threshold to determine the number of latent factors if using the Tracy-Widom statistic. By default, it uses 1e-3.
    - λ::Real=1e-5: Regularization parameter for the LFMM.
# Returns
  - A Vector{Float64} of length N with the RDA genomic offset values.

"""
function Geometric_GO(Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}, Xpred::AbstractMatrix{T2}, K::Union{Int, Nothing}=nothing; center=true,scale=true, tw_threshold::Real=1e-3, λ::Real=1e-5) where {T1<:Real,T2<:Real}
    if center
        Y = Y .- mean(Y, dims=1)
        mx = mean(X, dims=1)
        X = X .- mx
        Xpred = Xpred .- mx
    end
    if scale
        sx = std(X, dims=1)
        X = X ./ sx
        Xpred = Xpred ./ sx
    end
    if isnothing(K)
        eigenvalues = eigvals(Y*Y'/(size(Y, 1)-1))
        _, pvalues = GenomicOffsets.TracyWidom(eigenvalues)
        K = findfirst(pvalues .> tw_threshold) - 1
    end
    Bt = RidgeLFMM(Y, X, K, λ).Bt
    Cb = Bt * Bt' / size(Bt, 2)
    offsets = zeros(size(X, 1))
    for i in eachindex(offsets)
        diff = X[i, :] - Xpred[i, :]
        offsets[i] = diff' * Cb * diff
    end
    return offsets
end
  