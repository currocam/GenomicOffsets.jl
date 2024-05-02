"""
Geomtric genomic offset. TODO: Add reference
"""
struct GeometricGO{T<:Real}
  Cb::Matrix{T}
  K::Int
  λ::Real
  mx::Union{Nothing, Matrix{T}}
  sx::Union{Nothing, Matrix{T}}
end

"""
  fit(::Type{GeometricGO}, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}) where {T1<:Real,T2<:Real}

Fit the Geometric genomic offset model (that is, fit a LFMM) using the data `Y` and `X`.

# Arguments
- `Y::AbstractMatrix{T1}`: NxL matrix with individuals genotypes or allele frequencies.
- `X::AbstractMatrix{T}`: NxP matrix with environmental data.
- `K`: Number of latent factors to use when fitting the LFMM model. It can be determined using the elbow rule of thumb. If nothing is provided, the Tracy-Widom statistic to determine the number of latent factors (which might be overconservative).
- Optional arguments:
    - `center::Bool=true`: Center both the genotype and environmental data. Predicted environmental data is centered using the same mean as X. If false, data is assumed to be centered.
    - `scale::Bool=true`: Scale the environmental data. Predicted environmental data is scaled using the same standard deviation as X. If false, data is assumed to be scaled.
    - tw_threshold::Real=1e-3: Threshold to determine the number of latent factors if using the Tracy-Widom statistic. By default, it uses 1e-3.
    - λ::Real=1e-5: Regularization parameter for the LFMM.
# Returns
  - A GeometricGO model.
"""
function fit(::Type{GeometricGO}, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}, K::Union{Int, Nothing}=nothing; center=true,scale=true, tw_threshold::Real=1e-3, λ::Real=1e-5 ) where {T1<:Real,T2<:Real}
    mx = nothing; sx = nothing
    if center
        Y = Y .-  mean(Y, dims=1)
        mx = mean(X, dims=1)
        X = X .- mx
    end
    if scale
        sx = std(X, dims=1)
        X = X ./ sx
    end
    if isnothing(K)
        eigenvalues = eigvals(Y*Y'/(size(Y, 1)-1))
        _, pvalues = GenomicOffsets.TracyWidom(eigenvalues)
        K = findfirst(pvalues .> tw_threshold) - 1
    end
    Bt = RidgeLFMM(Y, X, K, λ).Bt
    Cb = Bt * Bt' / size(Bt, 2)
    return GeometricGO(Cb, K, λ, mx, sx)
end

"""
  genomic_offset(model::GeometricGO, X::AbstractMatrix{T}, Xpred::AbstractMatrix{T}) where T<:Real
  
  Compute the Geometric genomic offset.

# Arguments
- `model::GeometricGO`: A GeometricGO model.
- `X::AbstractMatrix{T}`: NxP matrix with environmental data. Data will be scaled and centered if the model was fitted with center and scale set to true.
- `Xpred::Matrix{T}`: NxP matrix with the predicted / altered environmental matrix. Data will be scaled and centered if the model was fitted with center and scale set to true.
# Returns
  - A Vector{Float64} of length N with the Geometric genomic offset values.
"""
function genomic_offset(model::GeometricGO, X::AbstractMatrix{T}, Xpred::AbstractMatrix{T}) where T<:Real
    if !isnothing(model.mx)
        X = X .- model.mx
        Xpred = Xpred .- model.mx
    end
    if !isnothing(model.sx)
        X = X ./ model.sx
        Xpred = Xpred ./ model.sx
    end
    offsets = zeros(size(X, 1))
    for i in eachindex(offsets)
        diff = X[i, :] - Xpred[i, :]
        offsets[i] = diff' * model.Cb * diff
    end
    return offsets
end