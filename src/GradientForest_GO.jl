using Interpolations
"""
Gradient genomic offset. TODO: Add reference
"""
struct GradientForestGO <: AbstractGO
  F::Vector{Interpolations.AbstractExtrapolation{Float64, 1}}
end

"""
  fit(::Type{GeometricGO}, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}) where {T1<:Real,T2<:Real}

Fit the Gradient forest genomic offset model (that is, fit a Gradient Forest) using the data `Y` and `X`.

# Arguments
- `Y::AbstractMatrix{T1}`: NxL matrix with individuals genotypes or allele frequencies.
- `X::AbstractMatrix{T}`: NxP matrix with environmental data.
- `ntrees`: Number of trees each forest (one per allele) will have.
# Returns
  - A GradientForestGO model.
"""
function fit(::Type{GradientForestGO}, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}; ntrees::Int=100) where {T1<:Real,T2<:Real}
    gf = GenomicOffsets.GradientForest.gradient_forest(Y, X; ntrees = ntrees)
    return GradientForestGO(gf.F)
end

"""
  genomic_offset(model::GradientForestGO, X::AbstractMatrix{T}, Xpred::AbstractMatrix{T}) where T<:Real
  
  Compute the Gradient Forest genomic offset.

# Arguments
- `model::GradientForestGO`: A GeometricGO model.
- `X::AbstractMatrix{T}`: NxP matrix with environmental data. Data will be scaled and centered if the model was fitted with center and scale set to true.
- `Xpred::Matrix{T}`: NxP matrix with the predicted / altered environmental matrix. Data will be scaled and centered if the model was fitted with center and scale set to true.
# Returns
  - A Vector{Float64} of length N with the Gradient forest genomic offset values.
"""
function genomic_offset(model::GradientForestGO, X::AbstractMatrix{T}, Xpred::AbstractMatrix{T}) where T<:Real
    N, P = size(X)
    distance = zeros(N)
    for p in 1:P
        F = model.F[p]
        distance .+= abs2.(F.(X[:, p]) .- F.(Xpred[:, p]))
    end
    sqrt.(distance)
end