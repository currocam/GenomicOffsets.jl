"""
Geometric genomic offset. TODO: Add reference
"""
struct GeometricGO{T<:Real} <: AbstractGO
    model::LFMM{T}
    K::Int
    λ::Real
    mx::Union{Nothing,Matrix{T}}
    sx::Union{Nothing,Matrix{T}}
end

"""
    fit(::Type{GeometricGO}, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}) where {T1<:Real,T2<:Real}

Fit the Geometric genomic offset model (that is, fit a LFMM) using the data `Y` and `X`.

# Arguments
- `Y::AbstractMatrix{T1}`: NxL matrix with individuals genotypes or allele frequencies.
- `X::AbstractMatrix{T}`: NxP matrix with environmental data.
- `K`: Number of latent factors to use when fitting th  e LFMM model. It can be determined using the elbow rule of thumb. If nothing is provided, the Tracy-Widom statistic to determine the number of latent factors (which might be overconservative).
- Optional arguments:
    - `center::Bool=true`: Center both the genotype and environmental data. Predicted environmental data is centered using the same mean as X. If false, data is assumed to be centered.
    - `scale::Bool=true`: Scale the environmental data. Predicted environmental data is scaled using the same standard deviation as X. If false, data is assumed to be scaled.
    - tw_threshold::Real=1e-3: Threshold to determine the number of latent factors if using the Tracy-Widom statistic. By default, it uses 1e-3.
    - λ::Real=1e-5: Regularization parameter for the LFMM.
# Returns
  - A GeometricGO model.
"""
function fit(::Type{GeometricGO}, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2},
             K::Union{Int,Nothing}=nothing; center=true, scale=true,
             tw_threshold::Real=1e-3, λ::Real=1e-5) where {T1<:Real,T2<:Real}
    mx = nothing
    sx = nothing
    if center
        Y = Y .- mean(Y; dims=1)
        mx = mean(X; dims=1)
        X = X .- mx
    end
    if scale
        sx = std(X; dims=1)
        X = X ./ sx
    end
    if isnothing(K)
        eigenvalues = eigvals(Y * Y' / (size(Y, 1) - 1))
        _, pvalues = GenomicOffsets.TracyWidom(eigenvalues)
        K = findfirst(pvalues .> tw_threshold) - 1
    end
    model = RidgeLFMM(Y, X, K, λ)
    return GeometricGO(model, K, λ, mx, sx)
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
function genomic_offset(model::GeometricGO, X::AbstractMatrix{T},
                        Xpred::AbstractMatrix{T}) where {T<:Real}
    if !isnothing(model.mx)
        X = X .- model.mx
        Xpred = Xpred .- model.mx
    end
    if !isnothing(model.sx)
        X = X ./ model.sx
        Xpred = Xpred ./ model.sx
    end
    offsets = zeros(size(X, 1))
    Bt = model.model.Bt
    Cb = Bt * Bt' / size(Bt, 2)
    for i in eachindex(offsets)
        diff = X[i, :] - Xpred[i, :]
        offsets[i] = diff' * Cb * diff
    end
    return offsets
end

"""
    genomic_offset(model::GeometricGO, X::AbstractMatrix{T}, Xpred::AbstractMatrix{T}) where T<:Real
  
Compute the Geometric genomic offset.

# Arguments
- `model::GeometricGO`: A GeometricGO model.
- `X::AbstractMatrix{T}`: NxP matrix with environmental data. Data will be scaled and centered if the model was fitted with center and scale set to true.
- `Xpred::Matrix{T}`: NxP matrix with the predicted / altered environmental matrix. Data will be scaled and centered if the model was fitted with center and scale set to true.
- `candidates::AbstractVector{Integer}`: Indices of the candidate loci.
# Returns
  - A Vector{Float64} of length N with the Geometric genomic offset values.
"""
function genomic_offset(model::GeometricGO{T}, X::AbstractMatrix{T},
                        Xpred::AbstractMatrix{T},
                        candidates::AbstractVector{T2}) where {T<:Real,T2<:Integer}
    if length(candidates) == 0
        return zeros(size(X, 1))
    end
    if !isnothing(model.mx)
        X = X .- model.mx
        Xpred = Xpred .- model.mx
    end
    if !isnothing(model.sx)
        X = X ./ model.sx
        Xpred = Xpred ./ model.sx
    end
    offsets = zeros(size(X, 1))
    Bt = model.model.Bt[:, candidates]
    Cb = Bt * Bt' / length(candidates)
    for i in eachindex(offsets)
        diff = X[i, :] - Xpred[i, :]
        offsets[i] = diff' * Cb * diff
    end
    return offsets
end

"""
    bootstrap_with_candidates(::Type{GeometricGO}, rng::Random.AbstractRNG,Y::AbstractMatrix{T1}, X::AbstractMatrix{T2},Xpred::AbstractMatrix{T2}, nboot::Int=500; candidates_threshold::Real=0.05, genomic_control::Bool=true, tw_threshold::Real=0.001) where {T1<:Real,T2<:Real}

Compute the Geometric genomic offset using a bootstrap approach. For every, bootstrap iteration, the model is fitted using a random subset of the columns of `Y`. The Tracy-Widom test is used to estimate the number of latent factors. The F-test is used to select putatively adaptative loci. The genomic offset is computed for the selected loci.

# Arguments
- `rng::Random.AbstractRNG`: Random number generator. If not provided, the global RNG is used.
- `Y::AbstractMatrix{T1}`: NxL matrix with individuals genotypes or allele frequencies.
- `X::AbstractMatrix{T2}`: NxP matrix with environmental data.
- `Xpred::AbstractMatrix{T2}`: NxP matrix with the predicted / altered environmental matrix.
- `nboot::Int=500`: Number of bootstrap iterations.
- Optional arguments:
  - `candidates_threshold::Real=0.05`: Threshold for the F-test. Better to be permisive.
  - `genomic_control::Bool=true`: Apply genomic control to the F-test.
  - `tw_threshold::Real=0.001`: Threshold for the Tracy-Widom test. Better to be conservative (as the number of latent factors is often overconservative).

# Returns
  - A matrix of size NxNboot with the genomic offset values. If no candidate loci are found, the row is filled with zeros.
"""
function bootstrap_with_candidates(::Type{GeometricGO}, rng::Random.AbstractRNG,
                                   Y::AbstractMatrix{T1}, X::AbstractMatrix{T2},
                                   Xpred::AbstractMatrix{T2}, nboot::Int=500;
                                   candidates_threshold::Real=0.05,
                                   genomic_control::Bool=true,
                                   tw_threshold::Real=0.001) where {T1<:Real,T2<:Real}
    Y = Y .- mean(Y; dims=1)
    mx = mean(X; dims=1)
    X = X .- mx
    Xpred = Xpred .- mx
    sx = std(X; dims=1)
    X = X ./ sx
    Xpred = X ./ sx
    _, L = size(Y)
    offsets = zeros(size(Y, 1), nboot)
    shared_seed = rand(rng, UInt)
    Threads.@threads for i in 1:nboot
        _rng = Random.seed!(copy(rng), shared_seed + i)
        Yboot = Y[:,sample(_rng, 1:L, L; replace=true)]
        eigenvalues = eigvals(Yboot * Yboot' / (size(Yboot, 1) - 1))
        _, pvalues = TracyWidom(eigenvalues)
        K = max(findfirst(pvalues .> tw_threshold) - 1, 1)
        model = GenomicOffsets.fit(GeometricGO, Yboot, X, K; center=false, scale=false)
        pvalues = LFMM_Ftest(model.model, Yboot, X; genomic_control=genomic_control,
                             center=false)
        qvalues = MultipleTesting.adjust(pvalues, BenjaminiHochberg())
        candidates = findall(qvalues .< candidates_threshold)
        offsets[:, i] = genomic_offset(model, X, Xpred, candidates)
    end
    return offsets
end
function bootstrap_with_candidates(::Type{GeometricGO}, Y::AbstractMatrix{T1},
                                   X::AbstractMatrix{T2}, Xpred::AbstractMatrix{T2},
                                   nboot::Int=500; candidates_threshold::Real=0.05,
                                   genomic_control::Bool=true,
                                   tw_threshold::Real=0.001) where {T1<:Real,T2<:Real}
    return bootstrap_with_candidates(GeometricGO, Random.GLOBAL_RNG, Y, X, Xpred, nboot;
                                     candidates_threshold=candidates_threshold,
                                     genomic_control=genomic_control,
                                     tw_threshold=tw_threshold)
end
