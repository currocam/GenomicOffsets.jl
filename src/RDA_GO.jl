using LinearAlgebra
using StatsBase
using MultivariateStats

"""
Redundancy Analysis (RDA) genomic offset. TODO: Add reference
"""
struct RDAGO{T<:Real} <: AbstractGO
    B::Matrix{T}
    PCA::PCA{T}
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
function fit(::Type{RDAGO}, Y::AbstractMatrix{T1},
             X::AbstractMatrix{T2}, ) where {T1<:Real,T2<:Real}
    X = hcat(ones(size(X, 1)), X)
    # Fit linear model
    B = X \ Y
    # Project predicted frequencies
    M = MultivariateStats.fit(PCA, (X * B)')
    return RDAGO(B, M)
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
function genomic_offset(model::RDAGO, X::AbstractMatrix{T},
                        Xpred::AbstractMatrix{T}) where {T<:Real}
    projX = predict(model.PCA, (hcat(ones(size(X, 1)), X) * model.B)')'
    projXpred = predict(model.PCA, (hcat(ones(size(Xpred, 1)), Xpred) * model.B)')'
    return sum((projX - projXpred) .^ 2; dims=2) / size(model.B, 2)
end

"""
    bootstrap_with_candidates(::Type{RDAGO}, rng::Random.AbstractRNG, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}, Xpred::AbstractMatrix{T2}, nboot::Int=500; candidates_threshold::Real=0.05, genomic_control::Bool=true, tw_threshold::Real=0.001) where {T1<:Real, T2<:Real}  
  
Compute the RDA genomic offset using a bootstrap approach. For every, bootstrap iteration, the model is fitted using a random subset of the columns of `Y`. The Tracy-Widom test is used to estimate the number of latent factors. The F-test is used to select putatively adaptative loci. The genomic offset is computed for the selected loci.

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
function bootstrap_with_candidates(::Type{RDAGO}, rng::Random.AbstractRNG,
                                   Y::AbstractMatrix{T1}, X::AbstractMatrix{T2},
                                   Xpred::AbstractMatrix{T2}, nboot::Int=500;
                                   candidates_threshold::Real=0.05,
                                   genomic_control::Bool=true,
                                   tw_threshold::Real=0.001) where {T1<:Real,T2<:Real}
    _, L = size(Y)
    # Center and scale data for LFMM
    mx = mean(X; dims=1)
    sdx = std(X; dims=1)
    Xscaled = (X .- mx) ./ sdx
    Yscaled = Y .- mean(Y; dims=1)
    # Allocate memory & initialize seed
    offsets = zeros(size(Y, 1), nboot)
    shared_seed = rand(rng, UInt)
    Threads.@threads for i in 1:nboot
        _rng = Random.seed!(copy(rng), shared_seed + i)
        sampled = sample(_rng, 1:L, L; replace=true)
        Yboot = Yscaled[:, sampled]
        # Fit LFMM for GEA
        lfmm = fit(GeometricGO, Yboot, Xscaled; center=false, scale=false,
                   tw_threshold=tw_threshold).model
        pvalues = LFMM_Ftest(lfmm, Yboot, Xscaled; genomic_control=genomic_control,
                             center=false)
        qvalues = adjust(pvalues, BenjaminiHochberg())
        candidates = findall(qvalues .< candidates_threshold)
        if length(candidates) > 0
            model = fit(RDAGO, Yboot[:, candidates], X)
            offsets[:, i] = genomic_offset(model, X, Xpred)
        end
    end
    return offsets
end

function bootstrap_with_candidates(::Type{RDAGO}, Y::AbstractMatrix{T1},
                                   X::AbstractMatrix{T2}, Xpred::AbstractMatrix{T2},
                                   nboot::Int=500; candidates_threshold::Real=0.05,
                                   genomic_control::Bool=true,
                                   tw_threshold::Real=0.001) where {T1<:Real,T2<:Real}
    return bootstrap_with_candidates(RDAGO, Random.GLOBAL_RNG, Y, X, Xpred, nboot;
                                     candidates_threshold=candidates_threshold,
                                     genomic_control=genomic_control,
                                     tw_threshold=tw_threshold)
end
