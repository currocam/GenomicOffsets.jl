using Random, MultipleTesting, StatsBase
"""
Risk of non-adaptnesss (RONA) genomic offset. First proposed by Rellstab et al. (2016).
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
function fit(::Type{RONA}, Y::AbstractMatrix{T1},
             X::AbstractMatrix{T2}) where {T1<:Real,T2<:Real}
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
function genomic_offset(model::RONA, X::AbstractMatrix{T},
                        Xpred::AbstractMatrix{T}) where {T<:Real}
    X = hcat(ones(size(X, 1)), X)
    Xpred = hcat(ones(size(Xpred, 1)), Xpred)
    return sum(abs.(X * model.B - Xpred * model.B); dims=2) / size(model.B, 2)
end

"""
    bootstrap_with_candidates(::Type{RONA}, rng::Random.AbstractRNG, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}, Xpred::AbstractMatrix{T2}, nboot::Int=500; candidates_threshold::Real=0.05, genomic_control::Bool=true, tw_threshold::Real=0.001) where {T1<:Real, T2<:Real}  
  
Compute the genomic offset for the RONA model using a bootstrap approach. For every, bootstrap iteration, the model is fitted using a random subset of the columns of `Y`. The Tracy-Widom test is used to estimate the number of latent factors. The F-test is used to select putatively adaptative loci. The genomic offset is computed for the selected loci.

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
  - λ::Real=1e-5: Regularization parameter for the LFMM.


# Returns
  - A matrix of size NxNboot with the genomic offset values. If no candidate loci are found, the row is filled with zeros.
"""
function bootstrap_with_candidates(::Type{RONA}, rng::Random.AbstractRNG,
                                   Y::AbstractMatrix{T1}, X::AbstractMatrix{T2},
                                   Xpred::AbstractMatrix{T2}, nboot::Int=500;
                                   candidates_threshold::Real=0.05,
                                   genomic_control::Bool=true,
                                   tw_threshold::Real=0.001,
                                   λ::Real=1e-5) where {T1<:Real,T2<:Real}
    _, L = size(Y)
    model = GenomicOffsets.fit(RONA, Y, X) # Fit entire model
    # Center and scale data for LFMM
    mx = mean(X; dims=1)
    sdx = std(X; dims=1)
    Xscaled = (X .- mx) ./ sdx
    Yscaled = Y .- mean(Y; dims=1)
    # Add intercept for RONA
    X = hcat(ones(size(X, 1)), X)
    Xpred = hcat(ones(size(Xpred, 1)), Xpred)
    # Allocate memory & initialize seed
    offsets = zeros(size(Y, 1), nboot)
    shared_seed = rand(rng, UInt)
    Threads.@threads for i in 1:nboot
        _rng = Random.seed!(copy(rng), shared_seed + i)
        sampled = sample(_rng, 1:L, L; replace=true)
        Yboot = Yscaled[:, sampled]
        # Fit LFMM for GEA
        lfmm = fit(GeometricGO, Yboot, Xscaled; center=false, scale=false,
                   tw_threshold=tw_threshold, λ=λ).model
        pvalues = LFMM_Ftest(lfmm, Yboot, Xscaled; genomic_control=genomic_control,
                             center=false)
        qvalues = adjust(pvalues, BenjaminiHochberg())
        candidates = findall(qvalues .< candidates_threshold)
        if length(candidates) > 0
            B = model.B[:, sampled[candidates]]
            offsets[:, i] = sum(abs.(X * B - Xpred * B); dims=2) / size(model.B, 2)
        end
    end
    return offsets
end

function bootstrap_with_candidates(::Type{RONA}, Y::AbstractMatrix{T1},
                                   X::AbstractMatrix{T2}, Xpred::AbstractMatrix{T2},
                                   nboot::Int=500; candidates_threshold::Real=0.05,
                                   genomic_control::Bool=true,
                                   tw_threshold::Real=0.001,
                                   λ::Real=1e-5) where {T1<:Real,T2<:Real}
    return bootstrap_with_candidates(RONA, Random.GLOBAL_RNG, Y, X, Xpred, nboot;
                                     candidates_threshold=candidates_threshold,
                                     genomic_control=genomic_control,
                                     tw_threshold=tw_threshold,
                                     λ=λ)
end
