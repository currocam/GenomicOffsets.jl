using Interpolations
"""
Gradient genomic offset. TODO: Add reference
"""
struct GradientForestGO <: AbstractGO
    F::Vector{Interpolations.AbstractExtrapolation{Float64,1}}
end

"""
    fit(::Type{GradientForestGO}, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}; ntrees::Int=100) where {T1<:Real,T2<:Real}

Fit the Gradient forest genomic offset model (that is, fit a Gradient Forest) using the data `Y` and `X`.

# Arguments
- `Y::AbstractMatrix{T1}`: NxL matrix with individuals genotypes or allele frequencies.
- `X::AbstractMatrix{T}`: NxP matrix with environmental data.
- `ntrees`: Number of trees each forest (one per allele) will have.
# Returns
  - A GradientForestGO model.
"""
function fit(::Type{GradientForestGO}, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2};
             ntrees::Int=100) where {T1<:Real,T2<:Real}
    nbins = Int(ceil(log2(size(Y, 1)/2)))
    gf = GenomicOffsets.GradientForest.gradient_forest(Y, X; ntrees=ntrees, nbins=nbins)
    return GradientForestGO(gf.F)
end

"""
    genomic_offset(model::GradientForestGO, X::AbstractMatrix{T}, Xpred::AbstractMatrix{T}) where T<:Real
  
Compute the Gradient Forest genomic offset. Attention! Our implementation does not perform extrapolation, which means that you may obtain unpredictable results if you try so. 

# Arguments
- `model::GradientForestGO`: A GeometricGO model.
- `X::AbstractMatrix{T}`: NxP matrix with environmental data. Data will be scaled and centered if the model was fitted with center and scale set to true.
- `Xpred::Matrix{T}`: NxP matrix with the predicted / altered environmental matrix. Data will be scaled and centered if the model was fitted with center and scale set to true.
# Returns
  - A Vector{Float64} of length N with the Gradient forest genomic offset values.
"""
function genomic_offset(model::GradientForestGO, X::AbstractMatrix{T},
                        Xpred::AbstractMatrix{T}) where {T<:Real}
    N, P = size(X)
    distance = zeros(N)
    for p in 1:P
        F = model.F[p]
        distance .+= abs2.(F.(X[:, p]) .- F.(Xpred[:, p]))
    end
    return sqrt.(distance)
end

"""
    bootstrap_with_candidates(::Type{GradientForestGO}, rng::Random.AbstractRNG,Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}, Xpred::AbstractMatrix{T2}, nboot::Int=500; ntrees::Int=100, candidates_threshold::Real=0.05, genomic_control::Bool=true, tw_threshold::Real=0.001) where {T1<:Real,T2<:Real}
  
Compute the Gradient Forest genomic offset using a bootstrap approach. For every, bootstrap iteration, the model is fitted using a random subset of the columns of `Y`. The Tracy-Widom test is used to estimate the number of latent factors. The F-test is used to select putatively adaptative loci. The genomic offset is computed for the selected loci.

# Arguments
- `rng::Random.AbstractRNG`: Random number generator. If not provided, the global RNG is used.
- `Y::AbstractMatrix{T1}`: NxL matrix with individuals genotypes or allele frequencies.
- `X::AbstractMatrix{T2}`: NxP matrix with environmental data.
- `Xpred::AbstractMatrix{T2}`: NxP matrix with the predicted / altered environmental matrix.
- `nboot::Int=500`: Number of bootstrap iterations.
- Optional arguments:
  - `ntrees::Int=100`: Number of trees each forest (one per allele) will have.
  - `candidates_threshold::Real=0.05`: Threshold for the F-test. Better to be permisive.
  - `genomic_control::Bool=true`: Apply genomic control to the F-test.
  - `tw_threshold::Real=0.001`: Threshold for the Tracy-Widom test. Better to be conservative (as the number of latent factors is often overconservative).

# Returns
  - A matrix of size NxNboot with the genomic offset values. If no candidate loci are found, the row is filled with zeros.
"""
function bootstrap_with_candidates(::Type{GradientForestGO}, rng::Random.AbstractRNG,
                                   Y::AbstractMatrix{T1}, X::AbstractMatrix{T2},
                                   Xpred::AbstractMatrix{T2}, nboot::Int=500;
                                   ntrees::Int=100, candidates_threshold::Real=0.05,
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
        Yboot = Y[:, sample(_rng, 1:L, L; replace=true)]
        eigenvalues = eigvals(Yboot * Yboot' / (size(Yboot, 1) - 1))
        _, pvalues = TracyWidom(eigenvalues)
        K = max(findfirst(pvalues .> tw_threshold) - 1, 1)
        pvalues = LFMM_Ftest(RidgeLFMM(Yboot, X, K; center=false), Yboot, X;
                             genomic_control=genomic_control, center=false)
        qvalues = adjust(pvalues, BenjaminiHochberg())
        candidates = findall(qvalues .< candidates_threshold)
        if length(candidates) > 0
            model = fit(GradientForestGO, Yboot[:, candidates], X; ntrees=ntrees)
            offsets[:, i] = genomic_offset(model, X, Xpred)
        end
    end
    return offsets
end

function bootstrap_with_candidates(::Type{GradientForestGO}, Y::AbstractMatrix{T1},
                                   X::AbstractMatrix{T2}, Xpred::AbstractMatrix{T2},
                                   nboot::Int=500; ntrees::Int=100,
                                   candidates_threshold::Real=0.05,
                                   genomic_control::Bool=true,
                                   tw_threshold::Real=0.001) where {T1<:Real,T2<:Real}
    return bootstrap_with_candidates(GradientForestGO, Random.GLOBAL_RNG, Y, X, Xpred,
                                     nboot; ntrees=ntrees,
                                     candidates_threshold=candidates_threshold,
                                     genomic_control=genomic_control,
                                     tw_threshold=tw_threshold)
end