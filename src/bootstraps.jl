using Random
function bootstrap(::Type{GO}, rng::Random.AbstractRNG, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}, Xpred::AbstractMatrix{T2}; nboot::Int=500, candidatesfn = LFMMCandidates) where {GO<:AbstractGO, T1<:Real, T2<:Real}
    Y = Y .- mean(Y, dims=1)
    _, L = size(Y)
    offsets = zeros(size(Y, 1), nboot)
    shared_seed = rand(rng, UInt)
    Threads.@threads for i in 1:nboot
      _rng = Random.seed!(copy(rng), shared_seed + i)
      idx = sample(_rng, 1:L, L, replace=true)
      Yboot = Y[:,idx]
      candidates = candidatesfn(Yboot, X)
      model = fit(GO, Yboot[:,candidates], X)
      offsets[:,i] = genomic_offset(model, X, Xpred)
    end
    return offsets
end

function bootstrap(::Type{GradientForestGO}, rng::Random.AbstractRNG, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}, Xpred::AbstractMatrix{T2}; ntrees = 100, nboot::Int=500, candidatesfn = LFMMCandidates) where {T1<:Real, T2<:Real}
    Y = Y .- mean(Y, dims=1)
    _, L = size(Y)
    offsets = zeros(size(Y, 1), nboot)
    shared_seed = rand(rng, UInt)
    Threads.@threads for i in 1:nboot
      _rng = Random.seed!(copy(rng), shared_seed + i)
      idx = sample(_rng, 1:L, L, replace=true)
      Yboot = Y[:,idx]
      candidates = candidatesfn(Yboot, X)
      model = fit(GradientForestGO, Yboot[:,candidates], X; ntrees=ntrees)
      offsets[:,i] = genomic_offset(model, X, Xpred)
    end
    return offsets 
end

bootstrap(::Type{GO}, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}, Xpred::AbstractMatrix{T2}; nboot::Int=500, candidatesfn = LFMMCandidates) where {GO<:AbstractGO, T1<:Real,T2<:Real} = bootstrap(GO, Random.GLOBAL_RNG, Y, X, Xpred; nboot=nboot, candidatesfn=candidatesfn)