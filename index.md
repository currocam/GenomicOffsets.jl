# GenomicOffsets.jl Documentation

The GenomicOffsets Julia package implements efficient versions of some of the most popular Genomic offset metrics. The package expects a genotype matrix, as well as the corresponding environmental matrix. Once you have fitted the model, you can compute the offset between some current population by providing a current environmental matrix (most frequently, the same environmental matrix you used to fit the model) and the modified one. 

```@contents
```

## Risk of non-adaptness

Fit the model:

```@docs
    fit(::Type{RONA}, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}) where {T1<:Real,T2<:Real}
```

Compute the genomic offset:

```@docs
    genomic_offset(model::RONA, X::AbstractMatrix{T}, Xpred::AbstractMatrix{T}) where T<:Real
```

Do bootstrapping:

```@docs
    bootstrap_with_candidates(::Type{RONA}, rng::Random.AbstractRNG, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}, Xpred::AbstractMatrix{T2}, nboot::Int=500; candidates_threshold::Real=0.05, genomic_control::Bool=true, tw_threshold::Real=0.001) where {T1<:Real, T2<:Real}
```


## Redundancy analysis genomic offset

Fit the model:

```@docs
    fit(::Type{RDAGO}, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}) where {T1<:Real,T2<:Real}
```

Compute the genomic offset:

```@docs
    genomic_offset(model::RDAGO, X::AbstractMatrix{T}, Xpred::AbstractMatrix{T}) where T<:Real
```

Do bootstrapping:

```@docs
    bootstrap_with_candidates(::Type{RDAGO}, rng::Random.AbstractRNG, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}, Xpred::AbstractMatrix{T2}, nboot::Int=500; candidates_threshold::Real=0.05, genomic_control::Bool=true, tw_threshold::Real=0.001) where {T1<:Real, T2<:Real}
```

## Gradient Forest genomic offset

Fit the model:

```@docs
    fit(::Type{GradientForestGO}, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}; ntrees::Int=100) where {T1<:Real,T2<:Real}
```

Compute the genomic offset:

```@docs
    genomic_offset(model::GradientForestGO, X::AbstractMatrix{T}, Xpred::AbstractMatrix{T}) where T<:Real
```

Do bootstrapping:

```@docs
    bootstrap_with_candidates(::Type{GradientForestGO}, rng::Random.AbstractRNG,Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}, Xpred::AbstractMatrix{T2}, nboot::Int=500; ntrees::Int=100, candidates_threshold::Real=0.05, genomic_control::Bool=true, tw_threshold::Real=0.001) where {T1<:Real,T2<:Real}
```

## Index

```@index
```