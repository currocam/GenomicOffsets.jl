# GenomicOffsets.jl Documentation

The GenomicOffsets Julia package implements efficient versions of some of the most popular Genomic offset metrics.

Disclaimer: this documentation is not finished and just provides a quick overview of the syntax. 

## Getting started

The package provides a miminal dataset you can use. You can load it as follos. 

```@example
using GenomicOffsets
Y, X, Xpred = data; 
```

On the one hand, `Y` is the genotype matrix (or allele frequencies) of the current locally adapted individuals (or populations). We have 100 individuals and 500 loci. 

```@example
using GenomicOffsets # hide
Y, X, Xpred = data; # hide
size(Y)
```

On the other hand, `X` is the current environmental matrix. We have 100 individuals and 2 ecological predictors.  `Xpred` are the same individuals (although this is not mandatory) for which we have a forecasted ecological predictors.  

```@example
using GenomicOffsets # hide
Y, X, Xpred = data; # hide
@assert size(X) == size(Xpred)
size(X)
```

We can compute the RONA (and the rest of supported methods) by first, fitting the model and then computing the predicted genomic offset. 

```@example
using GenomicOffsets # hide
Y, X, Xpred = data # hide
rona = fit(RONA, Y, X)
genomic_offset(rona, X, Xpred)
```

Other genomic offset metrics have additional parameters. For example, if computing the geometric genomic offset, you may want to provide the number of latent factors of the LFMM model (and not use the Tracy-Widom test to determined it). 

```@example
using GenomicOffsets # hide
Y, X, Xpred = data # hide
geometric = fit(GeometricGO, Y, X, 3)
genomic_offset(geometric, X, Xpred)
```

To see a more complicated use case, let's compute the 95% confidence niterval using  bootstrapping. We do it by computing the geometric genomic offset repeatly with different loci sampled with replacement. At every iteration, (a) the number of latent factors will be determined using the TracyWidom test, (b) a set of putatively adaptative loci will be identified with an F-test (corrected using genomic control) based on the LFMM model and (c) the geometric genomic offste is computed.  

```@example
using GenomicOffsets # hide
Y, X, Xpred = data # hide
using Random, StatsBase
nboots = 100
rng = Xoshiro(123) # random seed
offsets = bootstrap_with_candidates(GeometricGO, rng, 
    Y, X, Xpred, nboots;candidates_threshold=0.05)
@assert size(offsets) == (size(Y, 1), nboots) 
confint95 = [quantile(ind, [0.025, 0.975]) for ind in eachrow(offsets)]
```

## References

```@autodocs
Modules = [GenomicOffsets]
Order   = [:function, :type]
```