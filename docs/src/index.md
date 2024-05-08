# GenomicOffsets.jl Documentation

Genomic offsets (GO) statistics measure current individuals/populations' maladaptation in the event of a drastic climate change. These metrics are obtained from Genotype environment associations (GEA). 

This Julia package provides efficient implementations of the four most popular metrics: [RONA]( https://doi.org/10.1111/mec.13889), [RDA GO](https://doi.org/10.1111/2041-210X.13722), [Gradient Forest GO](https://doi.org/10.1111/ele.12376) and [Geometric GO](https://doi.org/10.1093/molbev/msad140). In addition, we have implemented an F-statistical test based on the [Linear](https://doi.org/10.1093/molbev/msz008)[ Mixed Latent](https://doi.org/10.1093/molbev/msz008) Factor model](https://doi.org/10.1093/molbev/msz008) (LFMM) to perform the identification of GEA candidates through hypothesis testing and a bootstrap approach that resamples loci, performs the GEA candidates identification and computes the desired genomic offset metric. 

We have documented several "recipes" in both Julia and R in the [README file of the repository](https://github.com/currocam/GenomicOffsets.jl) (so please take a look there and use this documentation as a reference for all available functions).

## TL;TR

You can compute different genomic offsets by fitting different models and computing the genomic offset.

```@example
# Julia
using GenomicOffsets
Y, X, Xpred = godataset
rona = fit(RONA, Y, X)
genomic_offset(rona, X, Xpred)
```

You can fit LFMM models to do hypothesis testing and identify GEA candidates. 

```@example
# Julia
using GenomicOffsets # hide
Y, X, Xpred = godataset # hide
lfmm = RidgeLFMM(Y, X, 1)
pvalues = LFMM_Ftest(lfmm, Y, X; genomic_control=true)
# Bonferroni correction at 0.05
candidates = findall(pvalues .< 0.05 / size(Y, 2))
```

Yo can compute bootstrapping confidence intervals for different genomic offset metrics using an approach we propose that sample *loci* and identifies GEA candidates at each iteration: 

```@example
using GenomicOffsets # hide
Y, X, Xpred = godataset # hide
using Random, StatsBase
nboots = 100
offsets = bootstrap_with_candidates(GeometricGO, Y, X, Xpred, nboots;candidates_threshold=0.05)
@assert size(offsets) == (size(Y, 1), nboots) 
confint95 = [quantile(ind, [0.025, 0.975]) for ind in eachrow(offsets)]
```

## References

```@autodocs
Modules = [GenomicOffsets]
Order   = [:type,:function]
```