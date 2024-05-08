# GenomicOffsets.jl

[![Build Status](https://github.com/currocam/GenomicOffsets.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/currocam/GenomicOffsets.jl/actions/workflows/CI.yml?query=branch%3Amain)

Genomic offsets (GO) statistics measure current individuals / populations maladaptation in the event of a drastic climate change. These metrics are obtained from Genotype x environment associations (GEA). 

This Julia package provides efficient implementations of the four most popular metrics: [RONA]( https://doi.org/10.1111/mec.13889), [RDA GO](https://doi.org/10.1111/2041-210X.13722), [Gradient Forest GO](https://doi.org/10.1111/ele.12376) and [Geometric GO](https://doi.org/10.1093/molbev/msad140). In addition, we have implemented a F-statistical test based on [Linear Mixed Latent Factor model](https://doi.org/10.1093/molbev/msz008) (LFMM) to perform the identification of GEA candidates through hypothesis testing and a bootstrap approach that resamples loci, performs the GEA candidates identification and computes the desired genomic offset metric. When fitting a LFMM without specifying a number of latent factors, we have also implemented an automatic Tracy-Widom test for doing so (although, in our experience, it is over-conservative). 

## Getting started (Julia users)

You can install the development package version from the REPL by entering in the package mode (press `]`) and executing

```julia
pkg> add https://github.com/currocam/GenomicOffsets.jl
```
Then, you can load the package by running
```julia
julia> using GenomicOffsets
```

## Getting started (R users)

The first step is to  [install Julia](https://julialang.org/downloads/). Hopefully, it should be straighforward. 
```bash
curl -fsSL https://install.julialang.org | sh
# For windows users
# winget install julia -s msstore
# Now, Julia should be in the system search PATH
julia -v
```
The second step is to install the [JuliaConnectoR](https://github.com/stefan-m-lenz/JuliaConnectoR) package. You can do it from the R console by running:
```R
install.packages("JuliaConnectoR")
```
Now, you can install the Julia package and loading into R (so you can use the functions) as:
```R
library(JuliaConnectoR)
stopifnot(juliaSetupOk())
juliaEval('using Pkg; Pkg.add("https://github.com/currocam/GenomicOffsets.jl"))')
GO <- juliaImport("GenomicOffsets")
```
## Recipes (both Julia and R)
### Loading data
To fit the GEA the package expects `Y`, a matrix of genotypes (encoded as integers) or allele frequencies and `X`, a matrix of environmental covariates. Then, to estimate genomic offsets the package expects two matrixes of current and altered/future environmental covariates, named `Xnow` and `Xpred`. However, in most cases, `X=Xnow`. We include a small toy dataset for you to play with the package. You can import it to the environment in Julia using
```julia
# Julia
using GenomicOffsets
Y, X, Xpred = godataset
```
Or, in R;
```R
# R
GO <- juliaImport("GenomicOffsets")
Y <- GO$read_godataset()[[1]]
X <- GO$read_godataset()[[2]]
Xpred <- GO$read_godataset()[[3]]
```
You may consider using external packages such as [SnpArrays.jl](https://openmendel.github.io/SnpArrays.jl/latest/#Summaries) or [https://github.com/knausb/vcfR](vcfR) to load your data. 

### Computing any genomic offset metric

All metrics share the same interface. First, you `fit` the GEA model and then you compute the genomic offset. 

### RONA
First, in Julia:
```julia
# Julia
model = fit(RONA, Y, X)
rona = genomic_offset(model, X, Xpred)
```
Or, in R:
```R
# R
model <- GO$fit(GO$RONA, Y, X)
rona <- GO$genomic_offset(model, X, Xpred)
```
Because you have to use the dollar sign operator to access to the functions in the Julia environment, more complicated expressions might look too verbose. Hereafter, we will use the base R `with` function to evaluate the expression *inside* the Julia environment. You can see the difference:
```R
# R
rona <- with(GO, {
    model <- fit(RONA, Y, X)
    genomic_offset(model, X, Xpred)
})
```

### Redundancy analysis GO (+ GEA candidates using F-test LFMM)

Now, we consider a slightly more complicated scenario where first we do hypothesis testing to identify a set of GEA candidates using an F-Test based on an LFMM model and then we compute the RDA genomic offset for that subset of alleles. 

First, we have to fit a LFMM model. For doing so, we also need to specify `K`, the number of latent factors (which can be done using a screeplot or, as we will see later, with the Tracy-Widom test).

In Julia:

```julia
# Julia
K = 3
lfmm = RidgeLFMM(Y, X, K)
pvalues = LFMM_Ftest(lfmm, Y, X; genomic_control=true)
# Bonferroni correction at 0.05
candidates = findall(pvalues .< 0.05 / size(Y, 2))
```

Or, in R:
```R
# R
model <- GO$fit(GO$RONA, Y, X)
rona <- GO$genomic_offset(model, X, Xpred)
```
Because you have to use the dollar sign operator to access to the functions in the Julia environment, more complicated expressions might look too verbose. Hereafter, we will use the base R `with` function to evaluate the expression *inside* the Julia environment. You can see the difference:
```R
# R
rona <- with(GO, {
    model <- fit(RONA, Y, X)
    genomic_offset(model, X, Xpred)
})
```

### Geometric Genomic offset (with GEA candidates)

When computing the Geometric genomic offset, we have to fit a LFMM model. For doing so, we need to determine the number of latent factors. This can be done using a screeplot. However, we provide a handy approach that uses Tracy-Widom tests to determined it automatically. 

First, in Julia:
```julia
# Julia
model = fit(GeometricGO, Y, X) # TW for p < 1e-3
model = fit(GeometricGO, Y, X, 5) # 5 latent factors
model = fit(GeometricGO, Y, X; tw_threshold=0.05) # TW for p < 0.05
offset = genomic_offset(model, X, Xpred)
```

Or, in R (notice the `5L` that denotes 5 is an integer and the `,` rather than the `;`):
```R
# R
offset <- with(GO, {
    model <- fit(GeometricGO, Y, X) # TW for p < 1e-3
    model <- fit(GeometricGO, Y, X, 5L) # 5 latent factors
    model <- fit(GeometricGO, Y, X, tw_threshold=0.05) # TW for p < 0.05
    genomic_offset(model, X, Xpred)
})
```

We can use this same model to do hypothesis testing and identify a set of GEA candidates. Therefore, we provide a simpler interface for this user case. 
