# GenomicOffsets.jl

[![Build Status](https://github.com/currocam/GenomicOffsets.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/currocam/GenomicOffsets.jl/actions/workflows/CI.yml?query=branch%3Amain)

Genomic offsets (GO) statistics measure current individuals/populations' maladaptation in the event of a drastic climate change. These metrics are obtained from Genotype environment associations (GEA). 

This Julia package provides efficient implementations of the four most popular metrics: [RONA]( https://doi.org/10.1111/mec.13889), [RDA GO](https://doi.org/10.1111/2041-210X.13722), [Gradient Forest GO](https://doi.org/10.1111/ele.12376) and [Geometric GO](https://doi.org/10.1093/molbev/msad140). In addition, we have implemented an F-statistical test based on the [Linear](https://doi.org/10.1093/molbev/msz008)[ Mixed Latent](https://doi.org/10.1093/molbev/msz008) Factor model](https://doi.org/10.1093/molbev/msz008) (LFMM) to perform the identification of GEA candidates through hypothesis testing and a bootstrap approach that resamples loci, performs the GEA candidates identification and computes the desired genomic offset metric. 

## Getting started (Julia users)

You can install the development package version from the REPL by entering the package mode (press `]`) and executing

```julia
pkg> add https://github.com/currocam/GenomicOffsets.jl
```
Then, you can load the package by running
```julia
julia> using GenomicOffsets
```

## Getting started (R users)

The first step is to [install Julia](https://julialang.org/downloads/). Hopefully, it should be straightforward. 
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
Now, you can install the Julia package and load it into R (so you can use the functions) as:
```R
library(JuliaConnectoR)
stopifnot(juliaSetupOk())
juliaEval('using Pkg; Pkg.add("https://github.com/currocam/GenomicOffsets.jl"))')
GO <- juliaImport("GenomicOffsets")
```
## Recipes (Julia and R)
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

#### RONA
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
Because you have to use the dollar sign operator to access the functions in the Julia environment, more complicated expressions might look too verbose. Hereafter, we will use the base R `with` function to evaluate the expression *inside* the Julia environment. You can see the difference:
```R
# R
rona <- with(GO, {
    model <- fit(RONA, Y, X)
    genomic_offset(model, X, Xpred)
})
```

#### Redundancy analysis GO (+ GEA candidates using F-test LFMM)

Now, we consider a slightly more complicated scenario where first we do hypothesis testing to identify a set of GEA candidates using an F-Test based on an LFMM model and then we compute the RDA genomic offset for that subset of alleles. 

First, we have to fit a LFMM model. For doing so, we also need to specify `K`, the number of latent factors (which can be done using a screeplot or, as we will see later, with the Tracy-Widom test). Then, we fit the RDA and compute the RDA genomic offset.

In Julia:
```julia
# Julia
K = 1 #From screeplot
lfmm = RidgeLFMM(Y, X, K)
pvalues = LFMM_Ftest(lfmm, Y, X; genomic_control=true)
# Bonferroni correction at 0.05
candidates = findall(pvalues .< 0.05 / size(Y, 2))
model = fit(RDAGO, Y[:,candidates], X)
offset = genomic_offset(model, X, Xpred)
```
Or, in R:
```R
# R
offset <- with(GO, {
    K <- 1L
    lfmm <- RidgeLFMM(Y, X, K)
    pvalues <- LFMM_Ftest(lfmm, Y, X, genomic_control=TRUE)
    # Bonferroni correction at 0.05
    candidates <-which(pvalues < 0.05 / ncol(Y))
    model <- fit(RDAGO, Y[,candidates], X)
    genomic_offset(model, X, Xpred)
})
```

When computing the RDA genomic offsets you have to decide how many canonical axes you are going to include and whether you want to weigh them based on their eigenvalues.

In Julia:
```julia
# Julia
model = fit(RDAGO, Y[:,candidates], X; K=5) #Using the 5 main canonical axes 
model = fit(RDAGO, Y[:,candidates], X; pratio=0.70) #As many axes so 70% of variances is preserved
offset = genomic_offset(model, X, Xpred)
# Weight distance with eigenvalues
offset = genomic_offset(model, X, Xpred; weighted=true)
```

In R:
```R
# R
with(GO, {
    #Using the 5 main canonical axes 
    model <- fit(RDAGO, Y[,candidates], X, K=5L) 
    #As many axes so 70% of variances is preserved
    model <- fit(RDAGO, Y[,candidates], X, pratio=0.70) 
    offset <- genomic_offset(model, X, Xpred)
    # Weight distance with eigenvalues
    genomic_offset(model, X, Xpred, weighted=TRUE)
})
```

#### Geometric GO (+ GEA candidates using F-test LFMM)

The geometric GO is based on the LFMM model. Therefore, it is very natural to use the same fitted model for both hypothesis testing and the geometric offset (and we provide a convenient interface for doing so). The F-test will find whether the environmental predictors improve the goodness of fit in comparison with just the fitted latent factors. After the set of candidates is determined, the geometric genomic offset might be computed from a subset of the effect sizes (the other part of the LFMM).

As mentioned above, when fitting a LFMM it is necessary to find the number of latent factors. We have implemented an approach based on the Tracy-Widom test to do so (although it might be over-conservative). 

In Julia:
```julia
# Julia
model = fit(GeometricGO, Y, X, 5) #Using 5 latent factors
model = fit(GeometricGO, Y, X) # By default using tw_threshold=1e-3
model = fit(GeometricGO, Y, X; tw_threshold=0.05) 
# You can check the chosen K
println(model.K)
pvalues = LFMM_Ftest(model, Y, X)
candidates = findall(pvalues .< 0.1 / size(Y, 2))
offset = genomic_offset(model, X, Xpred, candidates)
```

In R:
```R
# R
with(GO, {
    #Using 5 latent factors
    model <- fit(GeometricGO, Y, X, 5L)
    # By default using tw_threshold=1e-3
    model <- fit(GeometricGO, Y, X)
    # Now using 0.05 instead
    model <- fit(GeometricGO, Y, X, tw_threshold=0.05)
    # You can check the chosen K
    print(model$K)
    pvalues <- LFMM_Ftest(model, Y, X)
    candidates <- which(pvalues < 0.1 / ncol(Y))
    genomic_offset(model, X, Xpred, candidates)
})
```
#### Gradient Forest GO

```julia
# Julia
model = fit(GradientForestGO, Y[:,candidates], X)
model = fit(GradientForestGO, Y[:,candidates], X; ntrees=500)
offset = genomic_offset(model, X, Xpred)
```

In R:
```R
# R
with(GO, {
    #Using 5 latent factors
    model <- fit(GradientForestGO, Y[,candidates], X)
    # By default using tw_threshold=1e-3
    model <- fit(GradientForestGO,  Y[,candidates], X, ntrees=500L)
    genomic_offset(model, X, Xpred)
})
```

#### Bootstrapping genomic offsets

We provide a bootstrap implementation of all previous genomic offsets based on resampling of the different loci. At each bootstrap iteration we:

1. Take a sample with replacement of all given loci. 
2. We determine the number of latent factors, `K`,  using the Tracy-Widom test (for a given threshold). 
3. We fit an LFMM model (with `K` latent factors).
4. We compute p-values for each locus using the F-test (explained above). 
5. We adjust the p-values using the Benjamini–Hochberg procedure. 
6. We determine a set of candidates (such that the False discovery rate is below a certain threshold).
7. We compute the desired genomic offset. 

The output of the function is a matrix with as many rows as individuals/populations and as many columns as bootstrapping iterations.  

In Julia
```julia
# Julia
# Default values with RONA
offsets_boots = bootstrap_with_candidates(RONA, Y, X, Xpred)
# RDAGO with FDR < 0.1 and no genomic control
offsets_boots = bootstrap_with_candidates(RDAGO, Y, X, Xpred; candidates_threshold=0.1, genomic_control=false)
# Geometric with 200 bootstrapped samples
offsets_boots = bootstrap_with_candidates(GeometricGO, Y, X, Xpred, 200)
#  Gradient Forest with 200 trees and 100 boots
offsets_boots = bootstrap_with_candidates(GradientForestGO, Y, X, Xpred, 100; ntrees=200)
```

In R
```R
# R
with(GO, {
    # Default values with RONA
    offsets_boots <- bootstrap_with_candidates(RONA, Y, X, Xpred)
    # RDAGO with FDR < 0.1 and no genomic control
    offsets_boots <- bootstrap_with_candidates(RDAGO, Y, X, Xpred, candidates_threshold=0.1, genomic_control=FALSE)
    # Geometric with 200 bootstrapped samples
    offsets_boots <- bootstrap_with_candidates(GeometricGO, Y, X, Xpred, 20L)
    #  Gradient Forest with 200 trees and 100 boots
    offsets_boots = bootstrap_with_candidates(GradientForestGO, Y, X, Xpred, 10L, ntrees=200L)
})
```

When running this from R, you may want to take advantage of multithreading. If so, you can do it **before starting Julia** (that is, just after loading the `JuliaConnectoR` package):

```R
library(JuliaConnectoR)
Sys.setenv(JULIA_NUM_THREADS = 6)
```

#### And what more?
This README doesn't cover all details nor all optional arguments. You can read the docstrings in Julia by entering the help mode with the `?` operator. Sadly, this is not possible from R. Instead, you should take a look at the list of functions on the [documentation website](https://currocam.github.io/GenomicOffsets.jl/dev/). Please, notice that in Julia is common to overload the same functions for different *signatures*, so look for the signature you need carefully.