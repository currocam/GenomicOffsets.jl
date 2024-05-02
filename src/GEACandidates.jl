"""
    LFMM_Ftest(model::LFMM, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}; genomic_control::Bool=true, center=true) where {T1<:Real,T2<:Real}

Compute F-Statistic (F-test) for the Linear Fixed Effects Mixed Model (LFMM). TODO: Add reference & description.

# Arguments
- `model`: A `LFMM{T<:Real}` data structure (obtained from `RidgeLFMM`).
- `Y`: A centered genotype matrix of size NxL.
- `X`: A centered environmental matrix of size NxP.
- `genomic_control`: If true, apply genomic control to the F-statistic.
- `center`: If true, both the genotype matrix and the environmental matrix are centered. If false, both matrices are assumed to be centered.
# Returns
- A vector of p-values for loci in the genotype matrix.
"""
function LFMM_Ftest(model::LFMM, Y::AbstractMatrix{T1}, X::AbstractMatrix{T2}; genomic_control::Bool=true, center=true) where {T1<:Real,T2<:Real}
    if center
        Y = Y .- mean(Y, dims=1)
        X = (X .- mean(X, dims=1))
    end
    n = size(X, 1)
    d = size(X, 2)
    U = [ones(n) model.U]
    # Add the intercept
    ## partial regression
    res_Y = Y - U * (U \ Y)
    res_X = X - U * (U \ X)
    res_X = [ones(n) res_X]
    fitted = res_X * (res_X \ res_Y)
    # F-test
    mss = sum(abs2.(fitted), dims=1)
    resvar = sum(abs2.(res_Y - fitted), dims=1) / (n - d - 1)
    resvar = max.(resvar, 1e-10) # numerical problems
    Fscores = mss ./ d ./ resvar
    dist = FDist(d, n - d - 1)
    if genomic_control
        gif = median(Fscores) / quantile(dist, 0.5)
        Fscores = Fscores ./ gif
    end
    pvalues = 1 .- cdf.(dist, Fscores)
    return vec(pvalues)
end