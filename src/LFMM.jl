using LinearAlgebra, Statistics, RandomMatrices, Distributions, MultipleTesting

"""
    TracyWidom(eigenvalues::AbstractVector{T}; sort_eigenvalues::Bool=true) where T<:Real

Compute the Tracy-Widom statistics and p-values for a given set of eigenvalues. TODO: add reference. 

# Arguments
- `eigenvalues::AbstractVector{T}`: A vector of eigenvalues. 
- `sort_eigenvalues::Bool=true`: If true, sort the eigenvalues in descending order. If false, the eigenvalues are assumed to be sorted in descending order.
# Returns
- A tuple with the Tracy-Widom statistics and p-values. For a given threshold, it is advised not to keep all eigenvalues below the threshold, but up to the first one that is above the threshold (not included).
"""
function TracyWidom(eigenvalues::AbstractVector{T};
                    sort_eigenvalues::Bool=true) where {T<:Real}
    if sort_eigenvalues
        eigenvalues = sort(eigenvalues; rev=true)
    end
    L1 = reverse(cumsum(reverse(eigenvalues)))
    L2 = reverse(cumsum(reverse(eigenvalues .^ 2)))
    N = length(eigenvalues):-1:1
    S2 = N .^ 2 .* L2 ./ (L1 .^ 2)
    v = N .* (N .+ 2) ./ (S2 .- N)
    L = N .* eigenvalues ./ L1
    v_st = ifelse.(v .- 1 .> 0, v .- 1, NaN)
    v_st = .√v_st
    N_st = .√(N)
    μ = (v_st .+ N_st) .^ 2 ./ v
    σ = (v_st .+ N_st) ./ v .* (1 ./ v_st .+ 1 ./ N_st) .^ (1 / 3)
    twstat = (L .- μ) ./ σ
    pvalues = 1 .- cdf(RandomMatrices.TracyWidom(1), twstat)
    return twstat, pvalues
end

"""
    LFMM{T<:Real}

A struct to store the results of the Linear Fixed Effects Mixed Model (LFMM) as described by Caye et al. (2019).

```math
Y = X B^T + W = X B^T + U V^T
```
"""
struct LFMM{T<:Real}
    U::AbstractMatrix{T}
    Vt::AbstractMatrix{T}
    Bt::AbstractMatrix{T}
end

"""
    RidgeLFMM(Y::Matrix{T}, X::Matrix{T}, K::Int, λ::Float64) where T<:Real

Ridge solutions for the Linear Fixed Effects Mixed Model (LFMM) as described by Caye et al. (2019).

```math
Y = X B^T + W = X B^T + U V^T
```

# Arguments
- `Y`: A centered genotype matrix of size NxL.
- `X`: Environmental matrix of size NxP.
- `K`: Number of latent factors.
- `λ`: Regularization parameter (by 1e-5)
- `center`: If true, both the genotype matrix and the environmental matrix are centered. If false, both matrices are assumed to be centered.
# Returns
- A `LFMM{T<:Real}` data structure with the latent factors `U`, `Vt`, and the effect sizes `Bt`.

"""
function RidgeLFMM(Y::Matrix{T1}, X::Matrix{T2}, K::Int, λ=1e-5;
                   center=true) where {T1<:Real,T2<:Real}
    if size(Y, 1) != size(X, 1)
        throw(DimensionMismatch("The number of rows in Y and X must be equal."))
    end
    if center
        Y = Y .- mean(Y; dims=1)
        X = (X .- mean(X; dims=1))
    end
    n, p = size(X)
    Q, Σ, _ = svd(X; full=true)
    d = ones(n)
    @. d[1:p] = √(λ / (λ + Σ))
    D = Diagonal(d)
    # SVD of modified Y
    u, s, v = svd(D * Q'Y)
    U = Q * D^-1 * view(u, :, 1:K) * Diagonal(view(s, 1:K))
    Vt = view(v, :, 1:K)'
    Bt = (X'X + λ * I)^-1 * X' * (Y - U * Vt)
    return LFMM(U, Vt, Bt)
end