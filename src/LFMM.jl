using LinearAlgebra, Statistics, RandomMatrices
"""
    TracyWidom(eigenvalues::AbstractVector{T}; sort_eigenvalues::Bool=true) where T<:Real

Compute the Tracy-Widom statistics and p-values for a given set of eigenvalues. TODO: add reference. 

# Arguments
- `eigenvalues::AbstractVector{T}`: A vector of eigenvalues. 
- `sort_eigenvalues::Bool=true`: If true, sort the eigenvalues in descending order. If false, the eigenvalues are assumed to be sorted in descending order.
# Returns
- A tuple with the Tracy-Widom statistics and p-values. For a given threshold, it is advised not to keep all eigenvalues below the threshold, but up to the first one that is above the threshold (not included).
"""
function TracyWidom(eigenvalues::AbstractVector{T}; sort_eigenvalues::Bool=true) where T<:Real
    if sort_eigenvalues
        eigenvalues = sort(eigenvalues, rev=true)
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