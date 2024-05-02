module GradientForest
using Random, KernelDensity, StatsBase, DecisionTree, Interpolations, NumericalIntegration

# Density of predictors
##  Whitened Gaussian kernel with bandwidth given by Silverman’s rule-of-thumb p. 160
function predictors_density(X, npoints = 2^10, λ=0.9)
    densities = [kde(pred; npoints=npoints) for pred in eachcol(X)]
    for index in eachindex(densities)
        densities[index].density = λ.*densities[index].density .+ 1 .- λ
    end
    return densities
end

# Compute the importance of each feature as the increase in the OOB mean square prediction error
# https://scikit-learn.org/stable/modules/permutation_importance.html
function permutation_importance!(importances, y, X, tree, K=30, rng=Random.GLOBAL_RNG::Random.AbstractRNG)
    _, P = size(X)
    baseline = mean((y .- apply_tree(tree, X)).^2)
    for p in 1:P
        original = copy(X[:, p])
        for _ in 1:K
            X[:, p] = shuffle(rng, X[:, p])
            importances[p] += mean((y .- apply_tree(tree, X)).^2)
        end
        X[:, p] = original
        importances[p] = importances[p] / K - baseline
    end
end

# Compute the impurity of a node (for continuous valuqes)
function impurity(values)
    n = length(values)
    if n == 0
        return 0
    end
    return var(values)
end

# Assign a value to a bin
function bin(x, grid, low = 1, high = length(grid))
    while low < high
        mid = low + (high - low) ÷ 2
        if grid[mid] < x
            low = mid + 1
        else
            high = mid
        end
    end
    return high
end

# Traverse the tree and update the raw importances
function update_raw_importances!(Ipb, densities, node, inds, y, X)
    queue = [(node, inds)]
    while !isempty(queue)
        (node, inds) = popfirst!(queue)
        if node.featid < 1
            continue
        end
        condition = X[inds, node.featid] .< node.featval
        left_indexes = inds[condition]
        right_indexes = inds[.!condition]
        split_imp = impurity(view(y, inds)) - impurity(view(y, left_indexes)) - impurity(view(y, right_indexes))
        Ipb[node.featid, bin(node.featval, densities[node.featid].x)] += split_imp
        if typeof(node.left) != Leaf{Float64} && sum(condition) > 0
            push!(queue, (node.left, inds[condition]))
        end
        if typeof(node.right) != Leaf && length(inds) - sum(condition) > 0
            push!(queue, (node.right, inds[.!condition]))
        end
    end

end

# Compute R2 from OOB predictions
function compute_r2(y, yhat)
    return 1 - sum((y .- yhat).^2) / (var(y) * length(yhat))
end
# Compute the compositional turnover function
function composite_turnover(x, f)
    F = cumul_integrate(x, f)
    extrapolate(interpolate(x, F, FiniteDifferenceMonotonicInterpolation()), Flat())
end

function gradient_forest(Y, X; ntrees =500, nbins = 2^10, mtry = ceil(size(X, 2) / 3), npermiter = 1)
    N, P = size(X)
    _, L = size(Y)
    # goodness-of-fit measure for the forest for allele f
    # Computed using the out-of-bag samples
    R2f = zeros(L)
    # Accuracy importance for predictor p within the forest for allele f
    # Computed using the out-of-bag samples as permutation_importance
    Ipf = zeros(P, L)
    # Partitionated goodness-of-fit measure
    # Computed using R2f and Ipf
    R2pf = zeros(P, L)
    # Sum of the in-bag importances for predictor p within the forest for allele f
    # for each split the values is in bin(s) normalized by the density
    # Ifpb = \sum_{t, s \in BIN(x)} I_{fpts} / d(s)
    Ifpb = zeros(L, P, nbins)
    
    # Compute the density of each predictor
    densities = predictors_density(X, nbins)
    rng = Random.GLOBAL_RNG
    shared_seed = rand(rng, UInt)
    for f in 1:L
        y = view(Y,:, f)
        # Allocate yhat with NaN
        yhat = Vector{Float64}(undef, N)
        yhat_counter = zeros(Int64, N) # How many times a sample is out-of-bag
        permimp = zeros(P) # Permutation importance for each predictor
        forest = Vector{DecisionTree.Root}(undef, ntrees)
        inbags = Vector{Vector{Int}}(undef, ntrees)
        Threads.@threads for i in 1:ntrees
            # The Mersenne Twister (Julia's default) is not thread-safe.
            _rng = Random.seed!(copy(rng), shared_seed + i)
            inbag = sample(_rng, 1:N, N, replace=true) # In-bag samples
            inbags[i] = inbag
            # Train the tree
            forest[i] = build_tree(
                y[inbag],
                X[inbag, :],
                mtry, # n_subfeatures
                -1, # max_depth
                5, # min_samples_leaf
                2, # min_samples_split
                0.0; # min_purity_increase
                impurity_importance=false,
                rng=_rng
            )
        end
        for i in 1:ntrees
            inbag = inbags[i]
            outofbag = setdiff(1:N, inbag)
            tree = forest[i]
            # Out-of-bag predictions for R2f
            # we add, because we are going to divide by the counter later
            yhat[outofbag] += apply_tree(tree, view(X, outofbag, :))
            yhat_counter[outofbag] .+= 1
            # Compute the importance of each feature using permutation importance
            permutation_importance!(permimp, y[outofbag], X[outofbag, :], tree, npermiter)
            # Now the tricky part, we need to accumulate the in-bag importances
            # for each split (that uses the predictor p in the bin s)
            update_raw_importances!(view(Ifpb, f, :, :), densities, tree.node, inbag, y, X)
        end
        # Fix the out-of-bag predictions
        yhat[yhat_counter .> 0] ./= yhat_counter[yhat_counter .> 0]
        # Compute the R2f
        R2f[f] = compute_r2(y[yhat_counter .> 0], yhat[yhat_counter .> 0])
        # Registed the increase in the OOB mean square prediction error
        Ipf[:, f] = max.(permimp, 0)
        # Normalize the in-bag raw importances by the density
        for p in 1:P
            for s in 1:nbins
                Ifpb[f, p, s] /= densities[p].density[s]
            end
        end 
        # Compute the partitionated goodness-of-fit measure
        R2pf[:, f] = R2f[f] .* Ipf[:, f] / sum(Ipf[:, f])
    end
    
    # Overall importance of predictor p 
    R2p = sum(R2pf, dims=2) / L
    
    # Compute the sum across the bins
    sumIpf = reshape(sum(Ifpb, dims=3), P, L)
    # Ip(x)\Delta X / dp(x)
    Ip_dp = zeros(P, nbins)
    for p in 1:P
        for f in 1:L
            for s in 1:nbins
                Ip_dp[p, s] += R2pf[p, f] * Ifpb[f, p, s] / sumIpf[p, f]
            end
        end
        Ip_dp[p, :] ./= L * step(densities[p].x)
    end
    # Now we have to compute the composite cumulative importance
    F = [composite_turnover(densities[p].x, Ip_dp[p, :]) for p in 1:P]
    (R2p=R2p, F = F)  
end
end