using GenomicOffsets
using Test
using Serialization
using LinearAlgebra
using Statistics
using Random

function test_rona()
    data = open(deserialize, "data/example_data.jld")
    expected_output = open(deserialize, "data/rona.jld")
    model = fit(RONA, data.Y, data.X)
    rona = genomic_offset(model, data.X, data.Xpred)
    @test norm(rona - expected_output) < 1e-10
end

function test_rda()
    data = open(deserialize, "data/example_data.jld")
    expected_output = open(deserialize, "data/rda.jld")
    model = fit(RDAGO, data.Y, data.X)
    rda = genomic_offset(model, data.X, data.Xpred)
    @test norm(rda - expected_output) < 1e-10
    _ = genomic_offset(model, data.X, data.Xpred; weighted=true)
    return nothing
end

function test_tracy_widom()
    Y = open(deserialize, "data/example_data.jld").Y
    expected_output = open(deserialize, "data/tracywidom.jld")
    expected_twstat = expected_output[:statistic]
    Ycentered = Y .- mean(Y; dims=1)
    A = (Ycentered) * Ycentered' / (size(Ycentered, 1) - 1)
    eigenvalues = eigvals(A)
    twstat, pvalues = GenomicOffsets.TracyWidom(eigenvalues)
    # Check that there is only one NaN in the last index
    @test isnan(twstat[end])
    @test isnan(expected_twstat[end])
    @test norm(twstat[1:(end - 1)] - expected_twstat[1:(end - 1)]) < 1e-5
    # R implementation doesn't use p-values.  Wwe rely on RandomMatrices.jl to generate the Tracy-Widom distribution
    # rather than having a lookup table. Then, we just check the number of significant eigenvalues are the same.
    @test findfirst(pvalues .> 0.01) - 1 == Int(expected_output[:SigntEigenL])
end

function test_lfmm2()
    data = open(deserialize, "data/example_data.jld")
    expected_outputs = open(deserialize, "data/lfmm2.jld")
    Y = data.Y
    X = data.X
    function check_lfmm(expected, actual)
        @test norm(expected[:U] - actual.U) < 1e-10
        @test norm(expected[:V] - actual.Vt') < 1e-10
        @test norm(expected[:B] - actual.Bt') < 1e-10
    end
    check_lfmm(expected_outputs[:lfmm2_k2], RidgeLFMM(Y, X, 2))
    check_lfmm(expected_outputs[:lfmm2_k3], RidgeLFMM(Y, X, 3))
    check_lfmm(expected_outputs[:lfmm2_k25], RidgeLFMM(Y, X, 25))
    return nothing
end

function test_geometric()
    data = open(deserialize, "data/example_data.jld")
    expected_outputs = open(deserialize, "data/geometric.jld")
    models = [fit(GeometricGO, data.Y, data.X, k) for k in [2, 3, 25]]
    expected_offsets = [expected_outputs[:geometric_k2], expected_outputs[:geometric_k3],
                        expected_outputs[:geometric_k25]]
    for (model, offset) in zip(models, expected_offsets)
        @test norm(genomic_offset(model, data.X, data.Xpred) - offset) < 1e-10
    end
    # Check that runs
    genomic_offset(fit(GeometricGO, data.Y, data.X), data.X, data.Xpred)
    genomic_offset(fit(GeometricGO, data.Y, data.X, 2; center=false), data.X, data.Xpred)
    genomic_offset(fit(GeometricGO, data.Y, data.X, 2; scale=false), data.X, data.Xpred)
    genomic_offset(fit(GeometricGO, data.Y, data.X, 2; tw_threshold=0.1), data.X,
                   data.Xpred)
    genomic_offset(fit(GeometricGO, data.Y, data.X; tw_threshold=0.1), data.X, data.Xpred)
    genomic_offset(fit(GeometricGO, data.Y, data.X; tw_threshold=0.005), data.X,
                   data.Xpred)
    return nothing
end

function test_gradient_forest()
    data = open(deserialize, "data/example_data.jld")
    expected_output = open(deserialize, "data/gradientforest.jld")
    # TODO: random seeds
    model = fit(GradientForestGO, data.Y, data.X; ntrees=100)
    offset = genomic_offset(model, data.X, data.Xpred)
    @test cor(offset, expected_output) > 0.7
end

function test_Ftest()
    data = open(deserialize, "data/example_data.jld")
    expected_outputs = open(deserialize, "data/lfmm2.jld")
    @test norm(expected_outputs[:lfmm2_k2][:pvalues] .-
               GenomicOffsets.LFMM_Ftest(RidgeLFMM(data.Y, data.X, 2), data.Y, data.X)) <
          1e-10
    @test norm(expected_outputs[:lfmm2_k3][:pvalues] .-
               GenomicOffsets.LFMM_Ftest(RidgeLFMM(data.Y, data.X, 3), data.Y, data.X)) <
          1e-10
    @test norm(expected_outputs[:lfmm2_k25][:pvalues] .-
               GenomicOffsets.LFMM_Ftest(RidgeLFMM(data.Y, data.X, 25), data.Y, data.X)) <
          1e-10
    @test norm(expected_outputs[:lfmm2_k2][:pvalues] .-
               GenomicOffsets.LFMM_Ftest(fit(GeometricGO, data.Y, data.X, 2), data.Y,
                                         data.X)) <
          1e-4
    @test norm(expected_outputs[:lfmm2_k3][:pvalues] .-
               GenomicOffsets.LFMM_Ftest(fit(GeometricGO, data.Y, data.X, 3), data.Y,
                                         data.X)) <
          1e-4
    @test norm(expected_outputs[:lfmm2_k25][:pvalues] .-
               GenomicOffsets.LFMM_Ftest(fit(GeometricGO, data.Y, data.X, 25), data.Y,
                                         data.X)) <
          1e-4
end

function bootstraps()
    data = open(deserialize, "data/data.jld")
    neglogfitness = data[:neglogfitness]
    Y = data[:Y]
    X = data[:X]
    Xpred = data[:Xpred]
    for type in [RONA, RDAGO, GeometricGO]
        rng = Random.default_rng()
        boots = bootstrap_with_candidates(type, copy(rng), Y, X, Xpred,
                                          100)
        @test boots ≈
              bootstrap_with_candidates(type, copy(rng), Y, X, Xpred,
                                        100)
        @test median([cor(neglogfitness, col) for col in eachcol(boots)]) > 0.70
    end
    # TODO: check GradientForest
    boots = bootstrap_with_candidates(GradientForestGO, Y, X, Xpred, 100;
                                      ntrees=100)
    @test median([cor(neglogfitness, col) for col in eachcol(boots)]) > 0.40
    @test size(boots) == (size(Xpred, 1), 100)
end

@testset "GenomicOffsets.jl" begin
    test_rona()
    test_rda()
    test_tracy_widom()
    test_lfmm2()
    test_geometric()
    test_gradient_forest()
    test_Ftest()
    bootstraps()
end
