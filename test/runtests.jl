using GenomicOffsets
using Test
using Serialization
using LinearAlgebra
using Statistics

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
    rda = RDA_GO(data.Y, data.X, data.Xpred)
    @test norm(rda - expected_output) < 1e-10
  end

function test_tracy_widom()
  Y = open(deserialize, "data/example_data.jld").Y
  expected_output = open(deserialize, "data/tracywidom.jld")
  expected_twstat = expected_output[:statistic]
  Ycentered = Y .- mean(Y, dims=1)
  A = (Ycentered)*Ycentered'/(size(Ycentered, 1)-1)
  eigenvalues = eigvals(A)
  twstat, pvalues = GenomicOffsets.TracyWidom(eigenvalues)
  # Check that there is only one NaN in the last index
  @test isnan(twstat[end])
  @test isnan(expected_twstat[end])
  @test norm(twstat[1:end-1] - expected_twstat[1:end-1]) < 1e-5
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
end

function test_geometric()
  data = open(deserialize, "data/example_data.jld")
  expected_outputs = open(deserialize, "data/geometric.jld")
  @test norm(Geometric_GO(data.Y, data.X, data.Xpred, 2) - expected_outputs[:geometric_k2]) < 1e-10
  @test norm(Geometric_GO(data.Y, data.X, data.Xpred, 3) - expected_outputs[:geometric_k3]) < 1e-10
  @test norm(Geometric_GO(data.Y, data.X, data.Xpred, 25) - expected_outputs[:geometric_k25]) < 1e-10
  # Check that runs
  Geometric_GO(data.Y, data.X, data.Xpred)
  Geometric_GO(data.Y, data.X, data.Xpred, 2; center=false)
  Geometric_GO(data.Y, data.X, data.Xpred, 2; scale=false)
  Geometric_GO(data.Y, data.X, data.Xpred, 2; tw_threshold=0.1)
  Geometric_GO(data.Y, data.X, data.Xpred; tw_threshold=0.1)
  Geometric_GO(data.Y, data.X, data.Xpred; tw_threshold=0.005)
end

function test_gradient_forest()
  data = open(deserialize, "data/example_data.jld")
  expected_output = open(deserialize, "data/gradientforest.jld")
  # Check that, if using the same seed and same data, the model is the same
  # TODO: check seed
  #models = [GenomicOffsets.GradientForest.gradient_forest(data.Y, data.X; ntrees =10, rng=100).R2p for _ in 1:5]
  #for i in 2:5
  #  @test norm(models[1] - models[i]) < 1e-10
  #end  
  offset = GradientForest_GO(data.Y, data.X, data.Xpred; ntrees = 100)
  @test cor(offset,expected_output) > 0.7
 end

function test_Ftest()
  data = open(deserialize, "data/example_data.jld")
  expected_outputs = open(deserialize, "data/lfmm2.jld")
  @test norm(expected_outputs[:lfmm2_k2][:pvalues] .- GenomicOffsets.LFMM_Ftest(RidgeLFMM(data.Y, data.X, 2), data.Y, data.X)) < 1e-10
  @test norm(expected_outputs[:lfmm2_k3][:pvalues] .- GenomicOffsets.LFMM_Ftest(RidgeLFMM(data.Y, data.X, 3), data.Y, data.X)) < 1e-10
  @test norm(expected_outputs[:lfmm2_k25][:pvalues] .- GenomicOffsets.LFMM_Ftest(RidgeLFMM(data.Y, data.X, 25), data.Y, data.X)) < 1e-10
 end

@testset "GenomicOffsets.jl" begin
    # Write your tests here.
    test_rona()
    test_rda()
    test_tracy_widom()
    test_lfmm2()
    test_geometric()
    test_gradient_forest()
    test_Ftest()
end
