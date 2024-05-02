using GenomicOffsets
using Test
using Serialization
using LinearAlgebra
using Statistics

function test_rona()
  data = open(deserialize, "data/example_data.jld")
  expected_output = open(deserialize, "data/rona.jld")
  rona = RONA(data.Y, data.X, data.Xpred)
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

@testset "GenomicOffsets.jl" begin
    # Write your tests here.
    test_rona()
    test_rda()
    test_tracy_widom()
end
