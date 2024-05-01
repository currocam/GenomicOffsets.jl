using GenomicOffsets
using Test
using Serialization
using LinearAlgebra

function test_rona()
  data = open(deserialize, "data/example_data.jld")
  expected_output = open(deserialize, "data/rona.jld")
  rona = RONA(data.Y, data.X, data.Xpred)
  @test norm(rona - expected_output) < 1e-10
end

@testset "GenomicOffsets.jl" begin
    # Write your tests here.
    test_rona()
end
