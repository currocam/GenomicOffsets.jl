module GenomicOffsets
abstract type AbstractGO end

export RONA, fit, genomic_offset, RDAGO, TracyWidom, RidgeLFMM, GeometricGO,
       GradientForestGO, LFMM_Ftest, bootstrap_with_candidates, data

include("dataset.jl")
include("RONA.jl")
include("RDA_GO.jl")
include("LFMM.jl")
include("Geometric_GO.jl")
include("GradientForest.jl")
include("GradientForest_GO.jl")
include("GEACandidates.jl")
const data = dataset()
end
