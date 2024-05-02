module GenomicOffsets
abstract type AbstractGO end
export RONA
export fit
export genomic_offset
export RDAGO
export TracyWidom
export RidgeLFMM
export GeometricGO
export GradientForestGO
export bootstrap

include("RONA.jl")
include("RDA_GO.jl")
include("LFMM.jl")
include("Geometric_GO.jl")
include("GradientForest.jl")
include("GradientForest_GO.jl")
include("GEACandidates.jl")
include("bootstraps.jl")


end
