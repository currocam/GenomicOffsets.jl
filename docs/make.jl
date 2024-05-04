push!(LOAD_PATH,"../src/")
using Documenter
using GenomicOffsets

makedocs(
    sitename = "GenomicOffsets",
    format = Documenter.HTML(),
    modules = [GenomicOffsets]
)

deploydocs(
    repo = "github.com/currocam/GenomicOffsets.jl.git",
)