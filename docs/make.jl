push!(LOAD_PATH,"../src/")
using Documenter
using GenomicOffsets

makedocs(
    sitename = "GenomicOffsets",
    format = Documenter.HTML(),
    modules = [GenomicOffsets]
)