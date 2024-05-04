using Serialization

function dataset()
    project_path(parts...) = normpath(joinpath(@__DIR__, "..", parts...))
    deserialize(project_path("test/data/data.jld"))
end