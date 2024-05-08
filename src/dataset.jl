using Serialization

function read_godataset()
    project_path(parts...) = normpath(joinpath(@__DIR__, "..", parts...))
    data_dict = deserialize(project_path("test/data/data.jld"))
    return (data_dict[:Y], data_dict[:X], data_dict[:Xpred])
end