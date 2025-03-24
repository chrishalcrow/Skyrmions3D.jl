using Documenter, Skyrmions3D

makedocs(
    modules = [Skyrmions3D],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "Skyrmion3d.jl",
    authors  = "Chris Halcrow",
    pages = [
        "Home" => "index.md",
        "Create Skyrmion" => "making_skyrmions.md",
        "Transform and Combine" => "transform.md",
        "Compute properties" => "properties.md",
        "Flow" => "flow.md",
        "API" => "api.md"
    ]
)
