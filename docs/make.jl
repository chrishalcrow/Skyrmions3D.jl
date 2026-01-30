using Documenter, Skyrmions3D, DocumenterCitations

bib = CitationBibliography(
    joinpath(@__DIR__, "../paper/paper.bib");
    style=:numeric
)

makedocs(
    warnonly=:doctest,
    modules = [Skyrmions3D],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "Skyrmion3d.jl",
    authors = "Chris Halcrow",
    pages = [
        "Home" => "index.md",
        "Make, save and load" => "making_skyrmions.md",
        "Transform and Combine" => "transform.md",
        "Compute properties" => "properties.md",
        "Flow" => "flow.md",
        "Visualise" => "visualisation.md",
        "API" => "api.md",
        "Citations" => "citations.md",
        "Contribute" => "contributing.md",
        "References" => "references.md",
    ],
    repo = Documenter.Remotes.GitHub("chrishalcrow", "Skyrmions3D.jl"),
    checkdocs=:exports;
    plugins=[bib]
)

deploydocs(; repo = "github.com/chrishalcrow/Skyrmions3D.jl")
