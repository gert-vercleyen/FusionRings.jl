
using Documenter, FusionRings

makedocs(
    sitename = "FusionRings.jl",
    modules = [FusionRings],
    format = Documenter.HTML(prettyurls=false),
    pages = [
        "Home" => "index.md",
    ],
)
