using BoltzmannSolvers
using Documenter

DocMeta.setdocmeta!(BoltzmannSolvers, :DocTestSetup, :(using BoltzmannSolvers); recursive=true)

makedocs(;
    modules=[BoltzmannSolvers],
    authors="Jan Kuhfeld <jan.kuhfeld@rub.de> and contributors",
    repo="https://github.com/jqfeld/BoltzmannSolvers.jl/blob/{commit}{path}#{line}",
    sitename="BoltzmannSolvers.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jqfeld.github.io/BoltzmannSolvers.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jqfeld/BoltzmannSolvers.jl",
    devbranch="main",
)
