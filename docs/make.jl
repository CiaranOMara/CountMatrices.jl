using Pkg
using CountMatrices
using Documenter

DocMeta.setdocmeta!(CountMatrices, :DocTestSetup, :(using CountMatrices); recursive=true)

makedocs(;
    modules=[CountMatrices],
    authors=join(Pkg.TOML.parsefile("Project.toml")["authors"]),
    repo="https://github.com/CiaranOMara/CountMatrices.jl/blob/{commit}{path}#{line}",
    sitename="CountMatrices.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://CiaranOMara.github.io/CountMatrices.jl",
        edit_link="develop",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/CiaranOMara/CountMatrices.jl",
    devbranch="develop",
)
