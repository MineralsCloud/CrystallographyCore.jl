using CrystallographyCore
using Documenter

DocMeta.setdocmeta!(CrystallographyCore, :DocTestSetup, :(using CrystallographyCore); recursive=true)

makedocs(;
    modules=[CrystallographyCore],
    authors="singularitti <singularitti@outlook.com> and contributors",
    repo="https://github.com/MineralsCloud/CrystallographyCore.jl/blob/{commit}{path}#{line}",
    sitename="CrystallographyCore.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/CrystallographyCore.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/CrystallographyCore.jl",
    devbranch="main",
)
