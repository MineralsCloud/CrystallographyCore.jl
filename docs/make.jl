using CrystallographyCore
using Documenter

DocMeta.setdocmeta!(CrystallographyCore, :DocTestSetup, :(using CrystallographyCore, Unitful, UnitfulAtomic); recursive=true)

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
        "Manual" => [
            "Installation Guide" => "man/installation.md",
            "Definitions and conventions" => "man/definitions.md",
            "Migration Guide (`v0.3.x` to `v0.7.x`)" => "man/migration-v0.3-to-v0.7.md",
            "Troubleshooting" => "man/troubleshooting.md",
        ],
        "Reference" => Any[
            "Public API" => "lib/public.md",
            "Internals" => map(
                s -> "lib/internals/$(s)",
                sort(readdir(joinpath(@__DIR__, "src/lib/internals")))
            ),
        ],
        "Developer Docs" => [
            "Contributing" => "developers/contributing.md",
            "Style Guide" => "developers/style-guide.md",
            "Design Principles" => "developers/design-principles.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/CrystallographyCore.jl",
    devbranch="main",
)
