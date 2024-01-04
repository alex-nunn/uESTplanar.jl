using uESTplanar
using Documenter

ENV["GKSwstype"] = "100"

DocMeta.setdocmeta!(uESTplanar, :DocTestSetup, :(using uESTplanar); recursive=true)

makedocs(;
    modules=[uESTplanar],
    authors="Alex Nunn <alex.nunn@pm.me> and contributors",
    repo="https://github.com/alex-nunn/uESTplanar.jl/blob/{commit}{path}#{line}",
    sitename="uESTplanar.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://alex-nunn.github.io/uESTplanar.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => [
            "Convolutions" => "examples/01_convolutions.md",
            "Potential fields" => "examples/02_potential_fields.md",
            "Electric fields" => "examples/03_electric_fields.md",
            "Ion trajectories" => "examples/04_ion_trajectories.md"
        ],
        "Reference" => "reference.md"
    ],
)

deploydocs(;
    repo="github.com/alex-nunn/uESTplanar.jl",
    devbranch="main",
)
