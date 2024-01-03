using uESTplanar
using Documenter

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
    ],
)

deploydocs(;
    repo="github.com/alex-nunn/uESTplanar.jl",
    devbranch="main",
)
