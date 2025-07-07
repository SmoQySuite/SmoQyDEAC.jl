
using SmoQyDEAC
using Documenter
using Literate
using DocumenterCitations


example_names = ["fermion_greens","SmoQyDQMC","user_mutation"]
example_literate_sources = [joinpath(@__DIR__,"src/examples/$name.jl") for name in example_names]
example_script_destinations = [joinpath(@__DIR__,"../scripts") for name in example_names]
example_documentation_destination = joinpath(@__DIR__,"src/examples")
example_documentation_paths = [("examples/$name.md") for name in example_names]

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "references.bib");
    style=:numeric
)

DocMeta.setdocmeta!(SmoQyDEAC, :DocTestSetup, :(using SmoQyDEAC); recursive=true)

for i in eachindex(example_names)
    Literate.markdown(example_literate_sources[i], example_documentation_destination; 
                      execute = false,
                      documenter = true)
    Literate.script(example_literate_sources[i], example_script_destinations[i])
end

makedocs(;
    modules=[SmoQyDEAC],
    plugins=[bib],
    authors="James Neuhaus <jneuhau1@utk.edu>",
    repo="https://github.com/SmoQySuite/SmoQyDEAC.jl",
    sitename="SmoQyDEAC.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://SmoQySuite.github.io/SmoQyDEAC.jl/stable/",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => example_documentation_paths,
        "Derivations" => "derivations.md",
    ],
    draft = false
)

deploydocs(;
    repo="github.com/SmoQySuite/SmoQyDEAC.jl",
    devbranch="main",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#" ]
)

