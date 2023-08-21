
using SmoQyDEAC
using Documenter
using Literate

example_names = ["fermion_greens"]
example_literate_sources = [joinpath(@__DIR__,"src/examples/"*name*".jl") for name in example_names]
example_script_destinations = [joinpath(@__DIR__,"../scripts") for name in example_names]
example_documentation_destination = joinpath(@__DIR__,"src/examples")
example_documentation_paths = [joinpath(@__DIR__,"src/examples/$name.md") for name in example_names]



DocMeta.setdocmeta!(SmoQyDEAC, :DocTestSetup, :(using SmoQyDEAC); recursive=true)

for i in eachindex(example_names)
    Literate.markdown(example_literate_sources[i], example_documentation_destination; 
                      execute = false,
                      documenter = true)
    Literate.script(example_literate_sources[i], example_script_destinations[i])
end

makedocs(;
    modules=[SmoQyDEAC],
    authors="James Neuhaus <jneuhau1@utk.edu>",
    repo="https://github.com/sandimas/SmoQyDEAC.jl",
    sitename="SmoQyDEAC.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://sandimas.github.io/SmoQyDEAC.jl/",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => "examples/fermion_greens.md",
        "Derivations" => "derivations.md",
    ],
    draft = false
)

deploydocs(;
    repo="github.com/sandimas/SmoQyDEAC.jl",
    devbranch="main",
)

