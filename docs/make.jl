using Documenter, CPDSpatial

push!(LOAD_PATH, "../src/")

makedocs(;
         #modules=[CPDSpatial],
         sitename="CPDSpatial.jl",
         remotes=nothing,
         format=Documenter.HTML(; prettyurls=(ENV, "CI", "false") == "true"),
         pages = [
             "Home" => "index.md",
             "Chemical Percolation and Devolatilization" => "cpd.md",
             "Spatially embedded CPD" => "spatial.md",
             "Examples" => "examples.md",
             "Docstrings" => "docstrings.md",
         ]
         )

