using Documenter, CPDSpatial

push!(LOAD_PATH, "../src/")

makedocs(;
         #modules=[CPDSpatial],
         sitename="CPDSpatial.jl",
         remotes=nothing,
         format=Documenter.HTML(; prettyurls=(ENV, "CI", "false") == "true"),
         pages = [
             "Home" => "index.md",
             "CPD" => "cpd.md",
             "Docstrings" => "docstrings.md",
             "Spatial" => "spatial.md"]
         )

