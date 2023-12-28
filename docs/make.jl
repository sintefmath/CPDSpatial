using Documenter, CPDSpatial

push!(LOAD_PATH, "../src/")
makedocs(sitename="CPDSpatial.jl",
         remotes=nothing,
         format=Documenter.HTML(; prettyurls=(ENV, "CI", "false") == "true"),
         pages = [
             "Home" => "home.md",
             "CPD" => "cpd.md",
             "Spatial" => "spatial.md"]
         )

