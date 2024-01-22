using Documenter, CPDSpatial, Literate

push!(LOAD_PATH, "../src/")

# Literate.markdown("../examples/spatial_cpd_example.jl",
#                   "src/";
#                   execute=false)

# Literate.markdown("../examples/biochar_example.jl",
#                   "src/";
#                   execute=false)

makedocs(;
         sitename="CPDSpatial.jl",
         remotes=nothing,
         format=Documenter.HTML(repolink="https://github.com/Oddan/CPDSpatial";
                                prettyurls=(ENV, "CI", "false") == "true"),
         pages = [
             "Home" => "index.md",
             "Chemical Percolation and Devolatilization" => "cpd.md",
             "Spatially embedded CPD" => "spatial.md",
             "Examples" => ["CPD without spatial model" => "cpd_example.md",
                            "CPD in spatial model" => "spatial_cpd_example.md",
                            "Composite biomaterial case" => "biochar_example.md"],
             "Docstrings" => "docstrings.md",
         ]
         )

