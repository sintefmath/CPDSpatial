__precompile__(true)

module CPDSpatial

export JutulCPDSystem, select_primary_variables!, select_secondary_variables!
export compute_yield_curves, simple_test_domain, radial_test_domain, setup_state, compute_reattached_metaplast
export compute_internal_mass_fluxes
import Jutul
import JutulDarcy

include("ad_adaptations.jl")
include("percolation.jl")
include("cpd.jl")
include("flash.jl")
include("defaults.jl")
include("JutulCPDModel.jl")
include("tools.jl")

end # module
