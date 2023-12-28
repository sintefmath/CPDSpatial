# Copyright (c) 2023 Odd Andersen
# Copyright (c) 2023 SINTEF Digital

using Jutul
using JutulDarcy
using SpecialFunctions # for erfinv

export JutulCPDSystem, JutulCPDModel, select_primary_variables!
export select_secondary_variables!, select_default_darcy_parameters!
export ArrheniusReactionRate, values_per_entity, setup_state
export TemperatureProfile 

include("system_definition.jl")
include("variable_definitions.jl")
include("selects.jl")
include("equations.jl")
include("convergence.jl")

