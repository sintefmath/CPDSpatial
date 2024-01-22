# Copyright (c) 2023 Odd Andersen
# Copyright (c) 2023 SINTEF Digital

import Jutul: convergence_criterion
import JutulDarcy: convergence_criterion

# Define JutulCPDSystem
"""This is my docstring"""
struct JutulCPDSystem{T} <: MultiPhaseSystem where T<:Tuple
    phase::T # tuple of phases, for the moment only expected to contain a single phase
    num_tar_bins::Int
    imposed_global_temperature::Bool
    imposed_global_pressure::Float64  # set to a positive number to impose a global constant pressure
end

function JutulCPDSystem(;num_tar_bins=40,
                        imposed_global_temperature=false,
                        imposed_global_pressure=NaN)

    return JutulCPDSystem((JutulDarcy.VaporPhase(),),
                          num_tar_bins,
                          imposed_global_temperature,
                          imposed_global_pressure)
end

"""This is the docstring for JutulCPDModel"""
const JutulCPDModel = Jutul.SimulationModel{<:Any, <:JutulCPDSystem}

@inline JutulDarcy.get_phases(sys::JutulCPDSystem) = sys.phase
@inline JutulDarcy.number_of_phases(::JutulCPDSystem) = 1
@inline JutulDarcy.number_of_components(::JutulCPDSystem) = 1

"""This is the docstring for CPDFlowBoundaryCondition"""
struct CPDFlowBoundaryCondition{I, F, FUN} <: JutulForce
    cell::I
    pressure::F
    temperature::FUN
    trans_flow::F
    trans_thermal::F
end
