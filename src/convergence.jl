# Copyright (c) 2023 Odd Andersen
# Copyright (c) 2023 SINTEF Digital

# Override convergence criterion for mass conservation law in JutulDarcy, since
# we need to account for the possibility that a substance may have zero density
function Jutul.convergence_criterion(model::JutulCPDModel,
                                     storage, eq::ConservationLaw{:TotalCellMass}, eq_s, r;
                                     dt = 1, kwarg...)

    # compare error with initial cell mass
    cell_masses = storage.state.BulkDensity .* storage.state.BulkVolume

    rel_error = dt .* [maximum(abs.(r[:] ./ cell_masses))]

    names = ["Tot. Mass"]
    R = (Rel = (errors = rel_error, names = names), )
    return R

    
end

# Override convergence criterion for mass conservation law in JutulDarcy, since
# we need to account for the possibility that a substance may have zero density
function Jutul.convergence_criterion(model::JutulCPDModel,
                                     storage, eq::ConservationLaw{:ξ}, eq_s, r;
                                     dt = 1, kwarg...)
    M = global_map(model.domain)
    v = x -> as_value(Jutul.active_view(x, M, for_variables = false))
    Φ = v(storage.state.FluidVolume)
    ρ = v(storage.state.PhaseMassDensities)

    ρ = sum(ρ, dims=1) # sum over all components in the phase

    # handle the case where a substance has zero density
    ρ = Matrix(ρ)
    ρ[ρ .== 0] .= 1e-40
    
    nph = number_of_phases(model.system)
    cnv, mb = JutulDarcy.cnv_mb_errors(r, Φ, ρ, dt, Val(nph))

    names = JutulDarcy.phase_names(model.system)
    R = (CNV = (errors = cnv, names = names),
         MB = (errors = mb, names = names))
    return R
end

# Override convergence criterion for energy conservation law in JutulDarcy, since most of the
# energy is found in the solid, not the vapor phase
function Jutul.convergence_criterion(model::JutulCPDModel, storage,
                                     eq::ConservationLaw{:TotalThermalEnergy}, eq_s, r;
                                     dt = 1, kwarg...)

    # ensure no division by zero, and prevent convergence issues if there are
    # many orders of difference between the thermal energy in different cells
    # (primarily caused by different cell sizes)
    minimum_fac = 1e-2 # if the energy content of a cell is two orders of
                       # magnitude less than the mean we relax convergence
                       # criterion
    mean_energy = sum(storage.state.TotalThermalEnergy) / length(r)
    small = mean_energy == 0.0 ? eps() : minimum_fac * mean_energy

    # compute "CNV-like" error
    cnv_error = dt * [maximum(abs.(r[:] ./ max.(storage.state.TotalThermalEnergy, small)))]

    # compute "mb"-like error
    mb_error = dt * [abs(sum(r[:])) / sum(storage.state.TotalThermalEnergy)]

    names = ["Energy"]
    R = (CNV = (errors = cnv_error, names = names),
         MB = (errors = mb_error, names = names))
    return R
end

function Jutul.convergence_criterion(model::JutulCPDModel,
                                     storage, eq::CPDPressureEquation, eq_s, r;
                                     dt = 1, kwarg...)

    mean_pressure = sum(storage.state.Pressure) / length(r)
    min_fac = 1e-2 # if the pressure in a cell is two orders of magnitude less
                   # than the mean we relax relative convergence criterion

    # As above for energy conservation law
    small = mean_pressure == 0.0 ? eps() : min_fac * mean_pressure 

    rel_error = [maximum(abs.(r[:] ./ max.(storage.state.Pressure, small)))]
    

    names = ["Pressure"]
    R = (RelMax = (errors = rel_error, names = names), )
    return R
end
