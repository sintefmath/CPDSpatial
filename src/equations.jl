# Copyright (c) 2023 Odd Andersen
# Copyright (c) 2023 SINTEF Digital

export compute_boundary_mass_flux


struct CPDPressureEquation <: Jutul.JutulEquation end

function Jutul.number_of_equations_per_entity(model::SimulationModel, eq::CPDPressureEquation)
    return 1
end

struct £δcODEEquation <: Jutul.JutulEquation end

function Jutul.number_of_equations_per_entity(model::SimulationModel, eq::£δcODEEquation)
    return 3
end

function JutulDarcy.local_discretization(eq::CPDSpatial.CPDPressureEquation, i) return nothing end
function JutulDarcy.local_discretization(eq::CPDSpatial.£δcODEEquation, i) return nothing end

function Jutul.update_equation_in_entity!(target, cell, state, state0,
                                          eq::£δcODEEquation, model, dt, ldisc=nothing)

    # initial fraction of labile bridges equal initial fraction of intact
    # bridges that are not charred bridges
    £₀ = state.p₀[cell] - state.c₀[cell] 

    # starts at 0 and increases as £ decreases.  When £ = £₀, £frac = 1
    £frac = 1.0 - state.£δc[1, cell] / £₀ 

    # should always be the case, but nonlinear solver may throw different values
    # at it
    £frac = max(min(1.0, £frac), 0.0) 

    # fraction of light gas formed compared to total possible
    gfrac = state.g[cell] / (2 * (1 - state.c₀[cell])) 

    kb = compute_reaction_rate(state.b_rate[:, cell], state.Temperature[cell], £frac)
    kg = compute_reaction_rate(state.g_rate[:, cell], state.Temperature[cell], gfrac)
    kρ = compute_reaction_rate(state.ρ_rate[cell], 0.0, state.Temperature[cell])

    M = [-kb             0.0  0.0;
         2*kρ*kb/(kρ+1)  -kg  0.0;
         kb/(kρ+1)       0.0  0.0]         
    
    ∂£δc = state.£δc[:,cell] - state0.£δc[:, cell]

    target[:] = ∂£δc .- dt .*  (M * state.£δc[:, cell])
    
end

function Jutul.update_equation_in_entity!(target, cell, state, state0, eq::CPDPressureEquation,
                                          model, dt, ldisc=nothing)
    N = model.system.num_tar_bins

    light_gas_weight = state.mb[cell] / 2.0  # @@ would this really be the case??
    r = state.mb[cell] / state.ma[cell]
    £ = state.£δc[1, cell]
    δ = state.£δc[2, cell]
    c = state.£δc[3, cell]
    p = £ + c # all intact bridges

    weights = binned_molecular_weights(state.ma[cell], r, £ + c, state.σ[cell], δ, £, N)
        
    tot_moles = sum(state.ξvapor[:, cell] ./ [weights..., light_gas_weight])

    # Use ideal gas law to calculate pres[sure
    R = 8.31446261815324 # J/(mol K)
    pressure = tot_moles * R * state.Temperature[cell] / state.FluidVolume[cell]

    target[] = state[:Pressure][cell] - pressure
end



function compute_activation_energy(E₀, σ, frac)
    # compute the effective activation energy after the reaction process has run a
    # certain fraction `frac` of its course

    SMALL = 1e-3 # If we let `frac` reach 0.0 or 1.0, we get so strong nonlinear
                 # behavior that the nonlinear solver will have trouble
                 # converging.  Moreover, this part of the range is only
                 # relevant for extremely short intervals of time.
    res = E₀ + sqrt(2) * σ * erfinv(clamp(2 * frac - 1, -1 + SMALL, 1 - SMALL))
    return res
end
function compute_reaction_rate(A, E, T)
    # NB: 'R' here is the universal gas constant, which is used in computing
    # reaction rates.  Note that the universal gas constant uses a different
    # unit here (cal/mol/K) than in the transport part of the system, where the
    # SI unit is used (J/mol/K) in the computation of vapor density
    R = 1.9872036 # universal gas constant in cal/mol/K
    return A * exp(-E / (R * T))
end

function compute_reaction_rate(AEσ::Vector{V}, T) where V <: Real
    # If this version of the function was called, the activation energy is assumed to be
    # constant, and the progression of the reaction is not taken into account
    @assert AEσ[3] == 0.0 "This function should only be called for reactions with constant activation energy"
    return compute_reaction_rate(AEσ[1], AEσ[2], T)
end

function compute_reaction_rate(AEσ::Vector{V}, T, frac) where V <: Real
    E = compute_activation_energy(AEσ[2], AEσ[3], frac)
    return compute_reaction_rate(AEσ[1], E, T)
end

struct TrivialEquation{K, V} <: Jutul.JutulEquation where {K <: Symbol, V} end 

function JutulDarcy.local_discretization(eq::CPDSpatial.TrivialEquation{K, V}, i) where {K, V}
    return nothing
end
function JutulDarcy.update_equation_in_entity!(
    eq_buf, self_cell, state, state0,
    eq::CPDSpatial.TrivialEquation{K, V}, model, dt, ldisc=0) where {K, V}
    eq_buf[] = state[K][self_cell] - V
end

function Jutul.apply_forces_to_equation!(diag_part, storage, model,
                                         eq::TrivialEquation{:Temperature, V},
                                         eq_s, force::TemperatureProfile, time) where V
    # apply same temperature to all cells
    T = force.Tfun(time)
    @inbounds for cell in 1:length(diag_part)
        diag_part[cell] -= T
    end
end

function face_volume_flux(left, right, face, face_sign, state, model)
    # Specific version for tpfa flux
    kgrad = TPFA(left, right, face_sign)
    upw = SPU(left, right)
    ft = Jutul.flux_type(model.equations[:ξ_conservation])
    
    # Override the standard flux computation defined in JutulDarcy, since our
    # setting is different (e.g. we do not define densities)
    # flux_primitives removed in JutulDarcy from commit 58543bc9
    #∇p, T_f, gΔz  = JutulDarcy.flux_primitives(face, state, model, ft, kgrad, upw) 
    ∇p = JutulDarcy.pressure_gradient(state, kgrad)
    T_f = JutulDarcy.effective_transmissibility(state, face, kgrad)
    
    # volumetric vapor flux
    qtmp = - (T_f .* ∇p)
    λ = 1.0 / JutulDarcy.phase_upwind(upw, state.PhaseViscosities, 1, qtmp)
    q = λ * qtmp

    #q = q * 0; # @@@@@ Turn off all gas flux
    return q, upw
end

# Face flux for gas/vapor
@inline function Jutul.face_flux!(q_i, left, right, face, face_sign,
                                  eq::ConservationLaw{:ξ, <:Any},
                                  state, model::JutulCPDModel, dt,
                                  flow_disc::TwoPointPotentialFlowHardCoded)

    q, upw = face_volume_flux(left, right, face, face_sign, state, model)
    
    upwind_rho = JutulDarcy.upwind(upw, cell -> state.PhaseMassDensities[:, cell], q)
    
    ξflux = q .* upwind_rho
    
    return ξflux
end

# Face flux for total mass
@inline function Jutul.face_flux!(q_i, left, right, face, face_sign,
                                  eq::ConservationLaw{:TotalCellMass, <:Any},
                                  state, model::JutulCPDModel, dt,
                                  flow_disc::TwoPointPotentialFlowHardCoded)

    vapor_flux = Jutul.face_flux!(q_i, left, right, face, face_sign,
                                  model.equations[:ξ_conservation],
                                  state, model, dt, flow_disc)
    return JutulDarcy.setindex(q_i, sum(vapor_flux), 1)
end

# Face flux for thermal energy
@inline function Jutul.face_flux!(Q, left, right, face, face_sign,
                                  eq::ConservationLaw{:TotalThermalEnergy, <:Any},
                                  state, model::JutulCPDModel, dt,
                                  flow_disc::TwoPointPotentialFlowHardCoded)
    # Specific version for tpfa flux
    grad = TPFA(left, right, face_sign)
    upw = SPU(left, right)
    qv, upw_v = face_volume_flux(left, right, face, face_sign, state, model)
    q = thermal_heat_flux!(face, state, model, grad, upw, qv, upw_v)
    
    return JutulDarcy.setindex(Q, q, 1)
end

function thermal_heat_flux!(face, state, model, grad, upw, qv, upw_v)
    T = state.Temperature
    H = state.VaporEnthalpy

    λ = state.TotalThermalConductivity[face]
    
    conductive_flux = -λ * JutulDarcy.gradient(T, grad)

    #flow_common = JutulDarcy.kgrad_common(face, state, model, grad)
    H_face = JutulDarcy.phase_upwind(upw_v, H, 1, qv)

    convective_flux = H_face * qv
    return conductive_flux + convective_flux
end

function Jutul.apply_forces_to_equation!(acc, storage, model,
                                         eq::ConservationLaw{:TotalThermalEnergy}, 
                                         eq_s, force::V, time) where {V <: AbstractVector{<:CPDFlowBoundaryCondition}}

    hflux = compute_boundary_heat_flux(force, storage.state, time)
    
    for bc_ix in 1:length(force)
        acc_i = view(acc, :, force[bc_ix].cell)
        acc_i[] += hflux[bc_ix]
    end
end

function Jutul.apply_forces_to_equation!(acc, storage, model::JutulCPDModel,
                                         eq::ConservationLaw{:ξ}, eq_s,
                                         force::V, time) where {V <: AbstractVector{<:CPDFlowBoundaryCondition}}

    fluxes = compute_boundary_mass_flux(force, storage.state)

    for bc_ix in 1:length(force)
        acc_i = view(acc, :, force[bc_ix].cell)
        acc_i[:] += fluxes[bc_ix]
    end
end

function Jutul.apply_forces_to_equation!(acc, storage, model::JutulCPDModel,
                                         eq::ConservationLaw{:TotalCellMass}, eq_s,
                                         force::V, time) where {V <: AbstractVector{<:CPDFlowBoundaryCondition}}

    fluxes = compute_boundary_mass_flux(force, storage.state)

    for bc_ix in 1:length(force)
        acc_i = view(acc, :, force[bc_ix].cell)
        acc_i[] += sum(fluxes[bc_ix])
    end
end

function compute_boundary_heat_flux(force_bc, state, time)
    heat_fluxes = Vector{typeof(state[:Temperature][1])}(undef, length(force_bc))

    T = state[:Temperature]
    P = state[:Pressure]
    
    for bc_ix = 1:length(force_bc)

        # conductive flux
        bc = force_bc[bc_ix]
        T_f = bc.trans_thermal # thermal flow
        ΔT = T[bc.cell] - bc.temperature(time) 
        
        heat_fluxes[bc_ix] = T_f * ΔT

        # add in convective flux
        T_mf = bc.trans_flow # mass flow  
        Δp = P[bc.cell] - bc.pressure
        q = T_mf * Δp
        if q > 0
            mu = state[:PhaseViscosities][bc.cell]
            mflux = (q / mu) * state[:VaporEnthalpy][bc.cell]
            heat_fluxes[bc_ix] += mflux 
        end
    end
    return heat_fluxes
end

function compute_boundary_mass_flux(force_bc, state)
    fluxes = [Vector{typeof(state[:ξvapor][1, 1])}() for i = 1:length(force_bc)]

    p = state[:Pressure]
    
    for bc_ix in 1:length(force_bc)
        bc = force_bc[bc_ix]
        c = bc.cell
        T_f = bc.trans_flow
        Δp = p[c] - bc.pressure
        q = T_f * Δp

        mu = state[:PhaseViscosities][c]
        
        if q >= 0
            rhos = state[:PhaseMassDensities][:, c]
            fluxes[bc_ix] = (q / mu) * rhos
        else
            # influx of gas from surroundings into particle
            # set zero flux, but preserve possible AD structure
            fluxes[bc_ix] = 0.0 * (q / mu) * state[:ξvapor][:, c] / state[:FluidVolume][c]
        end
    end
    
    return fluxes
end

function Jutul.update_equation!(eq_s::ConservationLawTPFAStorage,
                                law::ConservationLaw{:TotalThermalEnergy},
                                storage,
                                model::JutulCPDModel, dt)
    # Zero out any sparse indices
    Jutul.reset_sources!(eq_s)

    # Next, update accumulation, "intrinsic" sources and fluxes
    Jutul.@tic "accumulation" update_accumulation!(eq_s, law, storage, model, dt)
    Jutul.@tic "fluxes" update_half_face_flux!(eq_s, law, storage, model, dt) 

end

function Jutul.update_equation!(eq_s::ConservationLawTPFAStorage,
                                law::ConservationLaw{:TotalCellMass},
                                storage,
                                model::JutulCPDModel, dt)
    # Zero out any sparse indices
    Jutul.reset_sources!(eq_s)

    # Next, update accumulation and fluxes
    Jutul.@tic "accumulation" update_accumulation!(eq_s, law, storage, model, dt)
    Jutul.@tic "fluxes" update_half_face_flux!(eq_s, law, storage, model, dt)  
end

function Jutul.update_equation!(eq_s::ConservationLawTPFAStorage,
                                law::ConservationLaw{:ξ},
                                storage, model::JutulCPDModel, dt) 
    # Zero out any sparse indices
    Jutul.reset_sources!(eq_s)

    # Next, update accumulation, "intrinsic" sources and fluxes
    Jutul.@tic "accumulation" update_accumulation!(eq_s, law, storage, model, dt)

    Jutul.@tic "fluxes" update_half_face_flux!(eq_s, law, storage, model, dt)  

    # add production of new detached fragments from CPD model
    Jutul.@tic "reaction" update_CPD_reaction!(eq_s, storage, model, dt) 
    
end

function update_CPD_reaction!(eq_s::ConservationLawTPFAStorage, storage, model, dt)
    # Update the accumulation term by adding the impact of reaction

    acc = get_entries(eq_s.accumulation) # target

    for cell in 1:model.domain.entities[Cells()]
        src = cpd_increment(cell, storage.state, storage.state0, model, dt)
        if any(isnan.(src))
            error("NaN in CPD reaction")
        end
        
        acc[:, cell] -= src
    end
    
end

function cpd_increment(cell, state, state0, model, dt)

    # initial cell mass  @@ could be precomputed
    init_cell_mass = state.BulkDensity[cell] * state.BulkVolume[cell] 

    # Compute reference values for production of new tar
    prev_ftar = state0.ftar_and_gas[1:end-1, cell]
    cur_ftar = state.ftar_and_gas[1:end-1, cell]

    # Compute reference values for production of new light gas
    prev_lightgas = state0.ftar_and_gas[end, cell]
    cur_lightgas = state.ftar_and_gas[end, cell]

    char_depletion = sum(cur_ftar) - sum(prev_ftar) + sum(cur_lightgas) - sum(prev_lightgas)
    char_depletion = max(char_depletion, 0.0) # Prevent roundoff errors
                                              # from causing negative values
    char_depletion *= init_cell_mass

    # compute new light gas produced.  The reference values are based on the
    # _original_ amount of mass, so we need to rescale the production based on the
    # _remaining_ mass
    new_lightgas = cur_lightgas - prev_lightgas
    new_lightgas = new_lightgas .* state.TotalCellMass[cell]

    # Compute amount of metaplast that is reattached into the matrix (@@ and becoming inert)
    reattached = min.(state.MetaplastCrossLinkRate[1:end-1, cell] * dt,
                      state.ξmetaplast[1:end-1, cell])
    @assert all(reattached .>= 0.0)
    
    # compute change in volatile tar masses, and append the light gas change
    ξ_increment = init_cell_mass * (cur_ftar .- prev_ftar)

    # Rescale increment in tar masses to balance with the actual char depletion
    # production of light gas and reattached metaplast
    ξ_increment = balance_tar_masses(ξ_increment, char_depletion, new_lightgas,
                                     state.ξmetaplast[:, cell], reattached)

    # Deducting reattached metaplast from increment 
    ξ_increment -= reattached  

    return vcat(ξ_increment, new_lightgas) / dt 
end

function balance_tar_masses(ξ_increment, char_depletion, new_lightgas,
                            ξ_metaplast, rettached)
    # We aim to balance ξ_increment (as computed from the basic CPD model) with the
    # actual char depletion, production of light gas and reattached metaplast, to ensure
    # mass balance.

    # Compute the depletable metaplast, i.e. the part of current metaplast that
    # will not be reattached
    ξ_depletable = ξ_metaplast[1:end-1] - rettached
    @assert all(ξ_depletable .+ eps() .>= 0.0)
    total_depletable = sum(ξ_depletable)
    @assert total_depletable + eps() >= 0.0

    tar_mass_delta = char_depletion - new_lightgas # production/depletion of tar should
                                                   # balance depletion of char and
                                                   # production of light gas
    if tar_mass_delta + total_depletable + eps() < 0
        @assert tar_mass_delta < 0
        # not enough metaplast left for the prescribed depletion.  
        tar_mass_delta = -total_depletable
    end

    # We have now determined the total tar increment/decrement value: tar_mass_delta.
    # We now need to rescale each bin in ξ_increment
    ξ_inc_sum = sum(ξ_increment)
    if abs(ξ_inc_sum) > eps() && (tar_mass_delta * ξ_inc_sum > 0)
        # Both quantities have same sign; apply a simple rescaling
        ξ_increment *= tar_mass_delta / ξ_inc_sum
    else
        # heuristic: add uniform vector
        alpha = (tar_mass_delta - ξ_inc_sum) / length(ξ_increment)
        ξ_increment .+= alpha * ones(length(ξ_increment))
    end
    #@assert abs(sum(ξ_increment) - tar_mass_delta) <= 1e-10 * max(abs(tar_mass_delta), abs(ξ_inc_sum))

    # # add to metaplast, ensure nonzero components, and rescale if necessary
    ξ_depletable_incremented = ξ_depletable + ξ_increment
    tot_remaining = sum(ξ_depletable_incremented)

    @assert tot_remaining + eps() >= 0.0

    ξ_depletable_incremented[ξ_depletable_incremented .< 0.0] .= 0.0 # ensure nonnegative components
    tmp = sum(ξ_depletable_incremented)
    if tmp > 0.0
        ξ_depletable_incremented = ξ_depletable_incremented * tot_remaining / tmp
    end

    # # the increment should now be of the right total value (tar_mass_delta),
    # # and ensure nonnegativity when combined with the remaining metaplast
    ξ_increment = ξ_depletable_incremented - ξ_depletable

    return ξ_increment
end
