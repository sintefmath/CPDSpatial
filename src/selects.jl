function Jutul.setup_forces(model::JutulCPDModel, sources = nothing, bc = nothing)

    return (sources = sources, bc = bc)
end

function Jutul.select_primary_variables!(S, sys::JutulCPDSystem, model::SimulationModel)

    # Transport primary variables
    Tmin = 273.15 # We need a value > 0K to avoid potential division by zero.
                  # We assume here that we won't need temperatures below the
                  # freezing point of water.
    Tmax = 3000.0 # We assume here that we won't need temperatures above 3000K.
    Pmin = 101325.0 / 2 # One half of an atmosphere as minimum, again to avoid
                        # potential division by 0.
    Pmax = 1e8 # 100 MPa as maximum (should not be attainable in practice anyway)
    S[:Temperature] = JutulDarcy.Temperature(Tmin, Tmax, nothing, nothing)
    S[:Pressure] = JutulDarcy.Pressure(; max_rel=1.0,
                                       scale=1e5, minimum = Pmin, maximum = Pmax)

    # Note: the last bin in :ξ and :ξvapor represent light gas, the other bins represent
    S[:ξ] = DetachedMass{:Total}() # everything that is no longer part of the
                                   # solid matrix, whether vapor or liquid
    # CPD-related primary variables
    S[:£δc] = BridgeStatus() # Labile, side chains, char bridges
end

function Jutul.select_secondary_variables!(S, sys::JutulCPDSystem, model::SimulationModel)

    # Transport secondary variables
    S[:PhaseMassDensities] = VaporDensity() 
    S[:VaporInternalEnergy] = VaporInternalEnergy()
    S[:VaporEnthalpy] = VaporEnthalpy()
    S[:NonvaporInternalEnergy] = NonvaporInternalEnergy()
    S[:TotalThermalEnergy] = TotalThermalEnergy()
    
    # Reference mass values, as produced by basic CPD
    S[:g1] = LightGas{:g1}()
    S[:g2] = LightGas{:g2}()
    S[:g] = LightGas{:g}()
    S[:ftar_and_gas] = DetachedMass{:Reference}()
    S[:fchar] = SolidCharMassFraction()

    # Mass values, taking full spatial model into account
    S[:ξvapor] = DetachedMass{:Vapor}() # vaporized mass (light gas + tar)
    S[:ξmetaplast] = DetachedMass{:Metaplast}() # metaplast (liquid phase)

    S[:CurrentCellMass] = TotalCellMass() # total current mass in cell (solid, vapor and liquid)
    S[:RemainingCharVolume] = RemainingCharVolume() # volume of remaining char
    S[:MetaplastCrossLinkRate] = DetachedMass{:MetaplastCrossLinkRate}() # cross-linked metaplast rate
end

function Jutul.select_parameters!(S, sys::JutulCPDSystem, model::SimulationModel)
    
    # Transport parameters
    # Solid (i.e. non-vapor) related heat transport parameters
    S[:CharHeatCapacity] = CharHeatCapacity()
    S[:BulkDensity] = BulkDensity() # NB: this is the _initial_ density of the input material
    S[:CharDensity] = CharDensity() # The density of just the char
    S[:BulkVolume] = JutulDarcy.BulkVolume() # _initial_ bulk volume

    # Vapor heat related heat transport parameters
    S[:VaporHeatCapacity] = VaporHeatCapacity()
    S[:FluidVolume] = JutulDarcy.FluidVolume()
    S[:TotalThermalConductivity] = ThermalConductivity()
         
    # Single-phase flow parameters (common with most models in JutulDarcy as well)
    S[:PhaseViscosities] = JutulDarcy.PhaseViscosities()
    S[:PhaseMassMobilities] = JutulDarcy.PhaseMassMobilities()
    S[:RelativePermeabilities] = JutulDarcy.RelativePermeabilitiesParameter()
    S[:Transmissibilities] = JutulDarcy.Transmissibilities()

    # CPD-specific material parameters
    S[:σ] = CoordinationNumber()
    S[:c₀] = InitialBridgeFraction() # charred bridges
    S[:p₀] = InitialBridgeFraction() # intact bridges (incl. charred bridges)
    S[:mb] = MolecularMass{:mb}() # bridge mass
    S[:ma] = MolecularMass{:ma}() # site mass

    # CPD-specific reaction parameters
    S[:b_rate] = ArrheniusReactionRate{:b_rate}() # bridge breaking reaction: £ →£⋆
    S[:g_rate] = ArrheniusReactionRate{:g_rate}() # formation of light gas: δ → g₁
    S[:ρ_rate] = FixedReactionRate{:ρ_rate}()     # ratio of side chain formation to
                                         # char bridge formation: kδ/kc
    S[:cross_rate] = ArrheniusReactionRate{:cross_rate}() # cross-linking rate 

    # CPD-specific flash parameters
    S[:flash_αβγ] = CPDFlashParams() # α, β and γ used in vapor pressure formula

end

function Jutul.select_minimum_output_variables!(out, sys::JutulCPDSystem, model)

    # Transport output variables
    push!(out, :TotalThermalEnergy)
    push!(out, :Temperature)
    push!(out, :Pressure)

    # CPD-specific output variables
    push!(out, :ξ)
    push!(out, :ξvapor)
    push!(out, :g1)
    push!(out, :g2)
    push!(out, :g)
    push!(out, :ftar_and_gas)
    push!(out, :fchar)

    # Values necessary for post-processing
    push!(out, :FluidVolume) # with evolving porosity, this value should also change
    #push!(out, :Porosity)
    push!(out, :MetaplastCrossLinkRate)

    push!(out, :PhaseMassDensities)
    push!(out, :NonvaporInternalEnergy)
    push!(out, :VaporInternalEnergy)
    push!(out, :PhaseViscosities) # might change with temperature, and necessary for
                                  # reconstructing boundary fluxes
end

function Jutul.select_equations!(eqs, sys::JutulCPDSystem, model::SimulationModel)

    # Temperature equation
    if sys.imposed_global_temperature
        # global temperature is imposed directly through force term
        eqs[:temperature_equation] = TrivialEquation{:Temperature, 0.0}()
    else
        # temperature determined by energy conservation
        hdisc = model.domain.discretizations.heat_flow
        eqs[:energy_conservation] = ConservationLaw(hdisc, :TotalThermalEnergy)
    end

    # # Pressure equation
    eqs[:pressure_equation] = CPDPressureEquation()
    
    # Conservation of vapor masses
    mdisc = model.domain.discretizations.mass_flow
    eqs[:mass_conservation] = ConservationLaw(mdisc, :ξ, model.system.num_tar_bins + 1)

    # CPD-specific equations
    eqs[:£δc] = £δcODEEquation()
end
