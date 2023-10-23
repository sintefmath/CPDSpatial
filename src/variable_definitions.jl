abstract type BridgeFraction <: Jutul.ScalarVariable end
abstract type BridgeFractionTriple <: Jutul.VectorVariables end
abstract type MassFraction <: Jutul.ScalarVariable end

Jutul.variable_scale(::BridgeFraction) = 1.0
Jutul.minimum_value(::BridgeFraction) = 0.0
Jutul.maximum_value(::BridgeFraction) = 2.0 # note that fraction of side chains to bridges may be > 1.
Jutul.associated_entity(::BridgeFraction) = Jutul.Cells()

# Define primary CPD variable, other than pressure and temperature
struct BridgeStatus <: BridgeFractionTriple end # 3-tuple: LabileBridges, SideChains, CharBridges
Jutul.minimum_value(::BridgeStatus) = 0.0
Jutul.maximum_value(::BridgeStatus) = 2.0 # note that fraction of side chains to bridges may be > 1.
struct DetachedMass{M} <: Jutul.VectorVariables where M <: Symbol end # binned mass of detached components 
Jutul.minimum_value(::DetachedMass) = 0.0
Jutul.maximum_value(::DetachedMass) = Inf # note that fraction of side chains to bridges may be > 1.
Jutul.variable_scale(::DetachedMass) = 1.0; 

#Jutul.absolute_increment(ξ::DetachedMass) = p.max_abs
#Jutul.relative_increment_limit(ξ::DetachedMass) = 1e1

struct TotalCellMass <: Jutul.ScalarVariable end
Jutul.minimum_value(::TotalCellMass) = 0.0

Jutul.associated_entity(::DetachedMass) = Jutul.Cells()

# Define secondary variables
struct LightGas{M} <: BridgeFraction where M <: Symbol end
struct LightGasMassFraction <: MassFraction end
#struct OrganicVolatilesMassFraction <: MassFraction end
struct SolidCharMassFraction <: MassFraction end

# Define parameters
struct CoordinationNumber <: Jutul.ScalarVariable end
struct FixedReactionRate{M} <: Jutul.ScalarVariable where M <: Symbol end
struct ArrheniusReactionRate{M} <: Jutul.VectorVariables where M <: Symbol end 
struct InitialBridgeFraction <: BridgeFraction end
struct MolecularMass{M} <: Jutul.ScalarVariable where M <: Symbol end
struct CPDFlashParams <:Jutul.VectorVariables end

function Jutul.values_per_entity(model, ::ArrheniusReactionRate{M}) where M
    return 3
end
function Jutul.values_per_entity(model, bs::BridgeStatus) 
    return Jutul.degrees_of_freedom_per_entity(model, bs)
end

function Jutul.degrees_of_freedom_per_entity(model, ::BridgeStatus) 
    return 3 # 3-tuple: LabileBridges, SideChains, CharBridges
end

function Jutul.values_per_entity(model, dm::DetachedMass) 
    # If N is the number of tar bins, then bin n < N represents the amount of
    # n-mers, whereas bin N represents the combined amount of N+k-mers where
    # k = [0 ... ∞)
    # Bin N+1 represents light gas.
    return model.system.num_tar_bins + 1 
end

function Jutul.degrees_of_freedom_per_entity(model, dm::DetachedMass)
    return Jutul.values_per_entity(model, dm)
end

Jutul.values_per_entity(model, ::CPDFlashParams) = 3 # α, β, γ

Jutul.associated_entity(::Union{CoordinationNumber,
                                FixedReactionRate,
                                ArrheniusReactionRate,
                                BridgeFraction,
                                BridgeFractionTriple,
                                MolecularMass}) = Cells()

@jutul_secondary function update_variable!(target, v::LightGas{:g1}, model, £δc, ix)
    @inbounds for i in ix
        £, δ, c = £δc[:, i]
        g1 = 2 * (1 - £ - c) - δ  # light gas from side chain conversion: 2 * (1-p) - δ
        target[i] = g1
    end
end

@jutul_secondary function update_variable!(target, v::LightGas{:g2}, model, £δc, c₀, ix)
    @inbounds for i in ix
        c = £δc[3, i]
        target[i] = 2 * (c - c₀[i])      # light gas from char bridge formation: 2 * (c - c0)
    end
end


@jutul_secondary function update_variable!(target, v::LightGas{:g}, model, g1, g2, ix)
    @inbounds for i in ix
        target[i] = g1[i] + g2[i] # return total light gas
    end
end

@jutul_secondary function update_variable!(target, v::DetachedMass{:MetaplastCrossLinkRate},
                                           model, cross_rate, Temperature, ξmetaplast, ix)
    @inbounds for i in ix
        c_rate = compute_reaction_rate(cross_rate[:, i], Temperature[i], 0.0)
        target[:, i] = c_rate .* ξmetaplast[:, i]
    end
end

@jutul_secondary function update_variable!(target, v::DetachedMass{:Reference}, model, mb, ma, σ, c₀, g, £δc, ix)
    # NB: DetachedMass{:Reference} represents a reference value, as computed by the basic CPD model
    @inbounds for i in ix
        r = mb[i] / ma[i]
        £ = min(£δc[1, i], 1.0) # nonlinear solver may cause this one to surpass 1
        δ = min(£δc[2, i], 2.0) # nonlinear solver may cause this one to surpass 2
        c = min(£δc[3, i], 1.0) # nonlinear solver may cause this one to surpass 1 

        p = min(£ + c, 1.0) #  (all intact bridges)
        
        target[1:end-1, i] = f_tar(r, p, σ[i], c₀[i], δ, £, bins=model.system.num_tar_bins)
        target[end, i] = f_gas(r, g[i], σ[i], c₀[i])
    end
end


@jutul_secondary function update_variable!(target, v::SolidCharMassFraction, model, ftar_and_gas, ix)
    # NB: SolidCharMassFraction represents a reference value, as computed by the basic CPD model
    @inbounds for i in ix
        target[i] = 1 - sum(ftar_and_gas[:, i]) # last entry is gas, so this represents 1 - ftar - fgas
    end
end

# specifying default reaction rate parameters for bridge-breaking, light gas
# formation and cross-linking: [preexponential factor, activation energy, variance]
# note: energy unit is in calories
function Jutul.default_values(model, var::ArrheniusReactionRate{M}) where M
    if M == :b_rate
        return repeat([2.6e15, 55400.0, 1800.0], 1, number_of_entities(model, var))
    elseif M == :g_rate
        #return repeat([3.0e15, 69000.0, 8100.0], 1, number_of_entities(model, var))
        return repeat([3.0e15, 60000.0, 8100.0], 1, number_of_entities(model, var))
    elseif M == :cross_rate
        return repeat([3.0e15, 65000.0, 0.0], 1, number_of_entities(model, var))
    end
    error("No default reaction rate values for the specified parameter");
end

function Jutul.default_values(model, var::FixedReactionRate{:ρ_rate})
    return 0.9
end

# specify default parameters for flash calculation (α, β, γ)
# The defaults are here given with with respect to SI units (Pascals, kg/mol),
# in contrast to the values used in related literature which are in relation to
# atm and g/mol ([87058.0, 299.0, 0.5903])
Jutul.default_values(model, var::CPDFlashParams) =
    repeat([8.82115e9, 299 * 10^(3*0.5903), 0.5903] , 1, number_of_entities(model, var)) 

# Secondary variables
struct NonvaporInternalEnergy <: Jutul.ScalarVariable end # char and metaplast
struct VaporInternalEnergy <: JutulDarcy.PhaseVariables end
struct VaporEnthalpy <: JutulDarcy.PhaseVariables end 
struct TotalThermalEnergy <: Jutul.ScalarVariable end

# Parameters
struct CharHeatCapacity <: Jutul.ScalarVariable end # also used for metaplast
struct BulkDensity <: Jutul.ScalarVariable end
struct CharDensity <: Jutul.ScalarVariable end # also used for metaplast
struct VaporHeatCapacity <: JutulDarcy.PhaseVariables end

struct ThermalConductivity <: Jutul.ScalarVariable end
Jutul.associated_entity(::ThermalConductivity) = Faces()

struct TemperatureProfile <: Jutul.JutulForce
    Tfun::Function  # t -> T(t), temperature as function of time
end

struct VaporDensity <: JutulDarcy.PhaseMassDensities end

function Jutul.values_per_entity(model, rho::VaporDensity) 
    return Jutul.degrees_of_freedom_per_entity(model, rho)
end

function Jutul.degrees_of_freedom_per_entity(model, ::VaporDensity)
    return model.system.num_tar_bins + 1
end


function Jutul.default_values(model::JutulCPDModel, ::JutulDarcy.BulkVolume)
    return copy(model.data_domain.data[:volumes][1]) # cell volumes
end

struct RemainingCharVolume <: Jutul.ScalarVariable end # volume of the char remaining in the cell

# methods to update secondary variables

@jutul_secondary function update_variable(U, var::RemainingCharVolume, model,
                                          BulkDensity, BulkVolume, CharDensity, fchar, ix)
    for i in ix
        remaining_char_mass = BulkDensity[i] * BulkVolume[i] * fchar[i]
        U[i] = remaining_char_mass / CharDensity[i]
    end
end

# The parameter name for vapor density is 'PhaseMassDensities', in line with
# usage in JutulDarcy.
@jutul_secondary function update_variable!(rho, density::VaporDensity, model, 
                                           Pressure, Temperature, ξvapor, FluidVolume, ix)
    for i in ix
        ρ = ξvapor[:, i] ./ FluidVolume[i]
        rho[:, i] .= ρ
    end
end

@jutul_secondary function update_variable!(U, var::VaporInternalEnergy, model,
                                           PhaseMassDensities, Temperature,
                                           VaporHeatCapacity, ix)
    for i in ix
        U[i] = VaporHeatCapacity[i] * sum(PhaseMassDensities[:, i]) * Temperature[i]
    end
end

@jutul_secondary function update_variable!(H, var::VaporEnthalpy, model,
                                           VaporInternalEnergy, Pressure, ix)
    for i in ix
        H[i] = VaporInternalEnergy[i] + Pressure[i]
    end
end

@jutul_secondary function update_variable!(U, var::NonvaporInternalEnergy, model,
                                           CharDensity, ξmetaplast, Temperature,
                                           RemainingCharVolume, CharHeatCapacity, ix)
    for i in ix
        # internal energy of char
        char_energy = CharHeatCapacity[i] * CharDensity[i] * Temperature[i]
        # internal energy of metaplast
        mplast_density = ξmetaplast[i] / RemainingCharVolume[i]
        metaplast_energy = CharHeatCapacity[i] * mplast_density * Temperature[i]
        
        U[i] = char_energy + metaplast_energy
    end
end

@jutul_secondary function update_variable!(U, var::TotalThermalEnergy, model,
                                           NonvaporInternalEnergy, RemainingCharVolume,
                                           VaporInternalEnergy, FluidVolume, ix)
    for i in ix
        U[i] = NonvaporInternalEnergy[i] * RemainingCharVolume[i] +
               VaporInternalEnergy[i] * FluidVolume[i]
    end
end

@jutul_secondary function update_variable!(U, var::TotalCellMass, model, fchar, RemainingCharVolume, CharDensity, ξ, ix)

    # Compute the total mass currently contained within the cell.  This equals
    # the mass remaining in the solid matrix, plus the mass of volatiles and gas
    # that remains within the cell.
    
    # The solid mass is computed as the remaining part after released gas and
    # tar has been subtracted.  The released gas and tar is computed by CPD, and
    # is here provided in unmodified form, i.e., the fraction of tar and gas
    # released based on the initial cell mass (not the current mass).

    for i in ix
        # The mass of the char only depends on the progression of the reference CPD reaction, 
        # regardless of how much vapor/liquid has been evacuted or not
        char_cell_mass = RemainingCharVolume[i] * CharDensity[i]

        nonchar_cell_mass = sum(ξ[:, i]) # sum of all detached mass (whether vapor or liquid)

        U[i] = char_cell_mass + nonchar_cell_mass
    end
    
end

@jutul_secondary function update_variable!(U, var::DetachedMass{:Metaplast}, model, ξ, ξvapor, ix)
    for i in ix
        metaplast = max.(ξ[:, i] - ξvapor[:, i], 0.0) # avoid potential negative values from numerical noise
        U[:, i] = metaplast
    end
end


@jutul_secondary function update_variable!(U, var::DetachedMass{:Vapor}, model, ma, mb, £δc, σ, 
                                           Temperature, Pressure, ξ, FluidVolume, flash_αβγ, ix)
    
    for i in ix

        r = mb[i] / ma[i]
        light_gas_weight = mb[i] / 2.0

        £ = £δc[1, i]
        δ = £δc[2, i]
        c = £δc[3, i]

        α = flash_αβγ[1, i]
        β = flash_αβγ[2, i]
        γ = flash_αβγ[3, i]
        
        molweights =
            vcat(binned_molecular_weights(ma[i], r, £ + c, σ[i], δ, £, model.system.num_tar_bins),
                 light_gas_weight)
        mplast, vapors = flash(ξ[:, i], molweights, Temperature[i], Pressure[i], α=α, β=β, γ=γ)

        U[:, i] = vapors
    end
end

function Jutul.default_parameter_values(data_domain, model, param::ThermalConductivity, symb)

    Uvapor = data_domain[:vapor_conductivity]
    Usolid = data_domain[:solid_conductivity] 

    poro = data_domain[:porosity]
    
    U = (Uvapor .* poro) .+ (Usolid .* (1 .- poro))
    
    g = physical_representation(data_domain)
    T = compute_face_trans(g, U)
    if any(x -> x < 0, T)
        c = count(x -> x < 0, T)
        @warn "$c negative transmissibilities detected."
    end
    return T
end
