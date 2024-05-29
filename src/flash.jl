# Copyright (c) 2023 Odd Andersen
# Copyright (c) 2023 SINTEF Digital

using Roots

export flash

# ----------------------------------------------------------------------------
"""
    flash(binned_component_masses, binned_weights, temperature, pressure)
   
This function performs a flash calculation. It takes in the following parameters:

# Arguments:
- `binned_component_masses` : tar mass binned by cluster size (from 1 to n)
- `binned_weights`    : molecular weights for each tar bin (kg / mol)
- `temperature`       : the temperature in Kelvin    
- `pressure`          : pressure in Pascal

# Optional parameters:

The vapor pressure is modeled by the formula `P = α * exp(-β M^γ/T)`, where
M is the molecular weight of the tar, and T is the temperature.  The factors
α, β, and γ are optional parameters.  The default values are:
- α = 8.82115e9 (in Pascal, which is equivalent to 87058 in atmospheres)
- β = 299    (in Kelvin)
- γ = 0.5903 (unitless)

# Returns
- `Lmass` : liquid mass for each tar bin
- `Vmass` : vapor mass for each tar bin
- `LmTot` : total liquid mass
- `VmTot` : total vapor mass

The flash calculation is based on the formula and arguments presented in 
"Chemical Percolation Model for Devolatilization. 3. Direct Use of E NMR 
Data To Predict Effects of Coal Type" (1992) by Fletcher and Kerstein.
"""
function flash(binned_component_masses,
               binned_weights,
               temperature,
               pressure;
               α = 8.82115e9,
               β = 299 * 10^(3*0.5903), # original β = 299 was for g/mol, not kg/mol
               γ = 0.5903)

    
    # handel degenerate case with no mass
    if sum(binned_component_masses) == 0
        dummy = zeros(length(binned_component_masses))
        return dummy, copy(dummy), 0.0, 0.0, 0.0
    end
    
    # avoid NaN issues when weights are 0 (which should only happen if the
    # corresponding bins are also empty).  
    weights = ifelse.(binned_weights .== 0, eps(), binned_weights)
    # Avoid NaN/complications if pressure is nonpositive
    pressure = max(pressure, eps())
    
    f = vec(binned_component_masses ./ weights) # moles of each tar bin
    F = sum(f)                                  # total moles
    z = f ./ F                                  # mole fractions

    # Compute the vapor pressure and the corresponding K_i for each tar bin
    # (where y_i = K_i * x_i, and y_i and x_i are the vapor and liquid
    # mole fractions, respectively)
    
    Pv = α * exp.(-β * weights.^γ ./ max(temperature, eps()))
    
    K = Pv ./ pressure
    
    # Compute the ratio V/F, where V is the combined vapor fraction and F is
    # the sum of vapor and liquid fractions.
    # (eps added to avoid divide-by-zero errors at V_over_F = 1 when K tends towards
    # zero)
    fun = (VoF, zz, Km1) -> sum(zz .* Km1 ./ (Km1 .* VoF .+ 1 .+ eps()))
    fun1 = (x) -> fun(x, z, K .- 1)

    V_over_F = (fun1(0.0) <= 0.0) ? 0.0 : # global pressure too high to be consistent Pv
               (fun1(1.0) >= 0.0) ? 1.0 : # global pressure too low -> evaporate all
                                find_zero_wrapper(fun, [z, K.-1], [0.0, 1.0])

    V = V_over_F * F # total vapor moles
    L = F - V        # total liquid moles

    # compute liquid and vapor compositions
    x = z ./ (V_over_F .* (K .- 1) .+ 1 .+ eps())
    y = K .* x 

    # compute liquid and vapor masses
    Lmass = L .* x .* weights # liquid mass for each tar bin
    Vmass = V .* y .* weights # vapor mass for each tar bin

    LmTot = sum(Lmass)        # total liquid mass
    VmTot = sum(Vmass)        # total vapor mass

    # return result
    return (Lmass, Vmass, LmTot, VmTot, V) # Liquid mass, vapor mass, total liquid mass,
                                           # total vapor mass, total vapor moles
end
