export prm_defaults

# Suggested default parameter values for the spatial CPD model

# Note that in line with usage in literature, the energies used to compute
# reaction rates (b_rate, g_rate, cross_rate and ρ_rate) are given in cal/mol,
# not J/mol!
prm_defaults = Dict([:c₀                          => 6.78e-2,                   #0.0,
                     :p₀                          => 0.53957,                   #0.61,
                     :σ                           => 3.98,           # nb: this is just σ, not σ+1
                     :mb                          => 67.095 / 1000,  # bridge mass (kg / mol)
                     :ma                          => 1.917e2 / 1000, # (mean) site mass (kg / mol)
                     :b_rate                      => [2.6e15, 55400.0, 1800.0], # A, E and variance (nb: E, σ is in cal/mol)
                     :g_rate                      => [3.0e15, 60000.0, 8100.0], # A, E and variance (nb: E, σ is in cal/mol)
                     :cross_rate                  => [3.0e15, 65000.0, 0.0],    # A, E and variance (nb: E, σ is in cal/mol)
                     :ρ_rate                      => 0.9,
                     :flash_αβγ                   => [8.82115e9, 299 * 10^(3*0.5903), 0.5903], # flash parameters (α, β, γ)
                     :Permeability                => 9.87e-13 * 1e-2, # 10 mD (m2)
                     :Porosity                    => 0.1,             # 10% porosity
                     :VaporHeatCapacity           => 1500.0,          # 1500.0, # J / (kg K)
                     :PhaseViscosities            => 1.5e-5,
                     :VaporThermalConductivity    => 1.0e-1,          # 
                     :CharThermalConductivity     => 0.2,             # 0.5e-1, # also covers metaplast
                     :BulkDensity                 => 350.0,           # kg / m3 (includes initial char, tar and gas)
                     :CharHeatCapacity            => 2000.0])         # J / (kg K) , also covers metaplast


# Note that in the literature, the flash parameters are typically given with
# respect to atmospheres for pressure and g/mol for molecular weights, whereas
# the defaults we use above are for Pa and kg/mol.
# In the original units, the vlaue for flash_αβγ above would be [87058.0, 299.0, 0.5903].

kcal2joule = 4184.0

default_biochar = Dict([:flash_αβγ => [8.82115e9, 299 * 10^(3*0.5903), 0.5903], # flash parameters (α, β, γ)
                       :Permeability => 3.16e-15, # geometric average of 1e-16 (wood) and 1e-13 (coal) (m2)
                       :Porosity => 0.66,
                       :VaporHeatCapacity => 1500.0, # J / (kg K)
                       :PhaseViscosities => 1.5e-5,
                       :VaporThermalConductivity => 0.25, # W/MK
                       :CharThermalConductivity => 0.25,  # W/MK
                       :BulkDensity => 500, # kg/m3 (Birch)
                       :CharHeatCapacity => 1750]) # J / (kg K) 
                       
default_lignin = Dict([:c₀ => 0.0,
                       :p₀ => 0.71,
                       :σ => 2.5,
                       :mb => 78.0/1000.0,
                       :ma => 208/1000.0,
                       :b_rate => [7.0e16, 55400.0,  500.0], # A, E (cal/mol) and variance (cal/mol)
                       :g_rate => [2.3e19, 69000.0, 2600.0],
                       :ρ_rate => 1.7,
                       :cross_rate => [3.0e15, 65000, 0.0]]);

default_cellulose = Dict([:c₀ => 0.0,
                          :p₀ => 1.0,
                          :σ => 2.0,
                          :mb => 45.4/1000.0,
                          :ma => 81.0/1000.0,
                          :b_rate => [2.0e16, 55400.0, 4100.0],
                          :g_rate => [3.0e15, 61200.0, 8100.0],
                          :ρ_rate => 100.0,
                          :cross_rate => [3.0e15, 65000, 0.0]]);

default_hemicellulose = Dict([:c₀ => 0.0,
                              :p₀ => 1.0,
                              :σ => 2.0,
                              :mb => 43.0/1000.0,
                              :ma => 77.5/1000.0,
                              :b_rate => [1.2e20, 51500.0, 100.0],
                              :g_rate => [3.0e15, 38200.0, 5000.0]],
                              :ρ_rate => 1.08,
                              :cross_rate => [3.0e15, 65000, 0.0]]);
