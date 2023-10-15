export prm_defaults

# Suggested default parameter values for the spatial CPD model

prm_defaults = Dict([:c₀                          => 6.78e-2,                   #0.0,
                     :p₀                          => 0.53957,                   #0.61,
                     :σ                           => 3.98,           # nb: this is just σ, not σ+1
                     :mb                          => 67.095 / 1000,  # bridge mass (kg / mol)
                     :ma                          => 1.917e2 / 1000, # (mean) site mass (kg / mol)
                     :b_rate                      => [2.6e15, 55400.0, 1800.0], # A, E and variance
                     :g_rate                      => [3.0e15, 60000.0, 8100.0], # A, E and variance
                     :cross_rate                  => [3.0e15, 65000.0, 0.0],    # A, E and variance
                     :ρ_rate                      => 0.9,
                     :flash_αβγ                   => [87058.0, 299.0, 0.5903],  # flash parameters (α, β, γ)
                     :Permeability                => 9.87e-13 * 1e-2, # 10 mD (m2)
                     :Porosity                    => 0.1,      # 10% porosity
                     :VaporHeatCapacity           => 1500.0, #1500.0,  # J / (kg K)
                     :PhaseViscosities            => 1.5e-5,
                     :VaporThermalConductivity    => 1.0e-1, # 
                     :CharThermalConductivity     => 0.2,#0.5e-1,  # also covers metaplast
                     :BulkDensity                 => 350.0,   # kg / m3 (includes initial char, tar and gas)
                     :CharHeatCapacity            => 2000.0]) # J / (kg K) , also covers metaplast
