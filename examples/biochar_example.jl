using Jutul
using JutulDarcy
using Plots
import PyPlot
using MAT

## Define the 'recombine' function, which we will use when zipping parameter
#  vectors further down
function recombine(a, b, c)
    @assert(length(a) == length(b) == length(c))
    N = length(a)
    @assert(N % 3 == 0)
    
    result = zeros(N)
    for i = 1:(N/3)
        result[3*(i-1) + 1] = a[3*(i-1) + 1]
        result[3*(i-1) + 2] = b[3*(i-1) + 2]
        result[3*(i-1) + 3] = c[3*(i-1) + 3]
    end
    return result
end


## ===========================================================================
#                 Setup pressure, heating and simulation time
# ============================================================================

P0 = 101325.0; # 1 atm pressure
Tfun = (t) -> 300.0 + max(20.0 * t / 60.0, 500.0); # heating rate: 20 K/min
duration = 55.0 * 60.0; # 55 minutes
num_timesteps = 10000;
exponent = 1.3;
timesteps = range(0, stop=1.0, length=num_steps+1).^exponent .* duration;
timesteps = diff(timesteps);

## ============================================================================
#                  Setup grid, domain and simulation model
# ============================================================================

# Define the domain.  The 'biochar' grid has each of its cell broken up into
# separate subcells for lignin, cellulose and hemicellulose.
G, bc = radial_test_domain("data/grids/G_1cm_40_biochar.mat", P0, Tfun);

# Define the system of equations to be solved
sys = JutulCPDSystem(num_tar_bins=20,
                     imposed_global_temperature=imposed_global_temperature);

# Define the model, which combines the domain and the equation system
model = SimulationModel(G, sys, context = DefaultContext());

# ----------------------------------------------------------------------------
## set up a parameter storage, with the default parameters for lignin, cellulose
# and hemicellulose.  The easiest way to do this at present is to generate
# separate parameter objects for each material type, and then fuse them by
# 'zipping' values.
prm = setup_parameters(model, default_biochar);
lignin = setup_parameters(model, default_lignin);
cellulose = setup_parameters(model, default_cellulose);
hemicellulose = setup_parameters(model, default_hemicellulose);

material_specific_parameters = [:b_rate, :g_rate, :ρ_rate,
                                :c₀, :p₀, :σ, :ma, :mb]

for par in material_specific_parameters
    prm[par] = recombine(lignin[par], cellulose[par], hemicellulose[par]);
end

# ----------------------------------------------------------------------------
## We set up the initial state.  Here too, we create separate states
# for each material, and fuse them to a single state in the end.
matstates = Dict();
CharDensity = Dict();
for material, key in zip([lignin, cellulose, hemicellulose],
                         [:lignin, :cellulose, :hemicellulose])
    matstates[key], CharDensity[key] =
        setup_state(model, material[:p₀][1], material[:c₀][1],
                    material[:ma][1], material[:mb][1], material[:σ][1],
                    material[:BulkDensity][1], P0, Tfun(0.0),
                    α=prm_defaults[:flash_αβγ][1],
                    β=prm_defaults[:flash_αβγ][2],
                    γ=prm_defaults[:flash_αβγ][3]);
end

state0 = Dict();
for k in keys(matstates[:lignin])
    state0[k] = recombine(matstates[:lignin][k], matstates[:cellulose][k],
                          matstates[:hemicellulose][k]);
end
state0[:CharDensity] = recombine(CharDensity[:lignin], CharDensity[:cellulose],
                                 CharDensity[:hemicellulose]);


# ----------------------------------------------------------------------------
## We set up the external 'forcing', which here consists of the boundary conditions
forces = setup_forces(model, sources=nothing, bc=bc);

## ============================================================================
#                           Run the simulation
# ============================================================================

# Set the proper tolerances for the simulation.  These should be low enough to
# ensure proper convergence, but not too low as this will increase computational
# time as well as leading to convergence issues for some types of convergence
# criteria.
tolerances = Dict((:default=>1e-8,
                   :pressure_equation=>1e-3, 
                   :energy_conservation=>1e-5,
                   :mass_conservation=>1e-9,
                   :£δc=>1e-12));

# Setup the simulator, with the model, the initial state and the parameters to
# be used
sim = Simulator(model, state0 = state0, parameters = prm);

# Run the simulation (this may take a while)
states, reports = simulate!(sim, timesteps, forces=forces, info_level=1,
                            max_residual=1e99, tolerances=tolerances,
                            max_timestep_cuts=20)

