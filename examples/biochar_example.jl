# Copyright (c) 2023 Odd Andersen
# Copyright (c) 2023 SINTEF Digital


using Jutul
using JutulDarcy
#using Plots
#import PyPlot
using MAT

## Define the 'recombine' function, which we will use when zipping parameter
#  vectors further down
function recombine(a::Matrix, b::Matrix, c::Matrix)
    N = size(a, 2)
    @assert(N % 3 == 0)
    
    result = zeros(size(a, 1), N)
    for i = 1:Int(N/3)
        result[:, 3*(i-1) + 1] = a[:, 3*(i-1) + 1]
        result[:, 3*(i-1) + 2] = b[:, 3*(i-1) + 2]
        result[:, 3*(i-1) + 3] = c[:, 3*(i-1) + 3]
    end
    return result
end

function recombine(a::Vector, b::Vector, c::Vector)
    @assert(length(a) == length(b) == length(c))
    N = length(a)
    @assert(N % 3 == 0)
    
    result = zeros(N)
    for i = 1:Int(N/3)
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
Tfun = (t) -> 300.0 + min(20.0 * t/60, 500.0); # heating rate: 20 K/min 

duration = 55.0 * 60.0; # 55 minutes
num_timesteps = 500; 

exponent = 1.0;
timesteps = range(0, stop=1.0, length=Int(ceil((num_timesteps+1)/2))).^exponent;
timesteps = [reverse(-timesteps[2:end])..., timesteps...] ./ 2 .+ 0.5;
timesteps = timesteps .* duration;

 timesteps = diff(timesteps);

## ============================================================================
#                  Setup grid, domain and simulation model
# ============================================================================

# Define the domain.  The 'biochar' grid has each of its cell broken up into
# separate subcells for lignin, cellulose and hemicellulose.
G, bc = radial_test_domain("data/grids/G_1cm_40_biochar.mat", P0, Tfun,
                           prm_defaults[:Permeability],
                           prm_defaults[:Porosity],
                           prm_defaults[:CharThermalConductivity],
                           prm_defaults[:VaporThermalConductivity]);

# Define the system of equations to be solved
sys = JutulCPDSystem(num_tar_bins=20, imposed_global_temperature=false); 
#sys = JutulCPDSystem(num_tar_bins=20, imposed_global_temperature=true); 

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
for (material, key) in zip([lignin, cellulose, hemicellulose],
                         [:lignin, :cellulose, :hemicellulose])
    matstates[key] = 
        CPDSpatial.setup_state(model, material[:p₀][1], material[:c₀][1],
                               material[:ma][1], material[:mb][1], material[:σ][1],
                               prm[:BulkDensity][1], P0, Tfun(0.0),
                               α=prm_defaults[:flash_αβγ][1],
                               β=prm_defaults[:flash_αβγ][2],
                               γ=prm_defaults[:flash_αβγ][3]);
end

state0 = Dict();
for k in keys(matstates[:lignin])
    state0[k] = recombine(matstates[:lignin][k], matstates[:cellulose][k],
                          matstates[:hemicellulose][k])
end

# ----------------------------------------------------------------------------
## We set up the external 'forcing', which here consists of the boundary conditions
forces = setup_forces(model, sources=nothing, bc=bc);
#forces = setup_forces(model, sources = TemperatureProfile(Tfun), bc=bc);

## ============================================================================
#                           Run the simulation
# ============================================================================

# Set the proper tolerances for the simulation.  These should be low enough to
# ensure proper convergence, but not too low as this will increase computational
# time as well as leading to convergence issues for some types of convergence
# criteria.
tolerances = Dict((:default=>1e-8,
                   :pressure_equation=>1e-3, 
                   :energy_conservation=>1e-3, #1e-5,
                   :ξ_conservation=>1e-3,#1e-6,
                   :total_mass_conservation=>1e-3,#1e-6,
                   :£δc=>1e-12));

# Setup the simulator, with the model, the initial state and the parameters to
# be used
sim = Simulator(model, state0 = state0, parameters = prm);

# Run the simulation (this may take a while)
states, reports = simulate!(sim, timesteps, forces=forces, info_level=1,
                            max_residual=1e99, tolerances=tolerances,
                            max_timestep_cuts=10)

# ============================================================================
#                  ANALYSE THE RESULTS
# ============================================================================

Plots.pyplot() # use the `pyplot` backend for plotting
cumtime = cumsum(timesteps); # vector with the exact time for each timestep
init_mass = prm[:BulkDensity]' * G.data[:volumes][1]; # total initial mass
N = length(states)
# utility function to wait for a keypress
function wait_for_key(info)
    print(stdout, info); read(stdin, 1); nothing;
end

# Post-processing of results

# compute tar and light gas yield over time, as well as the quantity of
# metaplast present in the material
lgas, tar, mplast = compute_yield_curves(bc, states, timesteps[1:N]);

# compute reattached metaplast
reattached = compute_reattached_metaplast(states, timesteps[1:N]);

# Plot the evolution of total char, gas and tar over time
lgasfrac = cumsum(lgas)./ init_mass;
tarfrac = cumsum(tar)./ init_mass;
mplastfrac = mplast ./ init_mass;
charfrac = 1 .- lgasfrac .- tarfrac;
plot(cumtime[1:N], lgasfrac, reuse=false, label="light gas")
plot!(cumtime[1:N], tarfrac, label = "volatile tar")
plot!(cumtime[1:N], mplastfrac, label = "metaplast")
plot!(cumtime[1:N], charfrac, label="char")
plot!(cumtime[1:N],
      cumsum([sum(x) for x in reattached[1:N]] ./ init_mass),
      label="reattached metaplast")

wait_for_key("Press any key to continue\n");

# We can create 2D arrays representing selected variables in space and time, e.g.
pmat = hcat([x[:Pressure] for x in states]...);
tmat = hcat([x[:Temperature] for x in states]...);
lgmat = hcat([x[:ξ][end,:] for x in states]...);
lgdens = hcat([x[:ξ][end,:]./ G[:volumes] for x in states]...);
liqdens = hcat([sum(x[:ξ] - x[:ξvapor], dims=1)[:] for x in states]...) ./ G[:volumes];
attached = hcat(reattached...) ./ G[:volumes];

# Plot a surface expressing the evolution of pressure in space and time.  Note
# that there are some spikes associated with the thermal shock on the boudnary
# at the start of the simulation, as well as when the last tar components
# are in the centre towards the end.  High-resolution in space and time is needed.
sfplot = Plots.surface(pmat', reuse=false, camera=(70, 30))
display(sfplot)

# Likewise, we can plot a surface expressing the light gas density in space
# and time.
sfplot2 = Plots.surface(lgdens', reuse=false, camera=(70, 30))
display(sfplot2)

# Results can also be saved as a matlab file, to benefit from its extensive
# plotting facilities:
using MAT
matwrite("result.mat", Dict("P"=>pmat, "T"=>tmat, "X"=>lgmat, "Xd"=>lgdens,
                            "attached" => attached, "Ld" => liqdens, "cumtime"=>cumtime))
