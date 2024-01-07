# Copyright (c) 2023 Odd Andersen
# Copyright (c) 2023 SINTEF Digital

using Jutul
using JutulDarcy
using CPDSpatial
using GLMakie

# ============================================================================
#             PARAMETERS PROVIDING OPTIONS TO RUNNING THE SCRIPT
# ============================================================================

# To turn off heat transport and impose a uniform temperature throughout
# the domain, set the following variable to `true`.
imposed_global_temperature = false;

# The default is to use a radial grid, representing a 1D row of cells from
# the centre to the edge of a spherical particle with radius 1 cm.
# Se the variable to `false` to use a cartesian grid instead of a radial.
use_radial_grid = true;

# We use atmospheric pressure, and make temperature a function of time,
# rising quickly from 300 K to a maximum of 1500 K.
P0 = 101325.0 # one atmosphere pressure, in Pascal
Tfun = (t) -> 300.0 + 1500.0 * min(1.0, t/1.5e-2);
#Tfun = (t) -> 300.0 + 1500.0 * min(1.0, t/15) # uncomment for a slow rise

# simulated duration
duration = 30.0; # in seconds
#duration = 45.0 # in seconds  # uncomment for a slow temperature rise

# ============================================================================
#                  SETUP THE DOMAIN AND THE SIMULATION MODEL
# ============================================================================

# setup the grid, `G` and the boundary conditions, `bc`.   (The default parameters
# are defined in the file `defaults.jl`).
if use_radial_grid
    G, bc = radial_test_domain("data/grids/G_1cm_100.mat", P0, Tfun, 
                               prm_defaults[:Permeability], 
                               prm_defaults[:Porosity],
                               prm_defaults[:CharThermalConductivity], 
                               prm_defaults[:VaporThermalConductivity]);
else
    G, bc = simple_test_domain((200, 1, 1), (0.01, 0.01, 0.01), P0, Tfun,
                               prm_defaults[:Permeability], 
                               prm_defaults[:Porosity],
                               prm_defaults[:CharThermalConductivity], 
                               prm_defaults[:VaporThermalConductivity],
                               isolated=[true, false, true, true, true, true]);
end

# Define the system of equations to be solved.  This is done by creating a
# JutulCPDSystem.  Currently, there are only two options to specify, namely the
# number of tar bins (we use 20, which is about the same as in the original CPD
# code), and whether to include internal heat transport or impose a global heat
# profile.
sys = JutulCPDSystem(num_tar_bins=20,
                     imposed_global_temperature=imposed_global_temperature);

# We now define the model, which consists of the system of equations `sys`, as
# well as the physical domain represented by the grid `G` and the boundary
# conditions `bc`.
model = SimulationModel(G, sys, context = DefaultContext());

# Set up a parameter storage for the model, with the default parameters
prm = setup_parameters(model, prm_defaults);

# We now compute the state, using a utility function found in `tools.jl`.  This
# function will compute the initial state.  The nontrivial parts of this
# computation is: (1) ensure pressure equilibrium between outside and the cell
# pore space; (2) compute the char density (which differs from bulk density
# as the latter also includes the non-evaporated tar).
state0 = CPDSpatial.setup_state(model, prm[:p₀][1], prm[:c₀][1], prm[:ma][1],
                                prm[:mb][1], prm[:σ][1], prm[:BulkDensity],
                                P0, Tfun(0),
                                α=prm_defaults[:flash_αβγ][1],
                                β=prm_defaults[:flash_αβγ][2],
                                γ=prm_defaults[:flash_αβγ][3]);

# We now define the timesteps to be used in the simulation.  We use a nonuniform
# timestep size with finer timesteps towards the start of the simulation, to better
# capture the effect of rapid heating at the interface when the tar particle is
# exposed to the external heat.
num_steps = 500;       
exponent = 1.7; 
timesteps = range(0, stop=1.0, length=num_steps+1).^exponent .* duration;
timesteps = diff(timesteps);

# We set the external 'forcing', which consists of the boundary conditions, and
# optionally an imposed global temperature profile (if `imposed_global_temperature`
# is set to `true`).
forces = setup_forces(model,
                      sources = imposed_global_temperature ? TemperatureProfile(Tfun) : nothing,
                      bc = bc);


# ============================================================================
#                  RUN THE SIMULATION
# ============================================================================

# Set the proper tolerances for the simulation.  These should be low enough to
# ensure proper convergence, but not too low as this will increase computational
# time as well as leading to convergence issues for some types of convergence
# criteria.
tolerances = Dict((:default=>1e-8,
                   :pressure_equation=>1e-3, 
                   :energy_conservation=>1e-5,
                   :mass_conservation=>1e-9,
                   :£δc=>1e-12))

# Setup the simulator, with the model, the initial state and the parameters to
# be used
sim = Simulator(model, state0 = state0, parameters = prm);

# The following lines allow specifying and applying an iterative linear solver.
# This is only useful if the number of cells becomes large, and is not needed
# for our modest-sized test grid, so the code is only provided here as reference.

# linear_solver = GenericKrylov(:gmres, max_iterations=3000,
#                               verbose=false, relative_tolerance=1e-6,
#                               preconditioner = ILUZeroPreconditioner(left=true, right=false))
linear_solver = nothing

# Run the simulation (this may take a while).  Note that if you run into
# convergence issues, it may be informative to specify a higher `info_level` (up
# to 4 possible).
states, reports = simulate!(sim, timesteps, forces=forces, info_level=1,
                            max_residual=1e99, tolerances=tolerances,
                            max_timestep_cuts=20,
                            linear_solver=linear_solver)

# ============================================================================
#                  ANALYSE THE RESULTS
# ============================================================================

cumtime = cumsum(timesteps); # vector with the exact time for each timestep
init_mass = prm[:BulkDensity]' * G.data[:volumes][1]; # total initial mass
N = length(states)

# compute tar and light gas yield over time, as well as the quantity of
# metaplast present in the material
lgas, tar, mplast = compute_yield_curves(bc, states, timesteps);

# compute reattached metaplast
reattached = compute_reattached_metaplast(states, timesteps);

# Plot the evolution of total char, gas and tar over time
lgasfrac = cumsum(lgas)./ init_mass;
tarfrac = cumsum(tar)./ init_mass;
mplastfrac = mplast ./ init_mass;
charfrac = 1 .- lgasfrac .- tarfrac;

f = GLMakie.Figure()
ax = GLMakie.Axis(f[1,1],
                  title="",
                  xlabel="Time",
                  ylabel="Fraction")
GLMakie.lines!(ax, cumtime[1:N], lgasfrac, label="light gas")
GLMakie.lines!(ax, cumtime[1:N], tarfrac, label="volatile tar")
GLMakie.lines!(ax, cumtime[1:N], mplastfrac, label="metaplast")
GLMakie.lines!(ax, cumtime[1:N], charfrac, label="char")
GLMakie.lines!(ax, cumtime[1:N],
               cumsum([sum(x) for x in reattached[1:N]] ./ init_mass),
               label="reattached metaplast")
GLMakie.axislegend()
display(GLMakie.Screen(), f)

# We can create 2D arrays representing selected variables in space and time, e.g.
pmat = hcat([x[:Pressure] for x in states]...);
tmat = hcat([x[:Temperature] for x in states]...);
lgmat = hcat([x[:ξ][end,:] for x in states]...);
lgdens = hcat([x[:ξ][end,:]./ G[:volumes] for x in states]...);
liqdens = hcat([sum(x[:ξ] - x[:ξvapor], dims=1)[:] for x in states]...) ./ G[:volumes];
attached = hcat(reattached...) ./ G[:volumes];


# Plot a surface expressing the evolution of pressure in space and time.  
fig, ax1 = GLMakie.surface(G[:cell_centroids][1,:], cumtime[1:N],
                           reverse(pmat, dims=1)./P0, axis=(type=Axis3,))
ax1.title="Pressure evolution"
ax1.xlabel="Distance from (spherical) particle center"
ax1.ylabel="Time (ms)"
ax1.zlabel="Pressure (atm)"
display(GLMakie.Screen(), fig)

# Likewise, we can plot a surface expressing temperature in space and time.
fig2, ax2 = GLMakie.surface(G[:cell_centroids][1,:], cumtime[1:N],
                           reverse(tmat, dims=1), axis=(type=Axis3,))
ax2.title="Temperature evolution"
ax2.xlabel="Distance from (spherical) particle center"
ax2.ylabel="Time (ms)"
ax2.zlabel="Temperature (K)"
display(GLMakie.Screen(), fig2)


# averaging light gas density over the three materials
fig3, ax3 = GLMakie.surface(G[:cell_centroids][1, :], cumtime[1:N],
                            reverse(lgdens, dims=1), axis=(type=Axis3,))
ax3.title="Light gas density"
ax3.xlabel="Distance from (spherical) particle center"
ax3.ylabel="Time (ms)"
ax3.zlabel="Density kg/m^3"
display(GLMakie.Screen(), fig3)

