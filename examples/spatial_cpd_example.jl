using Jutul
using JutulDarcy
using  Plots
import PyPlot
using MAT

# ============================================================================
#             PARAMETERS PROVIDING OPTIONS TO RUNNING THE SCRIPT
# ============================================================================

# To turn off heat transport and impose a uniform temperature throughout
# the domain, set the following variable to `true`.
imposed_global_temperature = false

# The default is to use a radial grid, representing a 1D row of cells from
# the centre to the edge of a spherical particle with radius 1 cm.
# Se the variable to `false` to use a cartesian grid instead of a radial.
use_radial_grid = true

# ============================================================================
#                  SETUP THE DOMAIN AND THE SIMULATION MODEL
# ============================================================================

# setup the grid, `G` and the boundary conditions, `bc`.   (The default parameters
# are defined in the file `defaults.jl`).
if use_radial_grid
    G, bc = radial_test_domain("data/grids/G_1cm_200.mat", P0, Tfun,
                               prm_defaults[:Permeability], 
                               prm_defaults[:Porosity],
                               prm_defaults[:CharThermalConductivity], 
                               prm_defaults[:VaporThermalConductivity]);
else
    G, bc = simple_test_domain([200, 1, 1], [1.0, 0.01, 0.01], P0, Tfun,
                               prm_defaults[:Permeability], 
                               prm_defaults[:Porosity],
                               prm_defaults[:CharThermalConductivity], 
                               prm_defaults[:VaporThermalConductivity]);
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
prm = setup_parameters(model, prm_defaults)

# We now compute the state, using a utility function found in `tools.jl`.  This
# function will compute the initial state.  The nontrivial parts of this
# computation is: (1) ensure pressure equilibrium between outside and the cell
# pore space; (2) compute the char density (which differs from bulk density
# as the latter also includes the non-evaporated tar).
state0, char_density = biochar_cpd.setup_state(model, prm[:p₀][1], prm[:c₀][1], prm[:ma][1],
                                               prm[:mb][1], prm[:σ][1], prm[:BulkDensity],
                                               P0, Tfun(0));

# We add the char density to the parameter store
prm[:CharDensity] = char_density;

# We now define the timesteps to be used in the simulation.  We use a nonuniform
# timestep size with finer timesteps towards the start of the simulation, to better
# capture the effect of rapid heating at the interface when the tar particle is
# exposed to the external heat.
duration = 0.03 * 1000; # Simulate 30 seconds
num_steps = 1000;       # Use 1000 timesteps (a fine time resolution is required
                        # to properly capture the CPD process and ensure global
                        # convergence.  Note that if timesteps are too long, the
                        # simulator will adapt by cutting them, but we ideally
                        # want to avoid this as it may increase total simulation
                        # computational time).  
exponent = 1.7; 
timesteps = range(0, stop=1.0, length=num_steps+1).^exponent .* duration;
timesteps = diff(timesteps);

# We set the external 'forcing', which consists of the boundary conditions, and
# optionally an imposed global temperature profile (if
# `imposed_global_temperature` is set to `true`).
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
tolerances = Dict((:default=>1e-8, #1e-9
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
