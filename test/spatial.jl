using Jutul
using JutulDarcy
using CPDSpatial

function spatial_simulation()
    # if simulation runs through without issues, the test is considered a success

    P0 = 101325.0 # one atmosphere pressure, in Pascal
    Tfun = (t) -> 300.0 + 1500.0 * min(1.0, t/1.5e-2); # temperature profile
    # simulated duration
    duration = 3.0; # in seconds

    # define simple cartesian domain
    G, bc = simple_test_domain((20, 1, 1), (0.01, 0.01, 0.01), P0, Tfun,
                               prm_defaults[:Permeability], 
                               prm_defaults[:Porosity],
                               prm_defaults[:CharThermalConductivity], 
                               prm_defaults[:VaporThermalConductivity],
                               isolated=[true, false, true, true, true, true])

    # define system
    sys = JutulCPDSystem(num_tar_bins=20, imposed_global_temperature=false);

    # define model
    model = SimulationModel(G, sys, context = DefaultContext());

    # use default parameters
    prm = setup_parameters(model, prm_defaults);

    # setup initial state
    state0 = CPDSpatial.setup_state(model, prm[:p₀][1], prm[:c₀][1], prm[:ma][1],
                                prm[:mb][1], prm[:σ][1], prm[:BulkDensity],
                                P0, Tfun(0),
                                α=prm_defaults[:flash_αβγ][1],
                                β=prm_defaults[:flash_αβγ][2],
                                    γ=prm_defaults[:flash_αβγ][3])

    num_steps = 50;       
    exponent = 1.7; 
    timesteps = range(0, stop=1.0, length=num_steps+1).^exponent .* duration;
    timesteps = diff(timesteps);

    # define boundary conditions/forcing terms
    forces = setup_forces(model, sources = TemperatureProfile(Tfun), bc = bc);

    # set tolerances
    tolerances = Dict((:default=>1e-8,
                   :pressure_equation=>1e-3, 
                   :energy_conservation=>1e-5,
                   :mass_conservation=>1e-9,
                       :£δc=>1e-12))

    # Setup the simulator, with the model, the initial state and the parameters to
    # be used
    sim = Simulator(model, state0 = state0, parameters = prm);
    
    # run simulation
    states, reports = simulate!(sim, timesteps, forces=forces, info_level=1,
                            max_residual=1e99, tolerances=tolerances,
                            max_timestep_cuts=20)

    # if everything went through without issues, we report success
    return length(states) == num_steps
end


@testset "Spatial simulation" begin
    @test spatial_simulation()
end
    
