using MAT

# P0 and T0 represents pressure (in Pascal) and temperature (in Kelvin).
#`p0` represents the fraction of intact bridges (incl. char bridges).  `c0`
# represents the fraction of char bridges.
"""
Setup an initial state for use with the spatial CPD model.

The initial state assumes pressure equilibrium with the surroundings.  Pressure
inside the material is computed from the moles of gas present, using the ideal
gas law.  Light gas is added until this pressure is reached.

The initial gas and tar is computed from the material parameters using the CPD
model.  The parameters to be provided are:
- `p0`: initial fraction of intact bridges 
- `c0`: initial fraction of char bridges
- `ma`: mean molecular mass per site
- `mb`: molecular mass of bridge
- `σ`: coordination number

Other parameters are:
- `bulk_density`: bulk density of the material, which includes the char and
                  initial tar taken together.
- `P0`: initial pressure
- `T0`: initial temperature
- `model`: the CPD model to use, which also contains information about the grid, 
           including cell volumes and pore volumes.

Returns:
- `state`: the initial state
- `char_density`: the initial char density.  This will normally be close to the
   bulk density, but does not include the contribution of initial tar.
"""
function setup_state(model::JutulCPDModel, p0::Real, c0::Real,
                     ma::Real, mb::Real, σ::Real, bulk_density, P0::Real, T0::Real;
                     α = 8.82115e9,
                     β = 299 * 10^(3*0.5903), # original β = 299 was for g/mol, not kg/mol
                     γ = 0.5903)

    @assert(0.0 <= c0 <= p0 <= 1.0)
    £ = p0 - c0
    δ = 2.0 * (1.0 - p0)
    c = c0

    r = mb / ma
    ftar = f_tar(r, p0, σ, c0, δ, £, bins = model.system.num_tar_bins)

    cell_masses = model.data_domain[:volumes] .* bulk_density

    # compute initial amount of gas in the pore space, estimated as the amount
    # of light gas necessary to fill the pore volume space at the initial pressure.
    porevolume = model.data_domain[:volumes] .* model.data_domain[:porosity]

    R = 8.31446261815324 # J/(mol K) Gas constant
    n_guess = P0 * porevolume / (R * T0)  # PV = nRT

    tar_component_masses = binned_molecular_weights(ma, r, p0, σ, δ, £, model.system.num_tar_bins)
    light_gas_weight = mb / 2.0
    ξ_comp_weights = vcat(tar_component_masses, light_gas_weight)
    
    # determine number of moles of light gas needed to correspond to the prescribed pressure
    # the moles of gas determined by the flash calculation should equal the number of moles required
    # to reproduce the pressure
    
    vapor_moles = (gmass, cell) -> flash(vcat(ftar * cell_masses[cell], gmass),
                                         ξ_comp_weights, T0, P0, α=α, β=β, γ=γ)[5]

    gmassvec = [find_zero(gmass -> vapor_moles(gmass, c) * R * T0 / porevolume[c] - P0,
                          [0.0, 10 * n_guess[c] * light_gas_weight])
                for c in 1:length(cell_masses)]

    ξ0 = vcat([ftar...] * cell_masses', gmassvec')

    ξ0vapor = [flash(ξ0[:, c], ξ_comp_weights, T0, P0, α=α, β=β, γ=γ)[2] for c in 1:length(cell_masses)]
    ξ0vapor = hcat(ξ0vapor...) # convert into matrix

    state = Jutul.setup_state(model, Dict(:Pressure => P0, 
                                         :Temperature => T0,
                                         :£δc => [£, δ, c],
                                         :ξ => ξ0,
                                         :ξvapor => ξ0vapor))
    
    char_density = (cell_masses - sum(ξ0, dims=1)[:]) ./ model.data_domain[:volumes]
    
    return state, char_density
end

"""
    Post-processing function to compute mass fluxes across the internal faces of the grid between
    each timestep.
"""
function compute_internal_mass_fluxes(sim, state)

    model = sim.model
    
    # compute internal fluxes
    disc = model.equations[:mass_conservation].flow_discretization
    num_internal_faces = maximum(x.face for x in disc.conn_data)

    num_components = model.system.num_tar_bins + 1
    result = zeros(num_components, num_internal_faces)
    dt_dummy = 0.0 # formal argument, but unused
    
    for conn in disc.conn_data
        if conn.face_sign == 1
            face, face_sign, other, self = conn

            state_named_tuple =
                (Pressure = state[:Pressure],
                 PhaseViscosities = state[:PhaseViscosities],
                 Transmissibilities = sim.storage.parameters.Transmissibilities,
                 TwoPointGravityDifference = sim.storage.parameters.TwoPointGravityDifference,
                 PhaseMassDensities = state[:PhaseMassDensities])
            
            result[:, face] .= Jutul.face_flux!([], self, other, face, face_sign,
                                                model.equations[:mass_conservation],
                                                state_named_tuple, model, dt_dummy, disc)
        end
    end

    return result
end

"""
    Post-processing function to compute mass fluxes across the external faces 
    of the grid between each timestep, and summarize these as lightgas, tar and
    metaplast (the volatiles that remain inside the cells)
"""
function compute_yield_curves(bc, states, dt)

    released_volatiles = [compute_boundary_mass_flux(bc, s) for s in states]

    # compute output across boundaries
    lightgas = [sum([y[end] for y in x]) for x in released_volatiles] .* dt
    tar = [sum([sum(y[1:end-1]) for y in x]) for x in released_volatiles] .* dt

    # compute internal metaplast amount
    metaplast = [sum(x[:ξ]) - sum(x[:ξvapor]) for x in states]
    
    return lightgas, tar, metaplast
end


"""
    Post-processing function to compute the amount of metaplast that has been
    reattached to the char surface in each cell.
"""
function compute_reattached_metaplast(states, timesteps)

    [sum(min.(states[i][:MetaplastCrossLinkRate] .* timesteps[i],
              states[i][:ξ] - states[i][:ξvapor]),
         dims=1)[:]
     for i in 1:length(states)]
end


"""
    Create a simulation domain and boundary condition object corresponding to a
    radial test domain with the given permeability, porosity, and solid and vapor.
"""
function radial_test_domain(gridname, P0, Tfun, perm, poro, 
                            solid_conductivity, vapor_conductivity)

    matfile = matread(gridname)
    wrapmesh = MRSTWrapMesh(matfile["G"])
    G = reservoir_domain(UnstructuredMesh(wrapmesh),
                         permeability=perm,
                         porosity=poro,
                         solid_conductivity=solid_conductivity,
                         vapor_conductivity=vapor_conductivity)

    # setup boundary conditions
    bc_trans = compute_boundary_trans(G, :permeability)
    bc_hs_trans = compute_boundary_trans(G, :solid_conductivity)
    bc_hv_trans = compute_boundary_trans(G, :vapor_conductivity)

    num_bc = length(bc_trans)

    tot_htrans = poro * bc_hv_trans + (1-poro) * bc_hs_trans

    # there should be only one open face, the outermost one.  This is the one
    # with the largest x-coordinate
    outer_bface = argmax(G[:boundary_centroids][1,:]) 

    bc = [CPDFlowBoundaryCondition(G[:boundary_neighbors][outer_bface],
                                   P0, Tfun, bc_trans[outer_bface], tot_htrans[outer_bface])]

    return G, bc    
    
end

"""
    Create a simulation domain and boundary condition object corresponding to a
    simple cartesian prism test domain with the given permeability, porosity, 
    and solid and vapor.  Each boundary may be designated as "isolated" or "open".
    Heat and mass fluxes will only pass through open boundaries.
"""
function simple_test_domain(resolution, phys_dim, P0, Tfun,
                            perm, poro, 
                            solid_conductivity, vapor_conductivity;
                            isolated=[false, false, false, false, false, false])
                            
    # setup grid
    G = reservoir_domain(CartesianMesh(resolution, phys_dim),
                         permeability=perm,
                         porosity=poro,
                         solid_conductivity=solid_conductivity,
                         vapor_conductivity=vapor_conductivity)
    
    # setup boundary conditions
    bc_trans = compute_boundary_trans(G, :permeability)
    bc_hs_trans = compute_boundary_trans(G, :solid_conductivity)
    bc_hv_trans = compute_boundary_trans(G, :vapor_conductivity)

    num_bc = length(bc_trans)

    tot_htrans = poro * bc_hv_trans + (1-poro) * bc_hs_trans

    normals = G[:boundary_normals]
    for dim in 1:3
        if isolated[2*dim-1] || isolated[2*dim]
            isolated_faces = ((normals[dim, :] .< 0.0) .& isolated[2*dim-1]) .|
                             ((normals[dim, :] .> 0.0) .& isolated[2*dim])
            bc_trans[isolated_faces] .= 0.0
            tot_htrans[isolated_faces] .= 0.0
        end
    end
    
    bc = [CPDFlowBoundaryCondition(G[:boundary_neighbors][i],
                                   P0, Tfun,
                                   bc_trans[i], tot_htrans[i]) for i in 1:num_bc]
    
    return G, bc
end
