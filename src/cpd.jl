export ReactionRateParams, MaterialParams, cpd

using DifferentialEquations
using SpecialFunctions

"""
    struct ReactionRateParams

A struct representing parameters for a reaction rate in Arrhenius form with
distributed activation energy.

# Fields
- `A`: Preexponential factor.
- `E`: Activation energy.
- `σ`: Variation in the activation energy.

"""
struct ReactionRateParams
    A::Float64 
    E::Float64 
    σ::Float64
end
# ----------------------------------------------------------------------------

"""
    struct CoalMaterialParams

A struct representing coal/biomass-specific material parameters.

# Fields
- `σp1`: Coordination number.
- `p₀`: Fraction of intact bridges (incl. charred bridges).
- `c₀`: Fraction of charred bridges (≤ p₀)
- `r`: Ratio of bridge mass `mb` to average site mass `ma`.
- `ma`: Average site mass.
"""
struct MaterialParams
    σp1::Float64 # coordination number (σ+1)
    p₀::Float64  # fraction of intact bridges (incl. charred bridges)
    c₀::Float64  # fraction of charred bridges (<= p₀)
    r::Float64   # ratio of bridge mass `mb` to site mass `ma`
    ma::Float64  # average site mass (NB: kg/mol, not g/mol!)
end

# The following constructor is useful when running the basic model, where
# average site masses are irrelevant
MaterialParams(σp1, p₀, c₀, r) = MaterialParams(σp1, p₀, c₀, r, 1.0)

# ----------------------------------------------------------------------------
"""
cpd(AEσb, AEσg, AEσV, mpar, duration, Tfun; basic_model=false)

This function performs a cpd simulation. It takes in the following parameters:

- AEσb: Reaction rate parameters for the bridge breaking reaction.
- AEσg: Reaction rate parameters for the formation of light gas.
- AEσρ: Reaction rate parameters representing the ratio of side chain formation 
        to char bridge formation.
- mpar: Material-specific parameters for coal or biomass.
- duration: Duration of the process in seconds.
- Tfun: A function representing the temperature as a function of time.
- basic_model: If true, limit computations to basic model from initial (1988) paper
"""
function cpd(AEσb::ReactionRateParams,     # bridge breaking reaction: £ →£⋆
             AEσg::ReactionRateParams,     # formation of light gas: δ → g₁
             AEσρ::ReactionRateParams,     # ratio of side chain formation to
                                           # char bridge formation: kδ/kc
             mpar::MaterialParams,         # coal/biomass-specific params
             duration::Float64,            # in seconds
             Tfun::Function;               # T(t): temperature as function of time
             basic_model = false,          # if true, ignore metaplast model
             Pfun::Function = t -> 101325, # P(t): pressure as function of time, in Pascal,
                                           # defaults to constant 1 atm = 101325Pa.  
                                           # Pressure is not used in the basic model.
             num_tar_bins = 20,            # number of tar bins (if basic_model = false)
             max_tstep = Inf               # 
             )
    # Arrhenius rate law
    R = 1.9872036 # universal gas constant in cal/mol/K
    kfun(A, E, T) = A * exp.(-E ./ (R*T)) # Arrhenius rate law

    # initial fractions of labile bridges (£) and side chains (δ)
    (£₀, δ₀) = (mpar.p₀ - mpar.c₀, 2 * (1 - mpar.p₀)) # initial fractions

    # distributed activation energy
    E_activation(dfrac, E, σ) =
        E + sqrt(2) * σ * erfinv(clamp(2 * dfrac - 1, -1+sqrt(eps()), 1-sqrt(eps())))

    # initial condition vector
    u₀ = [£₀, δ₀, mpar.c₀]

    # ODE integration with variable activation energies
    g1(u) = 2 * (1 - u[1] - u[3]) - u[2] # light gas from side chain conversion: 2 * (1-p) - δ  
    g2(u) = 2 * (u[3] - mpar.c₀)         # light gas from char bridge formation: 2 * (c - c0)
    g(u) = g1(u) + g2(u)                 # total light gas
    gfrac(u) = g(u) / (2 * (1-mpar.c₀))  # stage of process for light gas formation
    lfrac(u) = 1 - u[1] / £₀             # stage of process for labile bridge breaking

    # Activation energy functions for bridge breaking, light gas formation, and
    # side chain formation
    kb = (u, t) -> kfun(AEσb.A, E_activation(lfrac(u), AEσb.E, AEσb.σ), Tfun(t))
    kg = (u, t) -> kfun(AEσg.A, E_activation(gfrac(u), AEσg.E, AEσg.σ), Tfun(t))
    ρ =  (u, t) -> kfun(AEσρ.A, AEσρ.E, Tfun(t)) # The CPD model does not take variance
                                                 # in activation energy into account, so we use
                                                 # AEσρ.E directly
    # Specifying the differential equation:
    # - Line 1 of matrix M: destruction of labile bridges
    # - Line 2 of matrix M: formation of side chains
    # - Line 3 of matrix M: formation of char bridges
    M = (u, t) -> [-kb(u, t)                             0           0;
                   2ρ(u, t) * kb(u, t) / ( ρ(u, t) + 1)  (-kg(u, t)) 0;
                   kb(u, t) / (ρ(u, t) + 1)              0           0]
    
    Du(u, p, t) = M(u, t) * u # note that p (which Julia expects to be parameters) is not used

    # Solve the ODE
    prob = ODEProblem(Du, u₀, (0.0, duration))
    sol = solve(prob, alg_hints=[:stiff], dtmax=max_tstep)

    # Compute mass fractions
    £vec = [x[1] for x in sol.u]     # fraction of labile bridges as function of time
    δvec = [x[2] for x in sol.u]     # fraction of side chains as function of time
    cvec = [x[3] for x in sol.u]     # fraction of char bridges as function of time
    gvec = [g(u) for u in sol.u]     # stage of process for light gas formation

    pvec = £vec + cvec; # fraction of intact bridges as function of time

    input = (mpar.r, mpar.σp1-1, mpar.c₀, δvec, £vec, gvec, pvec, cvec,
             sol.t, g1.(sol.u), g2.(sol.u))

    return basic_model ? basic_percolation_model(input...) :
                         metaplast_percolation_model(Tfun, Pfun, mpar.ma, input...,
                                                     num_bins=num_tar_bins)
        
end

# ----------------------------------------------------------------------------
"""
    basic_percolation_model(r, σ, c₀, δvec, £vec, gvec, pvec, cvec, time, g1, g2)
    
Compute the basic percolation model from the initial (1988) paper.
"""
function basic_percolation_model(r, σ, c₀, δvec, £vec, gvec, pvec, cvec, time, g1, g2)

    fgas = f_gas(r, gvec, σ, c₀);             # mass fraction of light gas
    ftar = f_tar(r, pvec, σ, c₀, δvec, £vec); # mass fraction of tar
    fchar = 1 .- fgas .- ftar;                # mass fraction of char
    
    return (time = time,
            £vec =  £vec,
            δvec = δvec,
            cvec = cvec,
            g1 = g1,
            g2 = g2,
            gvec = gvec,
            fgas = fgas,
            ftar = ftar,
            fchar = fchar)
end

# ----------------------------------------------------------------------------
"""
    metaplast_percolation_model(ma, r, σ, c₀, δvec, £vec, gvec, pvec, cvec, time, g1, g2)
    
Compute the metaplast percolation model from the 1992 paper.  
"""
function metaplast_percolation_model(Tfun, Pfun, ma, r, σ, c₀, δvec, £vec, gvec,
                                     pvec, cvec, time, g1, g2;
                                     num_bins=20,
                                     acr=3.0e15, # pre-exponential factor for crosslinking
                                     ecr=65000.0 # activation energy for crosslinking
                                     )
    R = 1.9872036 # universal gas constant in cal/mol/K
    num_tsteps = length(time)
    light_gas_weight = r * ma / 2.0 # light gas molecular weight
    
    # return variables
    fgas = zeros(num_tsteps)       # mass fraction of light gas
    ftar = zeros(num_tsteps)       # mass fraction of tar (in volatile form, directly evacuated)
    fchar = ones(num_tsteps)       # mass fraction of char
    fmetaplast = zeros(num_tsteps) # mass fraction of metaplast
    fcross = zeros(num_tsteps)     # mass fraction of crosslinking
    T = Tfun(0) * ones(num_tsteps) # temperature

    bins_prev_tstep = zeros(num_bins) # mass fraction of tar in each bin at previous time step
    mplast_prev_tstep = zeros(num_bins+1) # mass fraction of metaplast in each bin at previous
                                          # time step, and an extra bin to store liquid part of
                                          # light gas
    mplast_bins = Vector{Vector{Float64}}(undef, num_tsteps) # mplast_bins number for each bin at
                                                             # each time step.  For reporting purposes
                                                             # only
    for i = 2:num_tsteps
        # compute temperature and pressure, current timestep
        T[i] = Tfun(time[i])
        P = Pfun(time[i])
        # Compute light gas mass fraction, adjusted for evacuated metaplast
        # (which does not contribute to light gas generation).
        # @@ Note: this is in line with the treatment of fgas in the original code, but
        # appears to be inaccurate, as fgas[i] now represents all the gas produced up to
        # and until timestep i, considering that the evacuated gas was _always_ absent.
        # When we take the difference with fgas[i-1], we then substract a quantity that
        # was computed based on the assumption that the evacuated gas _at_that_timestep_
        # was always absent.
        evacuated_tar = ftar[i-1]
        fgas[i] = f_gas(r, gvec[i], σ, c₀) * (1 - evacuated_tar)
                                                
        # compute binned tar mass fraction with no adjustments 
        bins = f_tar(r, pvec[i], σ, c₀, δvec[i], £vec[i], bins=num_bins)

        # compute corresponding molecular weights.  Divide by 1000 to get unit in kg/mol
        molweights = binned_molecular_weights(ma, r, pvec[i], σ, δvec[i], £vec[i], num_bins)

        # compute cross-linking
        dt = time[i] - time[i-1]

        cl_rate = acr * exp(-ecr / R / T[i]) # cross-linking rate

        clinkfrac = min(cl_rate * dt, 1) # fraction of metaplast that is re-attached
        
        fcross[i] = fcross[i-1] + clinkfrac * fmetaplast[i-1]
        
        # adjust bins for metaplast  

        # @@ Is the accounting here sound? It seems that we want to remove the
        # component that has been "vented" so far.  The way this is done is to
        # remove all finite clusters that were computed at last step, and then
        # "adding back in" the metaplast component.  In other words:
        # new_tar = total_produced_tar - vented
        #         = total_produced_tar - (total_produced_tar_prev - metaplast_prev)
        # This, however, seems to neglect the fact that total_produced_tar would
        # also depend on the reduction in total mass resulting from earlier venting.
        # On the other hand, the vented components would not produce new tar, only
        # change the composition of the already-produced tar.
        
        d_bins = max.(bins .- bins_prev_tstep, 0.0)  # presumably, removing what has been "vented"
        bins_prev_tstep = copy(bins)
        
        f = d_bins .+ ((1 .- clinkfrac) .* mplast_prev_tstep[1:end-1]) # adding back in what has not been "vented"
        d_gas  = max(fgas[i] - fgas[i-1], 0.0) # light gas produced this timestep

        # apply flash calculation
        (mplast_prev_tstep, vapors) = 
            flash([f..., d_gas], [molweights..., light_gas_weight], T[i], P)

        # the results above include the light gas, but we are only interested in the tar
        # so we do not include the final component of mplast_prev_tstep and vapors
        fmetaplast[i] = sum(mplast_prev_tstep[1:end-1])
        ftar[i] = sum(vapors[1:end-1]) + ftar[i-1]

        # compute char mass fraction as the remainder
        fchar[i] = 1 - ftar[i] - fgas[i]

        mplast_bins[i] = copy(mplast_prev_tstep)

    end

    fmetaplast[1] = fmetaplast[2] # we didn't know the initial value when
                                  # initializing the array above
    
    return (time = time,
            temp = T,
            £vec =  £vec,
            δvec = δvec,
            cvec = cvec,
            g1 = g1,
            g2 = g2,
            gvec = gvec,
            fgas = fgas,
            ftar = ftar,
            fchar = fchar,
            mplast_bins = mplast_bins,
            fmetaplast = fmetaplast,
            fcross = fcross)
end
