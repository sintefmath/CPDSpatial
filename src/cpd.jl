# Copyright (c) 2023 Odd Andersen
# Copyright (c) 2023 SINTEF Digital

export ReactionRateParams, MaterialParams, cpd

using DifferentialEquations
using SpecialFunctions


"""
    struct ReactionRateParams

A struct representing parameters for a reaction rate in Arrhenius form with
distributed activation energy.

# Fields
- `A::Float64`: Preexponential factor.
- `E::Float64`: Activation energy.
- `σ::Float64`: Variation in the activation energy.

"""
struct ReactionRateParams
    A::Float64 
    E::Float64 
    σ::Float64
end

function Base.show(io::IO, rrp::ReactionRateParams)
    println("Pre-exponential factor (A):      ", rrp.A)
    println("Activation energy (E):           ", rrp.E)
    println("Activation energy std. dev. (σ): ", rrp.σ)
end

# ----------------------------------------------------------------------------

"""
    struct CoalMaterialParams

A struct representing coal/biomass-specific material parameters.

# Fields
- `σp1::Float64`: Coordination number.
- `p₀::Float64`: Fraction of intact bridges (incl. charred bridges).
- `c₀::Float64`: Fraction of charred bridges (≤ p₀)
- `r::Float64`: Ratio of bridge mass `mb` to average site mass `ma`.
- `ma::Float64`: Average site mass.
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

function Base.show(io::IO, mpar::MaterialParams)
    println("Coordination number plus 1 (σ+1):       ", mpar.σp1)
    println("Fraction of intact bridges (p₀):        ", mpar.p₀,
            " (NB: incl. charred bridges)")
    println("Fraction of charred bridges (c₀):       ", mpar.c₀, "  (NB: c₀ ≤ p₀)")
    println("Ratio of bridge to site mass (r=mb/ma): ", mpar.r)
    println("Average site mass (ma):                 ", mpar.ma)
end


# ----------------------------------------------------------------------------
"""
    cpd(AEσb, AEσg, AEσV, mpar, duration, Tfun; metaplast_model, Pfun, num_tar_bins, max_tstep)

This function performs a cpd simulation. It takes in the following parameters:

# Arguments
- `AEσb::ReactionRateParams`: Reaction rate parameters for the bridge breaking reaction.
- `AEσg::ReactionRateParams`: Reaction rate parameters for the formation of light gas.
- `AEσρ::ReactionRateParams`: Reaction rate parameters representing the ratio of side chain formation 
                              to char bridge formation.
- `mpar::MaterialParams`: Material-specific parameters for coal or biomass.
- `duration::Float64`: Duration of the process in seconds.
- `Tfun::Function`: A function representing the temperature as a function of time.
- `metaplast_model::Symbol`: Should be one of the following:
    - `:none` - use the basic CPD model without metaplast
    - `:original` - use the metaplast model from the original CPD paper (Fletcher, 1992)
    - `:modified` - use a modified model that aims for mass conservation (experimental)
- `num_tar_bins::Integer` - Number of tar bins (if a metaplast model is used)
- `max_tstep::Float64` - maximum timestep to use in the ODE integration (default: `Inf`)

# Returns
- A tuple with result values.  The following items (vectors) are always returned:
    - `time`: timepoints
    - `£vec`: fraction of labile bridges as function of time
    - `δvec`: fraction of side chains as function of time
    - `cvec`: fraction of char bridges as function of time
    - `g1`: light gas from side chain conversion: 2 * (1-p) - δ  
    - `g2`: light gas from char bridge formation: 2 * (c - c0)
    - `gvec`: stage of process for light gas formation
    - `fgas`: mass fraction of light gas
    - `ftar`: mass fraction of tar
    - `fchar`: mass fraction of tar

- If a metaplast model is used, the following additional items (vectors) are also returned:
    - `temp`: temperature
    - `mplast_bins`: metaplast mass fraction (binned, last bin represents light gas)
    - `fmetaplast`: metaplast  mass fraction (summed over bins, excluding light gas)
    - `fcross`: crosslinking mass fraction
"""
function cpd(AEσb::ReactionRateParams,     # bridge breaking reaction: £ →£⋆
             AEσg::ReactionRateParams,     # formation of light gas: δ → g₁
             AEσρ::ReactionRateParams,     # ratio of side chain formation to
                                           # char bridge formation: kδ/kc
             mpar::MaterialParams,         # coal/biomass-specific params
             duration::Float64,            # in seconds
             Tfun::Function;               # T(t): temperature as function of time
             metaplast_model = :original,  # if :none, ignore metaplast model
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
    #sol = solve(prob, alg_hints=[:stiff], dtmax=max_tstep)
    sol = solve(prob, Rodas4(), alg_hints=[:stiff], dtmax=max_tstep)

    # Compute mass fractions
    £vec = [x[1] for x in sol.u]     # fraction of labile bridges as function of time
    δvec = [x[2] for x in sol.u]     # fraction of side chains as function of time
    cvec = [x[3] for x in sol.u]     # fraction of char bridges as function of time
    gvec = [g(u) for u in sol.u]     # stage of process for light gas formation

    gvec = max.(gvec, 0.0); # avoid negative values from integration errors
    pvec = max.(min.(£vec + cvec, 1.0), 0.0); # fraction of intact bridges as function of time
                                              # (use min and max to ward off roundoff
                                              # errors outside bounds) 
    input = (mpar.r, mpar.σp1-1, mpar.c₀, δvec, £vec, gvec, pvec, cvec,
             sol.t, g1.(sol.u), g2.(sol.u))

    if  metaplast_model == :none 
        return basic_percolation_model(input...)
    elseif metaplast_model == :original
        return metaplast_percolation_model_orig(Tfun, Pfun, mpar.ma, input...,  num_bins=num_tar_bins)
    elseif metaplast_model == :modified
        return metaplast_percolation_model_modif(Tfun, Pfun, mpar.ma, input...,num_bins=num_tar_bins)
    end
    error("Inexistant metaplast model specified")
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
function metaplast_percolation_model_modif(Tfun, Pfun, ma, r, σ, c₀, δvec, £vec, gvec,
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
    f_gas_ref_last = 0.0;
    f_tar_ref_last = zeros(num_bins);
    for ix = 1:num_tsteps
        cur_ix = ix;
        prev_ix = max(ix-1, 1); # first iteration uses zero-initialized arrays
                                # to refer to previous timestep
        
        # compute temperature and pressure, current timestep
        T[cur_ix] = Tfun(time[cur_ix])
        P = Pfun(time[cur_ix])

        evacuated_tar = ftar[prev_ix]
        remaining_mass = 1.0 - evacuated_tar
        if remaining_mass < 0.0
            remaining_mass = 0.0
            @warn "Remaining mass negative."
        end

        # computing percolation reference values, unadjusted for evacuated mass
        f_gas_ref = f_gas(r, gvec[cur_ix], σ, c₀)
        f_gas_ref = max(f_gas_ref, f_gas_ref_last) # current reference fgas cannot be less than previous,
                                                   # though it may happen from integration error in gvec
        
        f_tar_ref = f_tar(r, pvec[cur_ix], σ, c₀, δvec[cur_ix], £vec[cur_ix], bins=num_bins)
        
        # Compute light gas mass fraction, adjusted for evacuated metaplast
        # (which does not contribute to light gas generation).
        d_gas = (f_gas_ref - f_gas_ref_last) * remaining_mass # light gas produced this timestep,
                                                              # adjusted for remaining mass
                                                                              
        fgas[cur_ix] =  fgas[prev_ix] + d_gas
        # fgas[cur_ix] =  f_gas_ref * remaining_mass # @@ This is in line with the original formulation,
        #                                             # and gives a result closer to the original code,
        #                                             # but the conceptually right approach does not seem
        #                                             # to be obvious, since fgas[i] now represents all the
        #                                             # gas produced up to and until timestep i, under the
        #                                             # assumtion that the total mass equaled _remaining_mass_
        #                                             # throughout the process, which should underestimate
        #                                             # the actual gas produced since the total mass at
        #                                             # the previous timestep was higher than _remaining_mass_.
        #                                             # d_gas = fgas[cur_ix] - fgas[prev_ix]

        # compute char mass consumed over this timestep (tar is immobile, so we do
        # not have to make adjustments for evacuated mass
        d_char = sum(f_tar_ref) + f_gas_ref - sum(f_tar_ref_last) - f_gas_ref_last

        #@assert -2*eps() <= d_char <= remaining_mass (may happen due to inaccuracies in numerical integration)
        d_char = max(d_char, 0.0) # prevent negative mass
        
        # compute corresponding molecular weights.  Divide by 1000 to get unit in kg/mol
        molweights = binned_molecular_weights(ma, r, pvec[cur_ix], σ, δvec[cur_ix], £vec[cur_ix], num_bins)

        # compute cross-linking
        dt = time[cur_ix] - time[prev_ix] # will be zero on first iteration, but that's ok

        cl_rate = acr * exp(-ecr / R / T[cur_ix]) # cross-linking rate
        clinkfrac = min(cl_rate * dt, 1) # fraction of metaplast that is re-attached
        fcross[cur_ix] = fcross[prev_ix] + clinkfrac * fmetaplast[prev_ix]
        
        # compute actual tar production/depletion this step, and adjust for lost
        # mass, considering that total production/depletion of tar and gas
        # balances the actual change in tar mass.
        d_tarmass = d_char - d_gas
        d_bins = f_tar_ref - f_tar_ref_last
        @assert d_tarmass >= sum(d_bins) - eps()
        
        # apply rescaling to make summed bins match actual produced tar, and add 
        # to existing metaplast
        remaining_metaplast = (1 .- clinkfrac) * mplast_prev_tstep[1:end-1]
        f = increment_bins(d_tarmass, d_bins, remaining_metaplast)
        
        # apply flash calculation
        (mplast_prev_tstep, vapors) = 
            flash([f..., d_gas], [molweights..., light_gas_weight], T[cur_ix], P)

        # the results above include the light gas, but we are only interested in the tar
        # so we do not include the final component of mplast_prev_tstep and vapors
        fmetaplast[cur_ix] = sum(mplast_prev_tstep[1:end-1])
        ftar[cur_ix] = sum(vapors[1:end-1]) + ftar[prev_ix]

        # compute char mass fraction as the remainder
        fchar[cur_ix] = 1 - ftar[cur_ix] - fgas[cur_ix]

        mplast_bins[cur_ix] = mplast_prev_tstep

        f_gas_ref_last = f_gas_ref
        f_tar_ref_last = f_tar_ref

    end

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

function increment_bins(d_tarmass, d_bins, remaining_metaplast)

    # - The total increment/decrement should be d_tarmass
    # - The distribution of the increment should be as close as possible to d_bins
    # - the sum of d_bins and remaining_metaplast should be nonzero in each bin

    if sum(remaining_metaplast) + d_tarmass + eps() < 0
        # @warn "Tar-mass to deduct is less than what is left.  Truncating."
        d_tarmass = -sum(remaining_metaplast)
    end

    # rescale d_bins
    tmp = sum(d_bins)
    if tmp * d_tarmass > 0.0
        # same sign, simple rescaling possible
        d_bins = d_bins ./ tmp * d_tarmass
    else
        # heuristic: add uniform vector
        alpha = (d_tarmass - sum(d_bins)) ./ length(d_bins)
        d_bins .+= alpha .* ones(length(d_bins))
    end
    @assert abs(sum(d_bins) - d_tarmass) <= 1e-10 * max(abs(d_tarmass), abs(tmp), eps())
    
    # add to metaplast, ensure nonzero components, and rescale if necessary
    f = remaining_metaplast + d_bins
    tot_remaining = sum(f)

    @assert tot_remaining + eps() >= 0.0

    f[f .< 0.0] .= 0.0
    tmp = sum(f)
    if tmp > 0.0
        f = f ./ tmp * tot_remaining
    end

    return f
end

# @@@ Old version (closer in structure to original code, but mass conservation is not ensured)
function metaplast_percolation_model_orig(Tfun, Pfun, ma, r, σ, c₀, δvec, £vec, gvec,
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
