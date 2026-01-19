using Markdown
using CPDSpatial
using Interpolations
using DelimitedFiles
using DataStructures

isCI = get(ENV, "CI", "false")

if isCI == "true"
    using CairoMakie
    Mk = CairoMakie
else
    using GLMakie
    Mk = GLMakie
end

export cpd_benchmarking

# ----------------------------------------------------------------------------
"""
    cpd_benchmarking(choice)

This script provides multiple examples of the `cpd()` function applied on
different materials.  Most of the examples are taken from literature, with 
which results can be compared and benchmarked.

This function takes a single argument, `choice`, which specify which example you
want to run.

# Arguments
- `choice::Symbol` - The different choices are:
    - `:cpdheat` - Run cpd with the parameters from the original MATLAB CPDHeat code 
    - `:fcompare` - Run cpd to reproduce example presented in the [original CPD Fortran code](http://www.et.byu.edu/~tom/cpd/cpd92/cpdfiles.html)
    - `:three_coals` - Basic cpd simulation results for three different coal types
    - `:heating_rate` - This example illustrates the effect of the heating rate on tar yield.  
    - `:biochar` - Running CPD on three main components of biomaterials: *lignin*, *cellulose* and *hemicellulose* (xylan), using a *modified* metaplast model.
   - `:biochar_modif` - Running CPD on three main components of biomaterials: *lignin*, *cellulose* and *hemicellulose* (xylan), using a *modified* metaplast model with different parameters.
"""
function cpd_benchmarking(choice=nothing)

    dispatch = OrderedDict(:cpdheat => cpdheat_compare,
                           :fcompare => cpd_fortran_compare,
                           :three_coals => cpd_three_coals,
                           :heating_rate => cpd_heating_rate,
                           :biochar => cpd_biochar_constituents,
                           :biochar_modif => cpd_biochar_consituents_modif)
    get(dispatch, choice, () -> helptext(dispatch, choice))()
end

# ----------------------------------------------------------------------------
"""
Running CPD on three main components of biomaterials:
- lignin
- cellulose
- hemicellulose (here: xylan)
This example uses the *modified* (and experimental) metaplast model to ensure 
mass conservation.  This is the same model that is used in the spatial examples.
"""
function cpd_biochar_constituents()

    # Define reaction rate parameters for the three biomaterials
    AEσb_lignin = ReactionRateParams(7.0e16, 55400.0, 500.0)
    AEσg_lignin = ReactionRateParams(2.3e19, 69000.0, 2600.0)
    AEσρ_lignin = ReactionRateParams(1.7, 0.0, 0.0)

    AEσb_cellulose = ReactionRateParams(2.0e16, 55400.0, 4100.0)
    AEσg_cellulose = ReactionRateParams(3.0e15, 61200.0, 8100.0)
    AEσρ_cellulose = ReactionRateParams(100, 0.0, 0.0)

    AEσb_xylan = ReactionRateParams(1.2e20, 51500.0, 100.0)
    AEσg_xylan = ReactionRateParams(3.0e15, 38200.0, 5000.0)
    AEσρ_xylan = ReactionRateParams(100, 0.0, 0.0)

    # Define material parameters for the three biomaterials
    mpar_lignin = MaterialParams(3.5, 0.71, 0.0, 78.0/208.0, 0.208) # (σp1, p₀, c₀, r)
    mpar_cellulose = MaterialParams(3.0, 1.0, 0.0, 45.4/81.0, 0.081) # (σp1, p₀, c₀, r)
    mpar_xylan = MaterialParams(3.0, 1.0, 0.0, 43.0/77.5, 0.0775); # (σp1, p₀, c₀, r)

    # Define a heating profile and a duration for simulation
    start_temp = 300.0;
    end_temp = 800.0;
    heating_duration = 1.0;
    total_duration = 3.0;
    rate = (end_temp - start_temp) / heating_duration;
    Tfun = t -> start_temp + min(rate * t, end_temp);
    max_tstep=heating_duration/10000.

    # Running CPD on the three models
    res_lignin = cpd(AEσb_lignin, AEσg_lignin, AEσρ_lignin, mpar_lignin,
                     total_duration, Tfun, max_tstep = max_tstep, metaplast_model=:original)
    res_cellulose = cpd(AEσb_cellulose, AEσg_cellulose, AEσρ_cellulose, mpar_cellulose,
                        total_duration, Tfun, max_tstep = max_tstep, metaplast_model=:original)
    res_xylan = cpd(AEσb_xylan, AEσg_xylan, AEσρ_xylan, mpar_xylan,
                    total_duration, Tfun, max_tstep = max_tstep, metaplast_model=:original);


    # Print out information about parameters
    display(md"**------- LIGNIN PARAMETERS-------**")
    print_arguments(AEσb_lignin, AEσg_lignin, AEσρ_lignin, mpar_lignin)
    println()
    display(md"**------- CELLULOSE PARAMETERS-------**")
    print_arguments(AEσb_cellulose, AEσg_cellulose, AEσρ_cellulose, mpar_cellulose)
    println()
    display(md"**------- XYLAN PARAMETERS-------**")
    print_arguments(AEσb_xylan, AEσg_xylan, AEσρ_xylan, mpar_xylan)

    println()
    display(Markdown.parse("**duration simulated:** " * string(total_duration) * " seconds"))
    display(md"**Heating from 300K to 800K in 1s.**")
    display(md"**Modified, mass-conservative metaplast model with 20 tar bins.**")
    
    # plotting results
    plot_result(res_lignin, toptitle="Lignin")
    plot_result(res_cellulose, toptitle="Cellulose")
    plot_result(res_xylan, toptitle="Xylan")
    
end

"""
Same as `cpd_biochar_constituents`, but with modified parameters.
"""
function cpd_biochar_consituents_modif()


    # Define reaction rate parameters for the three biomaterials
    AEσb_lignin = ReactionRateParams(7.0e16, 55400.0, 500.0)
    AEσg_lignin = ReactionRateParams(2.3e19, 69000.0, 2600.0)
    AEσρ_lignin = ReactionRateParams(1.7, 0.0, 0.0)

    #AEσb_cellulose = ReactionRateParams(7.0e16, 55400.0, 500.0)
    #AEσg_cellulose = ReactionRateParams(2.3e19, 69000.0, 2600.0)
    #AEσρ_cellulose = ReactionRateParams(1.7, 0.0, 0.0)

    AEσb_cellulose = ReactionRateParams(1.2e20, 51500.0, 100.0)
    AEσg_cellulose = ReactionRateParams(3.0e15, 38200.0, 5000.0)
    AEσρ_cellulose = ReactionRateParams(100, 0.0, 0.0)

    #AEσb_cellulose = ReactionRateParams(2.0e16, 55400.0, 4100.0)
    #AEσg_cellulose = ReactionRateParams(3.0e15, 61200.0, 8100.0)
    #AEσρ_cellulose = ReactionRateParams(100, 0.0, 0.0)

    AEσb_xylan = ReactionRateParams(1.2e20, 51500.0, 100.0)
    AEσg_xylan = ReactionRateParams(3.0e15, 38200.0, 5000.0)
    AEσρ_xylan = ReactionRateParams(100, 0.0, 0.0)

    # Define material parameters for the three biomaterials
    mpar_lignin    = MaterialParams(3.5, 0.71, 0.0, 78.0/208.0, 0.28) # (σp1, p₀, c₀, r, ma)
    mpar_cellulose = MaterialParams(3.5, 0.71, 0.0, 78.0/208.0, 0.28) # (σp1, p₀, c₀, r, ma)
    mpar_xylan     = MaterialParams(3.5, 0.71, 0.0, 78.0/208.0, 0.28) # (σp1, p₀, c₀, r, ma)
    #mpar_cellulose = MaterialParams(3.0, 1.0, 0.0, 45.4/81.0, 0.081) # (σp1, p₀, c₀, r, ma)
    #mpar_xylan = MaterialParams(3.0, 1.0, 0.0, 43.0/77.5, 0.0775); # (σp1, p₀, c₀, r, ma)

    #(sigp1, p0, c0, r, ma) = (4.98, 0.53957, 6.78e-2, 0.35, 1.917e2 / 1000)
   # mpar_lignin = MaterialParams(sigp1, p0, c0, r, ma)
    #mpar_cellulose = MaterialParams(sigp1, p0, c0, r, ma)
    #mpar_xylan = MaterialParams(sigp1, p0, c0, r, ma)

    # Define a heating profile and a duration for simulation
    lfletcher=0;
    if (lfletcher == 0)
        start_temp = 300.0;
        end_temp = 1000.0;
        heating_duration = 1.4;
        total_duration = 1.4;
        rate = (end_temp - start_temp) / heating_duration;
        Tfun = t -> start_temp + min(rate * t, end_temp);
        max_tstep=heating_duration/10000.
    else
        total_duration = 0.03;
        max_tstep=total_duration/10000.
        Tfun = t -> 300 + 1500 * min(1.0, t/1.5e-2)
    end
    
    # Running CPD on the three models
    #res_lignin = cpd(AEσb_lignin, AEσg_lignin, AEσρ_lignin, mpar_lignin,
    #                 total_duration, Tfun, max_tstep = max_tstep, metaplast_model=:original)
    #res_cellulose = cpd(AEσb_cellulose, AEσg_cellulose, AEσρ_cellulose, mpar_cellulose,
    #                    total_duration, Tfun, max_tstep = max_tstep, metaplast_model=:modified)
    res_xylan = cpd(AEσb_xylan, AEσg_xylan, AEσρ_xylan, mpar_xylan,
                    total_duration, Tfun, max_tstep = max_tstep, metaplast_model=:original,
                    consistent_with_original_code=true);
                    #total_duration, Tfun, max_tstep = max_tstep, metaplast_model=:modified);

    #println("time=",res_xylan.time)
    #println("Temperature=",Tfun.(res_xylan.time))


    # Print out information about parameters
    #display(md"**------- LIGNIN PARAMETERS-------**")
    #print_arguments(AEσb_lignin, AEσg_lignin, AEσρ_lignin, mpar_lignin)
    #println()
    #display(md"**------- CELLULOSE PARAMETERS-------**")
    #print_arguments(AEσb_cellulose, AEσg_cellulose, AEσρ_cellulose, mpar_cellulose)
    #println()
    display(md"**------- XYLAN PARAMETERS-------**")
    print_arguments(AEσb_xylan, AEσg_xylan, AEσρ_xylan, mpar_xylan)

    println()
    display(Markdown.parse("**duration simulated:** " * string(total_duration) * " seconds"))
    display(md"**Heating from 300K to 800K in 1s.**")
    #display(md"**Modified, mass-conservative metaplast model with 20 tar bins.**")
    
    # plotting results
    #plot_result(res_lignin, toptitle="Lignin")
    #plot_result(res_cellulose, toptitle="Cellulose")
    plot_result(res_xylan, toptitle="Xylan")
    
end

# ----------------------------------------------------------------------------
"""
Using CPDSpatial to reproduce figure 6a from the paper "Chemical model for
devolatilization. 2. Temperature and Heating Rate Effects on Product Yields"
(1990) by Fletcher et al.  

This example illustrates the effect of the heating rate on tar yield.  
The curves produced display the same qualitative behavior and shape as in the
original paper, but with somewhat different values (the scaling differs a bit). 
Could this be due to slightly different parameters?
"""
function cpd_heating_rate()
    # Define common reaction rate parameters (A, E and σ)
    AEσb = ReactionRateParams(2.6e15, 55400, 1800)
    AEσg = ReactionRateParams(3.0e15, 69000, 8100)
    AEσρ = ReactionRateParams(2.0, 1510, 0); # should equal 0.9 at a temperature of 948 K

    # Define material parameters
    mpar = MaterialParams(4.6, 0.59, 0.11, 0.35); # (σp1, p₀, c₀, r) Illinois No. 6
    
    # Define the temperature profile (note that this will change in the loop below)
    start_temp = 500; # in Kelvin
    end_temp = 1500; # in Kelvin
    rate = 1.0; # in Kelvin per second
    duration = (end_temp - start_temp) / rate;
    Tfun = t -> start_temp + rate * t;

    f = Figure()
    ax1 = Axis(f[1,1],
                       xlabel="Temperature (K)",
                       ylabel="Mass Fraction Tar",
                       title="Tar Yield vs. Temperature",
                       width=500, height=500)
    ax2 = Axis(f[1,2],
                       xlabel="Temperature (K)",
                       ylabel="Mass Fraction Volatiles",
                       title="Volatiles vs. Temperature",
                       width=500, height=500)

    for i = 1:5
        res = cpd(AEσb, AEσg, AEσρ, mpar, duration, Tfun, metaplast_model=:none)
        lines!(ax1, Tfun.(res.time), res.ftar, label="rate = $rate K/s")
        lines!(ax2, Tfun.(res.time), res.ftar .+ res.fgas, label="rate = $rate K/s")
        rate *= 10.0;
        duration = (end_temp - start_temp) / rate;
        Tfun = t -> start_temp + rate * t;
    end
    
    axislegend(ax1, position=:rb)
    axislegend(ax2, position=:rb)
    if !(isCI == "true")
        Mk.activate!(title="Impact of heating rate")
        Mk.resize_to_layout!(f)
        Mk.display(Mk.Screen(), f)
    end

    # Print out information about parameters used to call cpd
    print_arguments(AEσb, AEσg, AEσρ, mpar)
    display(md"**No metaplast model used.**")
    display(md"Note that average site mass is not used since we do not have any 
               metaplast model.  Therefore, it is arbitrarily set to 1.0.")
    println()
    display(Markdown.parse("**duration simulated:** " * string(duration) * " seconds"))
    display(md"**Heating rates vary from 1K/sec to 10.000 K/sec.**")
    
end

# ----------------------------------------------------------------------------
"""
Reproducing CPD simulation results for three coals presented in the paper 
"Chemical model of coal devolatilization using percolation lattice statistics"
(1988) by Grant et al.  These coals are:
- zap lignite
- high volatile bituminous coal
- Montana Rosebud subbituminous coal
Note that there are some discrepancies, which could be partly due to a slightly
different thermal profile.
"""
function cpd_three_coals()

    ## Define common reaction rate parameters (A, E and σ)
    AEσb = ReactionRateParams(2.6e15, 55400, 1800)
    AEσg = ReactionRateParams(3.0e15, 69000, 8100)
    AEσρ = ReactionRateParams(0.9, 0, 0) # constant, equal to 0.9

    # Setting up material parameters
    mpar_zap_lignite = MaterialParams(4.5, 0.61, 0.30, 0.82) # (σp1, p₀, c₀, r)
    mpar_hv_coal = MaterialParams(4.6, 0.59, 0.11, 0.35) # (σp1, p₀, c₀, r)
    mpar_MR_subcoal = MaterialParams(5.8, 0.56, 0.12, 0.35) # (σp1, p₀, c₀, r)

    # duration of simulated time interval
    duration = 70e-3 # 70 milliseconds

    # setting up temperature function
    Tfun = interpolated_temperature_function("data/cpd/sample_temp_profile.txt")

    # calling cpd
    res_zap_lignite = cpd(AEσb, AEσg, AEσρ, mpar_zap_lignite, duration, t -> Tfun(t),
                          metaplast_model=:none)
    res_hv_coal = cpd(AEσb, AEσg, AEσρ, mpar_hv_coal, duration, t -> Tfun(t),
                      metaplast_model=:none)
    res_MR_subcoal = cpd(AEσb, AEσg, AEσρ, mpar_MR_subcoal, duration, t -> Tfun(t),
                         metaplast_model=:none)

    # Print out information about parameters used to call cpd
    print_arguments(AEσb, AEσg, AEσρ)
    display(md"**CPD material parameters, zap lignite**")
    println(mpar_zap_lignite)
    display(md"**CPD material parameters, high volatile bituminous coal**")
    println(mpar_hv_coal)
    display(md"**CPD material parameters, Montana Rosebud subbituminous coal**")
    println(mpar_MR_subcoal)
    println()
    display(Markdown.parse("**duration simulated:** " * string(duration) * " seconds"))
    display(md"**Heating from 300K to 1500K in 15ms.**")
    display(md"**Original metaplast model with 20 tar bins.**")

    # graphical plot of the result    
    plot_result(res_zap_lignite, toptitle="Zap lignite");
    plot_result(res_hv_coal, toptitle="High volatile bitumious coal");
    plot_result(res_MR_subcoal, toptitle="Montana Rosebud subbituminous coal");
end

# ----------------------------------------------------------------------------
"""
Run cpd with the parameters from the original CPDHeat code (Matlab version)
at http://www.et.byu.edu/~tom/cpd/cpdcodes.html
"""
function cpdheat_compare()
    println()
    display(md"
    We here use the parameters from the original CPDheat
    code (http://www.et.byu.edu/~tom/cpd/cpdcodes.html) and compare the results.
    The only parameter we change is the duration, which we shorten so as to
    focus on the interesting parts of the curves.

    The correspondence is relatively good, but there are some minor deviations, e.g. the
    value of the peak of the metaplast curve is slightly higher than in the CPDheat
    code.
    ")
    
    # dividing ma by 1000 as it is expressed in kg/mol rather than g/mol
    (sigp1, p0, c0, r, ma) = (4.98, 0.53957, 6.78e-2, 0.35, 1.917e2 / 1000)

    Tfun_CPDheat = t -> 300 + 1500 * min(1.0, t/1.5e-2)
    duration = 0.03 # NB: ensure the CPDheat model is set to run for the same duration!

    AEσb = ReactionRateParams(2.6e15, 55400, 1800)
    AEσg = ReactionRateParams(3.0e15, 60000, 8100)
    AEσρ = ReactionRateParams(0.9, 0, 0)
    mpar = MaterialParams(sigp1, p0, c0, r, ma)

    res_cpdheat = cpd(AEσb, AEσg, AEσρ, mpar, 
                      duration,
                      Tfun_CPDheat,
                      num_tar_bins=20,
                      max_tstep=duration/1000,
                      metaplast_model=:original)

    # Print out information about parameters used to call cpd
    print_arguments(AEσb, AEσg, AEσρ, mpar)
    display(Markdown.parse("**duration simulated:** " * string(duration) * " seconds"))
    display(md"**Heating from 298K to 1073K in 31ms.**")
    display(md"**Original metaplast model with 20 tar bins.**")

    # graphical plot of the result
    plot_result(res_cpdheat)

end

# ----------------------------------------------------------------------------
"""
Run cpd to reproduce example presented in the original CPD Fortran code at 
http://www.et.byu.edu/~tom/cpd/cpd92/cpdfiles.html.   Output is compared with
the 'cpdb.out' and 'cpd.out' files, which are available at the same website.
"""
function cpd_fortran_compare()

    println()
    display(md"
    In this example, we run a scenario using the parameters available at
    http://www.et.byu.edu/~tom/cpd/cpd92/cpdfiles.html.  We compare the results
    generated by the original CPD Fortran code, available at the same website.

    As can be noted on the generated graphs, the results generated by CPDSpatial
    correspond well with the online published results from the CPD Fortran code.
    ")
    
    # These parameters are taken from the webpage in the docstring abovev
    σp1 = 4.6 # coordination number plus 1
    mdel = 29-7 # mass of side chain (FORTRAN code subtract 7 internally)
    mw1 = 267   # average mass of site, incl. side chains
    mb = 2 * mdel # bridge mass is twice the side chain mass
    ma = mw1 - σp1 * mdel # average mass of site, excl. side chains

    p₀ = 0.61 # initial fraction of intact bridges, including charred bridges
    c₀ = 0.0  # initial fraction of charred bridges
    
    ## Define common reaction rate parameters (A, E and σ)
    AEσb = ReactionRateParams(2.6e15, 55400, 1800);
    AEσg = ReactionRateParams(3.0e15, 69000, 8100);
    AEσρ = ReactionRateParams(0.9, 0, 0); # constant, equal to 0.9
    
    # divide `ma` by 1000 to get it expressed in kg/mol rather than g/mol
    mpar = MaterialParams(σp1, p₀, c₀, mb/ma, ma/1000)
    
    duration = 45e-3 # 45 milliseconds

    # Read the temperature profile from the external file and set up a
    # corresponding temperature function
    Tfun = interpolated_temperature_function("data/cpd/sample_temp_profile.txt")

    # Run the CPD algorithm
    res = cpd(AEσb, AEσg, AEσρ, mpar, duration, t -> Tfun(t),
              num_tar_bins = 20, 
              metaplast_model=:original);

    # Print out information about parameters used to call cpd
    print_arguments(AEσb, AEσg, AEσρ, mpar)
    display(Markdown.parse("**duration simulated:** " * string(duration*1000) *
                           " milliseconds"))
    display(md"**Heating from 298K to 1073K in 31ms.**")
    display(md"**Original metaplast model with 20 tar bins.**")
    
    # Graphical plot of the result
    plot_result(res, toptitle="Result using CPDSpatial")
    
    # Plot the original data as well.  As we can see, the match is very good
    orig_result1 = readdlm("data/cpd/cpd02_cpdb_out.dat", comments=true, comment_char='#')
    orig_result2 = readdlm("data/cpd/cpd02_cpd_out.dat", comments=true, comment_char='#')

    res_fortran = (time = orig_result1[:,1],
                   £vec = orig_result1[:,2],
                   cvec = orig_result1[:,3],
                   δvec = orig_result1[:,4]*2,
                   g1 = orig_result1[:,5]*2,
                   g2 = orig_result1[:,6]*2,
                   fgas = orig_result2[1:end-1,6],
                   ftar = orig_result2[1:end-1,5],
                   fchar = orig_result2[1:end-1, 7],
                   fmetaplast = orig_result2[1:end-1,9])
    plot_result(res_fortran, toptitle="Result from original CPD Fortran code")
    
end

# ----------------------------------------------------------------------------
function interpolated_temperature_function(filename)
    temperature_table = readdlm(filename, ',',
                                comments=true, comment_char='#');

    temperature_table[:,1] *= 1e-3 # convert from milliseconds to seconds
    
    Tfun = linear_interpolation(temperature_table[:,1],
                                temperature_table[:,2],
                                extrapolation_bc=Flat());
end


# ----------------------------------------------------------------------------
function helptext(dispatch, choice)

    println("\nWhen calling this function, you should specify which benchmarking " *
        "case to use.  Choices are:\n")

    for (k, v) in dispatch
        println("----------")
        str1 = "`" * repr(k) * "` - "
        str2 = string(Base.doc(v))
        display(Markdown.parse(str1 * str2))
    end
    println("----------")
    if isnothing(choice)
        println("\nYou did not specify a choice.")
    else
        choice = repr(choice)
        println("\nYou specified the choice `$choice`, which is not in the list above.")
    end
    println()
    println("Note that choice is given as a Julia symbol, so don't forget the initial ':'.")
    
end

# ----------------------------------------------------------------------------
function plot_result(res; toptitle="")

    f = Figure()

    # Create the first subplot with (£vec, δvec, cvec)
    ax1 = Axis(f[1,1],
                       title="Labile Bridges, Side Chains, Char Bridges, gas release",
                       xlabel="Time",
                       ylabel="Fraction",
                       width=500, height=150)
    lines!(ax1, res.time, res.£vec, label="£vec", color = :blue)
    lines!(ax1, res.time, res.δvec/2, label="δvec/2", color = :red)
    lines!(ax1, res.time, res.cvec, label="cvec", color = :green)
    lines!(ax1, res.time, res.g1/2, label="g1/2", color = :coral1)
    lines!(ax1, res.time, res.g2/2, label="g2/2", color = :salmon4)
    axislegend(position=:lc)

    # Create the second subplot with (fgas, ftar, fchar)
    ax2 = Axis(f[2,1],
                       title="Gas Formation, Tar and Char",
                       xlabel="Time",
                       ylabel="Mass Fraction",
                       width=500, height=150)
    lines!(ax2, res.time, res.fgas, label="fgas")
    lines!(ax2, res.time, res.ftar, label="ftar")
    lines!(ax2, res.time, res.fchar, label="fchar")
    if in(:fmetaplast, keys(res))
        lines!(res.time, res.fmetaplast, label="fmetaplast")
    end

    axislegend(position=:lc)
    if !(isCI == "true")
        Mk.activate!(title=toptitle)
        Mk.resize_to_layout!(f)
        Mk.display(Mk.Screen(), f)
    end
end

# ----------------------------------------------------------------------------
function print_arguments(AEσb, AEσg, AEσρ, mpar=nothing)
    println()
    display(md"**Reaction rate parameters for bridge breaking reaction**")
    println(AEσb)
    display(md"**Reaction rate parameters for light gas formation**")
    println(AEσg)
    display(md"**Reaction rate parameters ratio of side chain formation**")
    println(AEσρ)
    if !isnothing(mpar)
        display(md"**CPD material parameters**")
        println(mpar)
    end
end

