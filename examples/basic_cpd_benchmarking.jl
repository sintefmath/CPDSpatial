using Markdown
using Plots
using Plots: plot 
using PyPlot # for pyplot backend to Plots
using Interpolations
using DelimitedFiles

# ----------------------------------------------------------------------------
function cpd_benchmarking(choice=nothing)

    dispatch = Dict(:cpdheat => cpdheat_compare,
                    :fcompare => cpd_fortran_compare)

    get(dispatch, choice, () -> helptext(dispatch, choice))()

    
end

# ----------------------------------------------------------------------------
"""
Run cpd with the parameters from the original CPDHeat code (Matlab version)
at http://www.et.byu.edu/~tom/cpd/cpdcodes.html
"""
function cpdheat_compare()
    # CPDheat compare:
    # We here use the parameters from the original CPDheat code and compare the results.
    # The only parameter we change is the duration, which we shorten so as to focus
    # on the interesting parts of the curves.
    # The correspondence is quite good, but there are some minor deviations, e.g. the
    # value of the peak of the metaplast curve is slightly higher than in the CPDheat
    # code.
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
    display(md"**Heating from 300K to 1500K in 15ms.**")
    display(md"**Original metaplast model with 20 tar bins.**")

    # graphical plot of the result
    plot_result(res_cpdheat)

end

# ----------------------------------------------------------------------------
function print_arguments(AEσb, AEσg, AEσρ, mpar)
    println()
    display(md"**Reaction rate parameters for bridge breaking reaction**")
    println(AEσb)
    display(md"**Reaction rate parameters for light gas formation**")
    println(AEσg)
    display(md"**Reaction rate parameters ratio of side chain formation**")
    println(AEσρ)
    display(md"**CPD material parameters**")
    println(mpar)
end

# ----------------------------------------------------------------------------
"""
Run cpd to reproduce example presented in the original CPD fortran code at 
http://www.et.byu.edu/~tom/cpd/cpd92/cpdfiles.html.   Output is compared with
the 'cpdb.out' file, which is available on that website.
"""
function cpd_fortran_compare()
    println("Inside cpd_fortran_compare")
    error("UNIMPLEMENTED")
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
    
    if isnothing(choice)
        println("\nYou did not specify a choice.")
    else
        choice = repr(choice)
        println("\nYou specified the choice `$choice`, which is not in the list above.")
    end
end

# ----------------------------------------------------------------------------
function plot_result(res)

    # Create the first subplot with (£vec, δvec, cvec)
    subplot1 = plot(res.time, res.£vec, label="£vec", reuse=false)
    plot!(res.time, res.δvec/2, label="δvec")
    plot!(res.time, res.cvec, label="cvec")
    plot!(res.time, res.g1/2, label="g1")
    plot!(res.time, res.g2/2, label="g2")
    xlabel!("Time")
    ylabel!("Fraction")
    title!("Labile Bridges, Side Chains, Char Bridges, gas release")

    # Create the second subplot with (fgas, ftar, fchar)
    subplot2 = plot(res.time, res.fgas, label="fgas", reuse=false)
    plot!(res.time, res.ftar, label="ftar")
    plot!(res.time, res.fchar, label="fchar")

    if in(:fmetaplast, keys(res))
        plot!(res.time, res.fmetaplast, label="fmetaplast")
    end
    
    xlabel!("Time")
    ylabel!("Mass Fraction")
    title!("Gas Formation, Tar, and Char")

    @layout [a b]
    p = plot(subplot1, subplot2, legend=:bottomright, size=(1200, 500))
    display(p)
end
