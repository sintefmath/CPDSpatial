using Plots
using Plots: plot 
using PyPlot # for pyplot backend to Plots
using Interpolations
using DelimitedFiles

pyplot() # set pyplot backend for Plots

# ============================================================================
# Define temperature profile encountered in some published work

temperature_file = "data/cpd/sample_temp_profile.txt";

# Read the matrix from the external file
temperature_table = readdlm(temperature_file, ',',
                            comments=true, comment_char='#');
temperature_table[:,1] *= 1e-3 # convert from milliseconds to seconds

Tfun = linear_interpolation(temperature_table[:,1],
                            temperature_table[:,2],
                            extrapolation_bc=Flat());

## Define common reaction rate parameters (A, E and σ)
AEσb = ReactionRateParams(2.6e15, 55400, 1800);
AEσg = ReactionRateParams(3.0e15, 69000, 8100);
AEσρ = ReactionRateParams(0.9, 0, 0); # constant, equal to 0.9


# ============================================================================
# Define plotting function

function plot_result(res)

    # Create the first subplot with (£vec, δvec, cvec)
    subplot1 = plot(res.time, res.£vec, label="£vec", reuse=false)
    plot!(res.time, res.δvec, label="δvec")
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


# ============================================================================
# CPDheat compare:
# We here use the parameters from the original CPDheat code and compare the results.
# The only parameter we change is the duration, which we shorten so as to focus
# on the interesting parts of the curves.
# The correspondence is quite good, but note that there are some deviations, e.g. the
# value of the peak of the metaplast curve is slightly higher than in the CPDheat
# code.

# dividing ma by 1000 as it is expressed in kg/mol rather than g/mol
(sigp1, p0, c0, r, ma) = (4.98, 0.53957, 6.78e-2, 0.35, 1.917e2 / 1000);
Tfun_CPDheat = t -> 300 + 1500 * min(1.0, t/1.5e-2);
duration = 0.03 # NB: ensure the CPDheat model is set to run for the same duration!

res_cpdheat = cpd(ReactionRateParams(2.6e15, 55400, 1800),
                  ReactionRateParams(3.0e15, 60000, 8100),
                  ReactionRateParams(0.9, 0, 0),
                  MaterialParams(sigp1, p0, c0, r, ma),
                  duration,
                  Tfun_CPDheat,
                  num_tar_bins=20,
                  max_tstep=duration/1000,
                  basic_model=false); 

plot_result(res_cpdheat)

# ============================================================================
# Reproducing the example presented with the original CPD FORTRAN code at
# http://www.et.byu.edu/~tom/cpd/cpd92/cpdfiles.html
# Output is compared with the 'cpdb.out' file available at that web page

mpar = MaterialParams(4.6, 0.61, 0.0, (29.0-7.0)/267.0, 267.0)
duration = 45e-3 # 45 milliseconds
res1_fortran = cpd(AEσb, AEσg, AEσρ, mpar, duration, t -> Tfun(t), basic_model=false);

plot_result(res1_fortran)

# Plot the original data as well.  As we can see, the match is very good
orig_result1 = readdlm("data/cpd/cpd02_cpdb_out.dat", comments=true, comment_char='#')
plot(orig_result1[:,1], orig_result1[:,2], label="£vec", reuse=false)
plot!(orig_result1[:,1], orig_result1[:,4], label="δ", reuse=true)
plot!(orig_result1[:,1], orig_result1[:,3], label="c", reuse=true)
plot!(orig_result1[:,1], orig_result1[:,5], label="g1 ", reuse=true)
plot!(orig_result1[:,1], orig_result1[:,6], label="g2", reuse=true)

orig_result2 = readdlm("data/cpd/cpd02_cpd_out.dat", comments=true, comment_char='#')
plot(orig_result2[:,1], orig_result2[:,5], label="ftar", reuse=false)
plot!(orig_result2[:,1], orig_result2[:,3], label="fcross", reuse=true)
plot!(orig_result2[:,1], orig_result2[:,7], label="fchar ", reuse=true)
plot!(orig_result2[:,1], orig_result2[:,6], label="fgas", reuse=true)
plot!(orig_result2[:,1], orig_result2[:,9], label="fmet", reuse=true)

# ============================================================================
# Reproducing 'zap lignite' example from paper "Chemical model of coal
# devolatilization using percolation lattice statistics" (1988) by Grant et al.
# Note that there are discrepancies, which may be partly due to a slightly
# different thermal profile(?)

mpar = MaterialParams(4.5, 0.61, 0.30, 0.82); # (σp1, p₀, c₀, r)
duration = 70e-3; # 70 milliseconds

res_zap_lignite = cpd(AEσb, AEσg, AEσρ, mpar, duration, t -> Tfun(t));

plot_result(res_zap_lignite);

# ============================================================================
# Reproducing 'high volatile bituminous coal' example from the same paper
# Note that there are discrepancies, which may be partly due to a slightly
# different thermal profile(?)

mpar = MaterialParams(4.6, 0.59, 0.11, 0.35); # (σp1, p₀, c₀, r)
duration = 70e-3; # 70 milliseconds

res_hvbc = cpd(AEσb, AEσg, AEσρ, mpar, duration, t -> Tfun(t));

plot_result(res_hvbc);

# ============================================================================
# Reproducing 'Montana Rosebud subbituminous coal' example from the same paper
# Note that there are discrepancies, which are likely due to the thermal profile
# used being slighly different.

mpar = MaterialParams(5.8, 0.56, 0.12, 0.35); # (σp1, p₀, c₀, r)
duration = 70e-3; # 70 milliseconds

res_rosebud = cpd(AEσb, AEσg, AEσρ, mpar, duration, t -> Tfun(t));

plot_result(res_rosebud);

# ============================================================================
# Reproducing figure 6a from the paper "Chemical model for
# devolatilization. 2. Temperature and Heating Rate Effects on Product Yields"
# (1990) by Fletcher et al.

AEσρ2 = ReactionRateParams(2.0, 1510, 0); # should equal 0.9 at a temperature of 948 K
mpar = MaterialParams(4.6, 0.59, 0.11, 0.35); # (σp1, p₀, c₀, r) Illinois No. 6

start_temp = 500; # in Kelvin
end_temp = 1500; # in Kelvin

pyplot()

p1 = plot(xlims=(start_temp, end_temp), ylims=(0, 0.3), xlabel="Temperature (K)",
         ylabel="Mass Fraction Tar", title="Tar Yield vs. Temperature")

p2 = plot(xlims=(start_temp, end_temp), ylims=(0, 0.8), xlabel="Temperature (K)",
         ylabel="Mass Fraction Volatiles", title="Volatiles vs. Temperature")

rate = 1.0; # in Kelvin per second
duration = (end_temp - start_temp) / rate;
Tfun = t -> start_temp + rate * t;

for i = 1:5
    res = cpd(AEσb, AEσg, AEσρ2, mpar, duration, Tfun, basic_model=true);
    plot!(p1, Tfun.(res.time), res.ftar, label="rate = $rate K/s", reuse=false)
    plot!(p2, Tfun.(res.time), res.ftar .+ res.fgas, label="rate = $rate K/s", reuse=false)
    rate *= 10.0;
    duration = (end_temp - start_temp) / rate;
    Tfun = t -> start_temp + rate * t;
end
@layout [a b]
p = plot(p1, p2, legend=:bottomright, size=(2000, 500))
display(p)

# Note: the curves display the same qualitative behavior and shape as in the
# paper, but with somewhat different values (the scaling differs a bit).  Could
# this be due to slightly different parameters?
