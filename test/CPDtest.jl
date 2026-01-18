using CPDSpatial
using Interpolations
using DelimitedFiles
using DataStructures


function reproduce_fortran_test()

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

    # temperature function
    Tfun = interpolated_temperature_function("../data/cpd/sample_temp_profile.txt")
    
    # Run the CPD algorithm
    
    res = cpd(AEσb, AEσg, AEσρ, mpar, duration, t -> Tfun(t),
              num_tar_bins = 20,
              consistent_with_original_code = false,
              metaplast_model=:original);

    verify_result(res, [0.32105, 0.16575, 0.28249, 0.55175]);
end

function verify_result(res, target)

    @test isapprox(res.cvec[end], target[1], atol=1e-4)
    @test isapprox(res.fgas[end], target[2], atol=1e-4)
    @test isapprox(res.ftar[end], target[3], atol=1e-4)
    @test isapprox(res.fchar[end], target[4], atol=1e-4)
end

function cpd_biochar_test()

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
    mpar_lignin = MaterialParams(3.5, 0.71, 0.0, 78.0/208.0, 0.28) # (σp1, p₀, c₀, r)
    mpar_cellulose = MaterialParams(3.0, 1.0, 0.0, 45.4/81.0, 0.081) # (σp1, p₀, c₀, r)
    mpar_xylan = MaterialParams(3.0, 1.0, 0.0, 43.0/77.5, 0.0775); # (σp1, p₀, c₀, r)

    # Define a heating profile and a duration for simulation
    start_temp = 300.0;
    end_temp = 800.0;
    heating_duration = 1.0;
    total_duration = 3.0;
    rate = (end_temp - start_temp) / heating_duration;
    Tfun = t -> start_temp + min(rate * t, end_temp);

    # Running CPD on the three models
    res_lignin = cpd(AEσb_lignin, AEσg_lignin, AEσρ_lignin, mpar_lignin,
                     total_duration, Tfun, max_tstep = 1e-2, metaplast_model=:modified,
                     consistent_with_original_code = false)
    res_cellulose = cpd(AEσb_cellulose, AEσg_cellulose, AEσρ_cellulose, mpar_cellulose,
                        total_duration, Tfun, max_tstep = 1e-2, metaplast_model=:modified,
                        consistent_with_original_code = false)
    res_xylan = cpd(AEσb_xylan, AEσg_xylan, AEσρ_xylan, mpar_xylan,
                    total_duration, Tfun, max_tstep = 1e-2, metaplast_model=:modified,
                    consistent_with_original_code = false)

    verify_result(res_lignin, [0.26296, 0.29029, 0.243073, 0.46663])
    verify_result(res_cellulose, [0.0099, 0.175608, 0.6645, 0.1599])
    verify_result(res_xylan, [0.0099, 0.2673, 0.55061, 0.1821])
end


function interpolated_temperature_function(filename)
    temperature_table = readdlm(filename, ',',
                                comments=true, comment_char='#');

    temperature_table[:,1] *= 1e-3 # convert from milliseconds to seconds
    
    Tfun = linear_interpolation(temperature_table[:,1],
                                temperature_table[:,2],
                                extrapolation_bc=Flat());
end



@testset "reproduce fortran test" begin
    reproduce_fortran_test()
end

@testset "biochar material test" begin
    cpd_biochar_test()
end


# @test cpd_biochar_test()
# @test heating_rate_test()

