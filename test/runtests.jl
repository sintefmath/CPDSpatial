# Copyright (c) 2023 Odd Andersen
# Copyright (c) 2023 SINTEF Digital

using CPDSpatial
using Test
using TestItems
using TestItemRunner

@testitem "CPD" begin
    include("CPDtest.jl")
end

@testitem "spatial" begin
    include("spatial.jl")
end

@run_package_tests
nothing



# @testset "Dummytest" begin
#     @test true
    
# end

# basic CPD tests
# - end values of the fcompare code (original metaplast)
# - end values of biochar (modified metaplast)
# - end values of heating_rate (no metaplast)
# - compare values of lgas, tar, mplast and reattached for spatial_cpd_example
