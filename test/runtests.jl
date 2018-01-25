#!/usr/bin/env julia

# Start Test Script
using FourierFlows
using Base.Test

# Run tests
tic()

println("-- Core tests --")

@testset "Grid tests" begin
    include("test_grid.jl")
end

@testset "FFT tests" begin
    include("test_fft.jl")
end

@testset "IFFT tests" begin
    include("test_ifft.jl")
end

@testset "Utils tests" begin
    include("test_utils.jl")
end

@testset "Timestepper tests" begin
    include("test_timesteppers.jl")
end

# @testset "BarotropicQG Stepper tests" begin
#     include("test_BarotropicQG_timesteppers.jl")
# end

println("-- Physics tests --")

@testset "Physics: TwoDTurb" begin
  include("test_twodturb.jl")
end

@testset "Physics: BarotropicQG" begin
    include("test_barotropicqg.jl")
end


println("Total test time: ", toq())
