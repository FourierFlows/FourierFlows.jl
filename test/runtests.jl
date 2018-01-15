#!/usr/bin/env julia

# Start Test Script
using FourierFlows
using Base.Test

# Run tests
tic()

@testset "Grid tests" begin
    include("test_grid.jl")
end

@testset "FFT tests" begin
    include("test_fft.jl")
end

@testset "IFFT tests" begin
    include("test_ifft.jl")
end

@testset "Stepper tests" begin
    include("test_timesteppers.jl")
end

#@time @testset "BarotropicQG and Timestepper tests" begin
#    include("test_BarotropicQG_timestep.jl")
#end

@testset "Utils tests" begin
    include("test_utils.jl")
end

println("Total test time: ", toq())
