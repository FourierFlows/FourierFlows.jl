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

@testset "TwoDTurb Stepper tests" begin
    include("test_timesteppers.jl")
end

@testset "BarotropicQG Stepper tests" begin
    include("test_BarotropicQG_timesteppers.jl")
end

@testset "BarotropicQG Rossby wave test" begin
    include("test_BarotropicQG_RossbyWave.jl")
end

@testset "Utils tests" begin
    include("test_utils.jl")
end

println("Total test time: ", toq())
