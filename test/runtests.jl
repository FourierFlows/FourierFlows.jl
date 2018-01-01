#!/usr/bin/env julia

#Start Test Script
using FourierFlows
using Base.Test

# Run tests

tic()
println("Starting Tests")
println(" ")

@time @testset "Grid tests" begin
    include("test_grid.jl")
end
@time @testset "FFT tests" begin
    include("test_fft.jl")
end
@time @testset "IFFT tests" begin
    include("test_ifft.jl")
end
@time @testset "TwoDTurb and Timestepper tests" begin
    include("test_twodturb.jl")
end
@time @testset "Utils tests" begin
    include("test_utils.jl")
end

toc()
