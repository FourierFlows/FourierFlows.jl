#!/usr/bin/env julia

#Start Test Script
using FourierFlows
using Base.Test

# Run tests

tic()
println("Test 1")
@time @test include("test_grid.jl")
toc()
