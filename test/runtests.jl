#!/usr/bin/env julia

# using CuArrays
using FourierFlows
using Requires
using Test

# Run tests

time = @elapsed begin

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

@testset "Diagnostics tests" begin
  include("test_diagnostics.jl")
end


println("-- Physics tests --")

@testset "Physics: Kuramoto-Sivashinsky" begin
  include("test_kuramotosivashinsky.jl")
end

@testset "Physics: TracerAdvDiff" begin
  include("test_traceradvdiff.jl")
end

@require CuArrays="3a865a2d-5b23-5a0f-bc46-62713ec82fae" begin

  println("-- CUDA tests --")

  @testset "CuGrid tests" begin
    using CuArrays
    include("test_cugrid.jl")
  end

end

end
println("Total test time: ", time)
