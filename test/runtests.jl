#!/usr/bin/env julia

#using CuArrays
using Requires
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


println("-- Physics tests --")

@testset "Physics: TwoDTurb" begin
  include("test_twodturb.jl")
end

@testset "Physics: BarotropicQG" begin
  include("test_barotropicqg.jl")
end

@testset "Physics: Kuramoto-Sivashinsky" begin
  include("test_kuramotosivashinsky.jl")
end

@require CuArrays begin

  println("-- CUDA tests --")

  @testset "CuGrid tests" begin
    using CuArrays
    include("test_cugrid.jl")
  end

end

println("Total test time: ", toq())
