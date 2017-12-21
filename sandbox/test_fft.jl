include("../src/FourierFlows.jl")

using FourierFlows, Base.Test

module FFTTests

""" Test how the real forward transform affect input. """
function test_input_rfft(grid)
  A = rand(grid.nx, grid.ny)
  B = zeros(Complex{Float64}, grid.nkr, grid.nl)

  Acopy = deepcopy(A)
  A_mul_B!(B, grid.rfftplan, A)

  A == Acopy
end


""" Test how the real forward transform affect input. """
function test_input_irfft(grid)
  A = rand(grid.nx, grid.ny)
  B = zeros(Complex{Float64}, grid.nkr, grid.nl)
  Acopy = deepcopy(A)

  A_mul_B!(B, grid.rfftplan, Acopy)
  Bcopy = deepcopy(B)

  A_mul_B!(A, grid.irfftplan, B)

  B == Bcopy
end


""" Test how the real forward transform affect input. """
function test_copying_irfft(grid)
  A = rand(grid.nx, grid.ny)
  B1 = zeros(Complex{Float64}, grid.nkr, grid.nl)
  B2 = zeros(Complex{Float64}, grid.nkr, grid.nl)
  Acopy = deepcopy(A)

  A_mul_B!(B1, grid.rfftplan, Acopy)
  Bcopy = deepcopy(B1)
  @. B2 = B1

  A_mul_B!(A, grid.irfftplan, B1)

  B2 == Bcopy
end





""" Test an rfft, irfft cycle. """
function test_cycle_rfft(grid)

  tol = 1e-15

  A = rand(grid.nx, grid.ny)
  B = zeros(Complex{Float64}, grid.nkr, grid.nl)
  Aorig = deepcopy(A)
  Anew  = zeros(grid.X)

  A_mul_B!(B, grid.rfftplan, A)
  A_mul_B!(Anew, grid.irfftplan, B)

  maximum(abs.(Anew-Aorig)) < tol
end


end # end FFTTests


import FFTTests: test_input_rfft, test_input_irfft, test_cycle_rfft,
  test_copying_irfft


# Run tests
grid = TwoDGrid(64, 2π)

@test test_input_rfft(grid)

println("test_input_rfft:    ", @test test_input_rfft(grid))
println("test_input_irfft (does irfft leave its input intact?):   ", test_input_irfft(grid))
println("test_cycle_rfft:    ", @test test_cycle_rfft(grid))
println("test_copying_irfft: ", @test test_copying_irfft(grid))
