include("../src/fourierflows.jl")

using Base.Test
using FourierFlows
import FourierFlows.TwoDTurb


function integralsquare(func, grid)
    sum(abs2.(func))*grid.dx*grid.dy
end

function testparsevalsum(func, grid; realvalued=true)
    # compute the integral in physical space
    integral = integralsquare(func, grid)

    # compute the integral in wavenumber space using Parseval's theorem
    if realvalued==true
        funch = rfft(func)
    elseif realvalued==false
        funch = fft(func)
    end
    parsevalsum = FourierFlows.parsevalsum(abs2.(funch), grid)

    # compare the two results
    isapprox(integral, parsevalsum; atol=1.0e-14)
end

function testparsevalsum2(func, grid; realvalued=true)
    # compute the integral in physical space
    integral = integralsquare(func, grid)

    # compute the integral in wavenumber space using Parseval's theorem
    if realvalued==true
        funch = rfft(func)
    elseif realvalued==false
        funch = fft(func)
    end
    parsevalsum2 = FourierFlows.parsevalsum2(funch, grid)

    # compare the two results
    isapprox(integral, parsevalsum2; atol=1.0e-14)
end


# Test on a rectangular grid
nx, ny = 64, 128             # number of points
Lx, Ly = 2π, 3π             # Domain width
g = TwoDGrid(nx, Lx, ny, Ly)

x, y = g.X, g.Y
k0, l0 = 2π/Lx, 2π/Ly
sig = 0.5

# a real-valued gaussian
f1 = exp.(-(x.^2 + y.^2)/(2*sig^2))

# a complex-valued function
f2 = exp.( im*(2k0*x + 3l0*y.^2) ).*( exp.(-(x.^2 + y.^2)/(2sig^2))
                                        + 2im*exp.(-(x.^2 + y.^2)/(5sig^2)) )

@testset "Parsevalsum Tests" begin
  @test testparsevalsum(f1, g; realvalued=true)   #real valued f with rfft
  @test testparsevalsum(f1, g; realvalued=false)  #real valued f with fft
  @test testparsevalsum(f2, g; realvalued=false)  #complex valued f with fft
  @test testparsevalsum2(f1, g; realvalued=true)  #real valued f with rfft
  @test testparsevalsum2(f1, g; realvalued=false) #real valued f with fft
  @test testparsevalsum2(f2, g; realvalued=false) #complex valued f with fft
end
