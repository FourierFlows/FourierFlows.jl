# import FourierFlows.TwoDTurb

test_fftwavenums() = FourierFlows.fftwavenums(6; L=2π) == [0, 1, 2, 3, -2, -1]  

function test_domainaverage(n; xdir=true)
  g = TwoDGrid(n, 2π)
  if xdir; c = cos.(g.X).^2
  else;    c = cos.(g.Y).^2
  end
  0.5 ≈ FourierFlows.domainaverage(c, g)
end

function test_rms(n)
  g = TwoDGrid(n, 2π)
  q = cos.(g.X)
  isapprox(FourierFlows.rms(q), sqrt(1/2))
end

# This test could use some further work.
function test_radialspectrum(n, ahkl, ahρ; debug=false, atol=0.1)
  g = TwoDGrid(n, 2π)
  ah = ahkl.(g.K, g.L)
  ah[1, 1] = 0.0

  ρ, ahρ_estimate = FourierFlows.radialspectrum(ah, g; refinement=16)

  if debug
    println(sum(ahρ.(ρ) - ahρ_estimate))
    return ah, ahρ_estimate, ahρ.(ρ)
  else
    normalizeddiff = sum(abs.(ahρ.(ρ) - ahρ_estimate)) / length(ρ)
    return isapprox(normalizeddiff, 0.0, atol=atol)
  end
end

function integralsquare(func, grid)
    sum(abs2.(func))*grid.dx*grid.dy
end

function test_parsevalsum(func, grid; realvalued=true)
  # Compute integral in physical space
  integral = integralsquare(func, grid)
  # Compute integral in wavenumber space using Parseval's theorem
  if realvalued==true;      funch = rfft(func)
  elseif realvalued==false; funch = fft(func)
  end
  parsevalsum = FourierFlows.parsevalsum(abs2.(funch), grid)
  isapprox(integral, parsevalsum; atol=1.0e-14)
end

function test_parsevalsum2(func, grid; realvalued=true)
  # Compute integral in physical space
  integral = integralsquare(func, grid)
  # Compute integral in wavenumber space using Parseval's theorem
  if realvalued==true;      funch = rfft(func)
  elseif realvalued==false; funch = fft(func)
  end
  parsevalsum2 = FourierFlows.parsevalsum2(funch, grid)
  isapprox(integral, parsevalsum2; atol=1.0e-14)
end


"""
Compute the J(a,b) and compare with analytic_answer.
"""
function test_jacobian(a, b, analytic_answer, grid)
    isapprox(norm(FourierFlows.jacobian(a, b, grid)), norm(analytic_answer); 
             atol=g.nx*g.ny*1e-14)
end

function test_createarrays(T=Float64, dims=(13, 45))
  a1, b1 = zeros(T, dims), zeros(T, dims)
  FourierFlows.@createarrays T dims a2 b2
  a1 == a2 && b1 == b2
end


# Run tests -------------------------------------------------------------------
# Test on a rectangular grid
nx, ny = 64, 128   # number of points
Lx, Ly = 2π, 3π    # Domain width
g = TwoDGrid(nx, Lx, ny, Ly)
x, y = g.X, g.Y
k0, l0 = 2π/Lx, 2π/Ly

# Real and complex-valued functions
σ = 0.5
f1 = exp.(-(x.^2 + y.^2)/(2σ^2))
f2 = exp.( im*(2k0*x + 3l0*y.^2) ).*( 
      exp.(-(x.^2 + y.^2)/(2σ^2)) + 2im*exp.(-(x.^2 + y.^2)/(5σ^2)) )
                                        
# Sine waves
k1, l1 = 2*k0, 6*l0
k2, l2 = 3*k0, -3*l0
s1 = sin.(k1*x + l1*y)
s2 = sin.(k2*x + l2*y)

# Analytical expression for the Jacobian of s1 and s2
Js1s2 = (k1*l2-k2*l1)*cos.(k1*x + l1*y).*cos.(k2*x + l2*y)

@test test_parsevalsum(f1, g; realvalued=true)   # Real valued f with rfft
@test test_parsevalsum(f1, g; realvalued=false)  # Real valued f with fft
@test test_parsevalsum(f2, g; realvalued=false)  # Complex valued f with fft
@test test_parsevalsum2(f1, g; realvalued=true)  # Real valued f with rfft
@test test_parsevalsum2(f1, g; realvalued=false) # Real valued f with fft
@test test_parsevalsum2(f2, g; realvalued=false) # Complex valued f with fft

@test test_jacobian(s1, s1, 0*s1, g)  # Test J(a,a) = 0
@test test_jacobian(s1, s2, Js1s2, g) # Test J(s1,s2) = Js1s2

@test test_createarrays()
@test test_fftwavenums()
@test test_rms(32)
@test test_domainaverage(32; xdir=true)
@test test_domainaverage(32; xdir=false)

# Radial spectrum tests. Note that ahρ = ∫ ah ρ dθ.
n = 128; δ = n/10                 # Parameters
ahkl(k, l) = exp(-(k^2+l^2)/2δ^2) #  a = exp(-ρ²/2δ²)
    ahρ(ρ) = 2π*ρ*exp(-ρ^2/2δ^2)  # aᵣ = 2π ρ exp(-ρ²/2δ²)
@test test_radialspectrum(n, ahkl, ahρ)

ahkl(k, l) = exp(-(k^2+l^2)/2δ^2) * k^2/(k^2+l^2) #  a = exp(-ρ²/2δ²)*cos(θ)²
    ahρ(ρ) = π*ρ*exp(-ρ^2/2δ^2)                   # aᵣ = π ρ exp(-ρ²/2δ²)
@test test_radialspectrum(n, ahkl, ahρ)
