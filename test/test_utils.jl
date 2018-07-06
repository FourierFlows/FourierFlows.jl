
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
Compute the J(a, b) and compare with analytic_answer.
"""
function test_jacobian(a, b, analytic_answer, grid)
    # it's important to use atol for this test since when analytic_answer=0
    # the rtol is not conclusive (e.g., isapprox(1e-5, 0, rtol=1e-10) is false)
    isapprox(FourierFlows.jacobian(a, b, grid), analytic_answer;
             atol=g.nx*g.ny*1e-14)
end

"""
Test the createarrays macro.
"""
function test_createarrays(T=Float64, dims=(13, 45))
  a1, b1 = zeros(T, dims), zeros(T, dims)
  FourierFlows.@createarrays T dims a2 b2
  a1 == a2 && b1 == b2
end

"""
Test the peakedisotropicspectrum function.
"""
function testpeakedisotropicspectrum()
  n, L = 128, 2π
  gr = TwoDGrid(n, L)
  k0, E0 = 6, 0.5
  qi = FourierFlows.peakedisotropicspectrum(gr, k0, E0; allones=true)
  ρ, qhρ = FourierFlows.radialspectrum(rfft(qi).*gr.invKKrsq, gr)

  ρtest = ρ[ (ρ.>15.0) .& (ρ.<=17.5)]
  qhρtest = qhρ[ (ρ.>15.0) .& (ρ.<=17.5)]

  isapprox(abs.(qhρtest)/abs(qhρtest[1]), (ρtest/ρtest[1]).^(-2), rtol=5e-3)
end



abstract type TestVars <: AbstractVars end

physicalvars = [:a, :b]
transformvars = [ Symbol(var, :h) for var in physicalvars ]
 forcedvar = [:Fh]

varspecs = cat(1,
  FourierFlows.getfieldspecs(physicalvars, :(Array{T,2})),
  FourierFlows.getfieldspecs(transformvars, :(Array{Complex{T},2})))

eval(FourierFlows.structvarsexpr(:Vars, varspecs; parent=:TestVars))

function Vars(g)
  @createarrays typeof(g.Lx) (g.nx, g.ny) a b
  @createarrays Complex{typeof(g.Lx)} (g.nkr, g.nl) ah bh
  Vars(a, b, ah, bh)
end

function teststructvarsexpr(g)
  v = Vars(g)
  (
    typeof(v.a)==Array{Float64,2} &&
    typeof(v.ah)==Array{Complex{Float64},2} &&
    size(v.a)==size(v.b) &&
    size(v.ah)==size(v.bh) &&
    size(v.a)==(g.nx, g.ny) &&
    size(v.ah) == (g.nkr, g.nl)
  )
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

# Sine/Exp waves
k1, l1 = 2*k0, 6*l0
k2, l2 = 3*k0, -3*l0

sin1 = sin.(k1*x + l1*y)
sin2 = sin.(k2*x + l2*y)
exp1 = exp.(im*(k1*x + l1*y))
exp2 = exp.(im*(k2*x + l2*y))

# Analytical expression for the Jacobian of sin1 and sin2 and of exp1 and exp2
Jsin1sin2 = (k1*l2-k2*l1)*cos.(k1*x + l1*y).*cos.(k2*x + l2*y)
Jexp1exp2 = (k2*l1-k1*l2)*exp.(im*((k1+k2)*x + (l1+l2)*y))

@test test_parsevalsum(f1, g; realvalued=true)   # Real valued f with rfft
@test test_parsevalsum(f1, g; realvalued=false)  # Real valued f with fft
@test test_parsevalsum(f2, g; realvalued=false)  # Complex valued f with fft
@test test_parsevalsum2(f1, g; realvalued=true)  # Real valued f with rfft
@test test_parsevalsum2(f1, g; realvalued=false) # Real valued f with fft
@test test_parsevalsum2(f2, g; realvalued=false) # Complex valued f with fft

@test test_jacobian(sin1, sin1, 0*sin1, g)  # Test J(a, a) = 0
@test test_jacobian(sin1, sin2, Jsin1sin2, g) # Test J(sin1, sin2) = Jsin1sin2
@test test_jacobian(exp1, exp2, Jexp1exp2, g) # Test J(exp1, exp2) = Jexp1exps2

@test test_createarrays()
@test test_fftwavenums()
@test test_rms(32)
@test test_domainaverage(32; xdir=true)
@test test_domainaverage(32; xdir=false)
@test teststructvarsexpr(g)

# Radial spectrum tests. Note that ahρ = ∫ ah ρ dθ.
n = 128; δ = n/10                 # Parameters
ahkl(k, l) = exp(-(k^2+l^2)/2δ^2) #  a = exp(-ρ²/2δ²)
    ahρ(ρ) = 2π*ρ*exp(-ρ^2/2δ^2)  # aᵣ = 2π ρ exp(-ρ²/2δ²)
@test test_radialspectrum(n, ahkl, ahρ)

ahkl(k, l) = exp(-(k^2+l^2)/2δ^2) * k^2/(k^2+l^2) #  a = exp(-ρ²/2δ²)*cos(θ)²
    ahρ(ρ) = π*ρ*exp(-ρ^2/2δ^2)                   # aᵣ = π ρ exp(-ρ²/2δ²)
@test test_radialspectrum(n, ahkl, ahρ)

@test testpeakedisotropicspectrum()
