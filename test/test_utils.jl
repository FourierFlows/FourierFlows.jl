function test_domainaverage(n; xdir=true)
  g = TwoDGrid(n, 2π)
  X, Y = gridpoints(g)
  if xdir; c = @. cos(X)^2
  else;    c = @. cos(Y)^2
  end
  0.5 ≈ FourierFlows.domainaverage(c, g)
end

# This test could use some further work.
function test_radialspectrum(n, ahkl, ahρ; debug=false, atol=0.1)
  g = TwoDGrid(n, 2π)
  ah = @. ahkl(g.k, g.l)
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

abstract type TestVars <: AbstractVars end

physicalvars = [:a, :b]
transformvars = [ Symbol(var, :h) for var in physicalvars ]
 forcedvar = [:Fh]

varspecs = cat(FourierFlows.getfieldspecs(physicalvars, :(Array{T,2})),
  FourierFlows.getfieldspecs(transformvars, :(Array{Complex{T},2})), dims=1)

eval(FourierFlows.structvarsexpr(:VarsFields, physicalvars, transformvars))
eval(FourierFlows.structvarsexpr(:VarsFieldsParent, physicalvars, transformvars; parent=:TestVars))
eval(FourierFlows.structvarsexpr(:VarsSpecs, varspecs))
eval(FourierFlows.structvarsexpr(:VarsSpecsParent, varspecs; parent=:TestVars))

function VarsFields(g)
  @createarrays typeof(g.Lx) (g.nx, g.ny) a b
  @createarrays Complex{typeof(g.Lx)} (g.nkr, g.nl) ah bh
  VarsFields(a, b, ah, bh)
end

function VarsFieldsParent(g)
  @createarrays typeof(g.Lx) (g.nx, g.ny) a b
  @createarrays Complex{typeof(g.Lx)} (g.nkr, g.nl) ah bh
  VarsFieldsParent(a, b, ah, bh)
end

function VarsSpecs(g)
  @createarrays typeof(g.Lx) (g.nx, g.ny) a b
  @createarrays Complex{typeof(g.Lx)} (g.nkr, g.nl) ah bh
  VarsSpecs(a, b, ah, bh)
end

function VarsSpecsParent(g)
  @createarrays typeof(g.Lx) (g.nx, g.ny) a b
  @createarrays Complex{typeof(g.Lx)} (g.nkr, g.nl) ah bh
  VarsSpecsParent(a, b, ah, bh)
end

function test_structvarsexprFields(g)
  v1 = VarsFields(g)
  (
    typeof(v1.a)==Array{Float64,2} &&
    typeof(v1.ah)==Array{Complex{Float64},2} &&
    size(v1.a)==size(v1.b) &&
    size(v1.ah)==size(v1.bh) &&
    size(v1.a)==(g.nx, g.ny) &&
    size(v1.ah) == (g.nkr, g.nl)
  )
  v2 = VarsFieldsParent(g)
  (
    typeof(v2.a)==Array{Float64,2} &&
    typeof(v2.ah)==Array{Complex{Float64},2} &&
    size(v2.a)==size(v2.b) &&
    size(v2.ah)==size(v2.bh) &&
    size(v2.a)==(g.nx, g.ny) &&
    size(v2.ah) == (g.nkr, g.nl)
  )
end

function test_structvarsexprSpecs(g)
  v1 = VarsSpecs(g)
  (
    typeof(v1.a)==Array{Float64,2} &&
    typeof(v1.ah)==Array{Complex{Float64},2} &&
    size(v1.a)==size(v1.b) &&
    size(v1.ah)==size(v1.bh) &&
    size(v1.a)==(g.nx, g.ny) &&
    size(v1.ah) == (g.nkr, g.nl)
  )
  v2 = VarsSpecsParent(g)
  (
    typeof(v2.a)==Array{Float64,2} &&
    typeof(v2.ah)==Array{Complex{Float64},2} &&
    size(v2.a)==size(v2.b) &&
    size(v2.ah)==size(v2.bh) &&
    size(v2.a)==(g.nx, g.ny) &&
    size(v2.ah) == (g.nkr, g.nl)
  )
end
