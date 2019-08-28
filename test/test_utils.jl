function test_fltype()
  Tf = Float16
  Tc = Complex{Tf}
  Tf == fltype(Tf) && Tf == fltype(Tc)
end

function test_cxtype()
  Tf = Float16
  Tc = Complex{Tf}
  Tc == cxtype(Tf) && Tc == cxtype(Tc)
end

function test_innereltype(T=Float32)
  a = [zeros(T, 3), zeros(T, 2, 5)]
  T == innereltype(a)
end

function test_superzeros()
  T = Float64
  dims1 = (2, 3)
  dims2a = (7,)
  dims2b = dims1
  dims2 = (dims2a, dims2b)

  a1 = rand(T, dims1)
  a2 = [(1+im)*rand(T, dims2a), (2+im)*rand(T, dims2b)]

  a1a = zeros(T, dims1)
  a2a = [zeros(Complex{T}, dims2a), zeros(Complex{T}, dims2b)]

  a1z = superzeros(a1)
  a2z = superzeros(a2)
  a1d = superzeros(dims1)
  a2d = superzeros(Complex{T}, dims2)

  @superzeros innereltype(a1) a1 a1ma dum
  @superzeros innereltype(a2) a2 a2ma dum
  @superzeros innereltype(a1) dims1 a1md dum
  @superzeros innereltype(a2) dims2 a2md dum

  ( a1a == a1z && a1a == a1d && a1a == a1ma && a1a == a1md &&
    a2a == a2z && a2a == a2d && a2a == a2ma && a2a == a2md )
end

function test_supertuplezeros(; T1=Float64, T2=Complex{Float64}, dims1=(1,), dims2=(3, 3))
  T = (T1, T2)
  dims = (dims1, dims2)
  a = superzeros(T, dims)
  ( eltype(a[1]) == T1 && eltype(a[2]) == T2 &&
    size(a[1]) == dims1 && size(a[2]) == dims2 ) 
end

function test_domainaverage(dev::Device, n)
  g = TwoDGrid(dev, n, 2π)
  X, Y = gridpoints(g)
  cx = @. cos(X)^2
  cy = @. cos(Y)^2
  0.5 ≈ FourierFlows.domainaverage(cx, g) && 0.5 ≈ FourierFlows.domainaverage(cy, g)
end

# This test could use some further work.
function test_radialspectrum(dev::Device, n, ahkl, ahρ; debug=false, atol=0.1)
  g = TwoDGrid(dev, n, 2π)
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
  isapprox(integral, parsevalsum; rtol=rtol_utils)
end

function test_parsevalsum2(func, grid; realvalued=true)
  # Compute integral in physical space
  integral = integralsquare(func, grid)
  # Compute integral in wavenumber space using Parseval's theorem
  if realvalued==true;      funch = rfft(func)
  elseif realvalued==false; funch = fft(func)
  end
  parsevalsum2 = FourierFlows.parsevalsum2(funch, grid)
  isapprox(integral, parsevalsum2; rtol=rtol_utils)
end

"""
Compute the jacobian J(a, b) and compare the result with `analytic`. Use `atol` for the comparison
to ensure validity when `analytic=0`.
"""
test_jacobian(a, b, analytic, grid) = isapprox(FourierFlows.jacobian(a, b, grid), analytic; 
                                               atol=grid.nx*grid.ny*10*eps())

"""
Test the zeros macro.
"""
function test_zeros(T=Float64, dims=(13, 45))
  a1, b1 = zeros(T, dims), zeros(T, dims)
  @zeros T dims a2 b2
  a1 == a2 && b1 == b2
end

abstract type TestVars <: AbstractVars end

physicalvars = [:a, :b]
 fouriervars = [:ah, :bh]

varspecs = cat(FourierFlows.getfieldspecs(physicalvars, :(Array{T,2})),
  FourierFlows.getfieldspecs(fouriervars, :(Array{Complex{T},2})), dims=1)

eval(FourierFlows.varsexpression(:VarsFields, physicalvars, fouriervars))
eval(FourierFlows.varsexpression(:VarsFieldsParent, physicalvars, fouriervars; parent=:TestVars))

eval(FourierFlows.varsexpression(:VarsSpecs, varspecs; typeparams=:T))
eval(FourierFlows.varsexpression(:VarsSpecsParent, varspecs; parent=:TestVars, typeparams=:T))

function test_varsexpression_fields(g::AbstractGrid{T}) where T
  @zeros T (g.nx, g.ny) a b
  @zeros Complex{T} (g.nkr, g.nl) ah bh
  v1 = VarsFields(a, b, ah, bh)

  (
     typeof(v1.a) == Array{T,2} &&
    typeof(v1.ah) == Array{Complex{T},2} &&
       size(v1.a) == size(v1.b) &&
      size(v1.ah) == size(v1.bh) &&
       size(v1.a) == (g.nx, g.ny) &&
      size(v1.ah) == (g.nkr, g.nl)
  )
end

function test_varsexpression_fields_parent(g)
  @zeros T (g.nx, g.ny) a b
  @zeros Complex{T} (g.nkr, g.nl) ah bh
  v2 = VarsFieldsParent(a, b, ah, bh)
  (
     typeof(v2.a) == Array{T,2} &&
    typeof(v2.ah) == Array{Complex{T},2} &&
       size(v2.a) == size(v2.b) &&
      size(v2.ah) == size(v2.bh) &&
       size(v2.a) == (g.nx, g.ny) &&
      size(v2.ah) == (g.nkr, g.nl)
  )
end

function test_varsexpression_specs(g::AbstractGrid{T}) where T
  @zeros T (g.nx, g.ny) a b
  @zeros Complex{T} (g.nkr, g.nl) ah bh
  v1 = VarsSpecs(a, b, ah, bh)

  (
     typeof(v1.a) == Array{T,2} &&
    typeof(v1.ah) == Array{Complex{T},2} &&
       size(v1.a) == size(v1.b) &&
      size(v1.ah) == size(v1.bh) &&
       size(v1.a) == (g.nx, g.ny) &&
      size(v1.ah) == (g.nkr, g.nl)
  )
end

function test_varsexpression_specs_parent(g)
  @zeros T (g.nx, g.ny) a b
  @zeros Complex{T} (g.nkr, g.nl) ah bh
  e2 = VarsSpecsParent(a, b, ah, bh)
  (
     typeof(v2.a) == Array{T,2} &&
    typeof(v2.ah) == Array{Complex{T},2} &&
       size(v2.a) == size(v2.b) &&
      size(v2.ah) == size(v2.bh) &&
       size(v2.a) == (g.nx, g.ny) &&
      size(v2.ah) == (g.nkr, g.nl)
  )
end

function test_supersize()
  a = rand(16, 16)
  dimsa = size(a)

  b = [rand(1), rand(3, 34)]
  dimsb = ((1,), (3, 34))

  dimsa == supersize(a) && dimsb == supersize(b)
end
