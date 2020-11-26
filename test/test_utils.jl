function test_fltype()
  Tf = Float16
  Tc = Complex{Tf}
  Tf == fltype(Tf) && Tf == fltype(Tc) && Tf == fltype((Tc, Tf)) && Tf == fltype((Tf, Tc))
end

function test_cxtype()
  Tf = Float16
  Tc = Complex{Tf}
  Tc == cxtype(Tf) && Tc == cxtype(Tc)
end

function test_innereltype(T=Complex{Float32})
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

# This test could use some further work.
function test_radialspectrum(dev::Device, n, ahkl, ahρ; debug=false, atol=0.1, rfft=false)
  grid = TwoDGrid(dev, n, 2π)
  if rfft==true
    ah = @. ahkl(grid.kr, grid.l)
  else   
    ah = @. ahkl(grid.k, grid.l)
  end
  CUDA.@allowscalar ah[1, 1] = 0.0

  ρ, ahρ_estimate = FourierFlows.radialspectrum(ah, grid; refinement=16)

  if debug
    println(sum(ahρ.(ρ) - ahρ_estimate))
    return ah, ahρ_estimate, ahρ.(ρ)
  else
    normalizeddiff = sum(abs.(ahρ.(ρ) - ahρ_estimate)) / length(ρ)
    return isapprox(normalizeddiff, 0.0, atol=atol)
  end
end

function integralsquare(func, grid::OneDGrid)
    sum(abs2.(func))*grid.dx
end

function integralsquare(func, grid::TwoDGrid)
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

"""
Test supersize().
"""
function test_supersize()
  a = rand(16, 16)
  dimsa = size(a)

  b = [rand(1), rand(3, 34)]
  dimsb = ((1,), (3, 34))

  dimsa == supersize(a) && dimsb == supersize(b)
end

"""
Test arraytype().
"""
function test_arraytype(dev::Device)
  dev==CPU() ? ArrayType(dev)==Array : ArrayType(dev)==CuArray
end

function test_arraytypeTdim(dev::Device, T=Float64, dim=1)
  dev==CPU() ? ArrayType(dev, T, dim)<:Array{T, dim} : ArrayType(dev, T, dim)<:CuArray{T, dim}
end
