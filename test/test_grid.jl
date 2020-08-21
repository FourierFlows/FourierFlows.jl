testnx(g, nx) = isapprox(g.nx, nx)
testny(g, ny) = isapprox(g.ny, ny)
testnz(g, nz) = isapprox(g.nz, nz)

# Physical grid tests
function testdx(dev, g::Union{OneDGrid{T}, TwoDGrid{T}, ThreeDGrid{T}}) where T
  dxgrid = @. g.x[2:end] - g.x[1:end-1]
  dxones = ArrayType(dev)(g.dx*ones(T, size(dxgrid)))
  isapprox(dxgrid, dxones)
end

function testdy(dev, g::Union{TwoDGrid{T}, ThreeDGrid{T}}) where T
  dygrid = @. g.y[2:end] - g.y[1:end-1]
  dyones = ArrayType(dev)(g.dy*ones(T, size(dygrid)))
  isapprox(dygrid, dyones)
end

function testdz(dev, g::ThreeDGrid{T}) where T
  dzgrid = @. g.z[2:end] - g.z[1:end-1]
  dzones = ArrayType(dev)(g.dz*ones(T, size(dzgrid)))
  isapprox(dzgrid, dzones)
end

testx(g) = isapprox(g.x[end]-g.x[1], g.Lx-g.dx)
testy(g) = isapprox(g.y[end]-g.y[1], g.Ly-g.dy)
testz(g) = isapprox(g.z[end]-g.z[1], g.Lz-g.dz)

testdk(g::Union{OneDGrid, TwoDGrid, ThreeDGrid}) = isapprox(g.k[2], 2π/g.Lx)
testdl(g::Union{TwoDGrid, ThreeDGrid}) = isapprox(g.l[2], 2π/g.Ly)
testdm(g::ThreeDGrid) = isapprox(g.m[2], 2π/g.Lz)

# Test proper arrangement of fft wavenumbers. k = [0:nx/2, -nx/2+1:-1].
flippednegatives(k, mid) = -reverse(k[mid+2:end])
function testwavenumberalignment(k, nx)
  mid = Int(nx/2)
  positives = k[2:mid]
  negatives = flippednegatives(k, mid)
  isapprox(reshape(positives, mid-1), reshape(negatives, mid-1))
end

testk(g) = testwavenumberalignment(g.k, g.nx)
testl(g::Union{TwoDGrid, ThreeDGrid}) = testwavenumberalignment(g.l, g.ny)
testm(g::ThreeDGrid) = testwavenumberalignment(g.m, g.nz)
testkr(g) = isapprox(cat(g.k[1:g.nkr-1], abs(g.k[g.nkr]), dims=1), g.kr)

function testgridpoints(dev::Device, g::TwoDGrid{T}) where T
  X, Y = gridpoints(g)
  dXgrid = @. X[2:end, :] - X[1:end-1, :]
  dYgrid = @. Y[:, 2:end] - Y[:, 1:end-1]
  dXones = ArrayType(dev)(g.dx*ones(T, size(dXgrid)))
  dYones = ArrayType(dev)(g.dy*ones(T, size(dYgrid)))
  isapprox(dXgrid, dXones) && isapprox(dYgrid, dYones)
end

function testgridpoints(dev::Device, g::ThreeDGrid{T}) where T
  X, Y, Z = gridpoints(g)
  dXgrid = @. X[2:end, :, :] - X[1:end-1, :, :]
  dYgrid = @. Y[:, 2:end, :] - Y[:, 1:end-1, :]
  dZgrid = @. Z[:, :, 2:end] - Z[:, :, 1:end-1]
  dXones = ArrayType(dev)(g.dx*ones(T, size(dXgrid)))
  dYones = ArrayType(dev)(g.dy*ones(T, size(dYgrid)))
  dZones = ArrayType(dev)(g.dz*ones(T, size(dZgrid)))
  isapprox(dXgrid, dXones) && isapprox(dYgrid, dYones) && isapprox(dZgrid, dZones)
end

function testdealias(g::OneDGrid)
  T = typeof(g.Lx)
  fh = ones(Complex{T}, size(g.kr))
  dealias!(fh, g)
  kmax = round(maximum(g.kr)*2/3)
  isapprox(sum(abs.(fh[ g.kr .>= kmax ])), 0)
end

function testdealias(g::TwoDGrid)
  T = typeof(g.Lx)
  fh = ones(Complex{T}, size(g.Krsq))
  dealias!(fh, g)
  kmax = round(maximum(g.kr)*2/3)
  lmax = floor(maximum(g.l)*2/3)

  temp = 0.0
  for j = 1:g.nl, i = 1:g.nkr
    if ( g.kr[i] >= kmax ) & ( g.l[j] >= lmax || g.l[j] < -lmax )
      temp += abs.(fh[i, j]) #temp = sum of |fh| for aliased wavenumbers
    end
  end
  isapprox(temp, 0)
end

function testdealias(g::ThreeDGrid)
  T = typeof(g.Lx)
  fh = ones(Complex{T}, size(g.Krsq))
  dealias!(fh, g)
  kmax = round(maximum(g.kr)*2/3)
  lmax = floor(maximum(g.l)*2/3)
  mmax = floor(maximum(g.m)*2/3)

  temp = 0.0
  for k = 1:g.nm, j = 1:g.nl, i = 1:g.nkr
    if ( g.kr[i] >= kmax ) & ( g.l[j] >= lmax || g.l[j] < -lmax ) & ( g.m[k] >= mmax || g.m[k] < -mmax )
      temp += abs.(fh[i, j, k]) #temp = sum of |fh| for aliased wavenumbers
    end
  end
  isapprox(temp, 0)
end

function testtypedonedgrid(dev::Device, nx, Lx; T=Float64)
  gr = OneDGrid(nx, Lx, T=T)
  typeof(gr.dx)==T && typeof(gr.x[1])==T && typeof(gr.Lx)==T && eltype(gr) == T
end

function testtypedtwodgrid(dev::Device, nx, Lx, ny=Lx, Ly=Lx; T=Float64)
  gr = TwoDGrid(nx, Lx, ny, Ly, T=T)
  typeof(gr.dx)==T && typeof(gr.dy)==T && typeof(gr.x[1])==T && typeof(gr.y[1])==T && typeof(gr.Lx)==T && typeof(gr.Ly)==T && eltype(gr) == T
end

function testtypedthreedgrid(dev::Device, nx, Lx, ny=nx, Ly=Lx, nz=nx, Lz=Lx; T=Float64)
  gr = ThreeDGrid(nx, Lx, ny, Ly, nz, Lz, T=T)
  typeof(gr.dx)==T && typeof(gr.dy)==T && typeof(gr.dz)==T && typeof(gr.x[1])==T && typeof(gr.y[1])==T && typeof(gr.z[1])==T && typeof(gr.Lx)==T && typeof(gr.Ly)==T && typeof(gr.Lz)==T && eltype(gr) == T
end

function testmakefilter(dev::Device, g::AbstractGrid)
  T = eltype(g)
  filter = FourierFlows.makefilter(g)
  G = typeof(g)
  nofilter = G<:OneDGrid ? filter[@. g.kr*g.dx/π < 0.65] : ( G<:TwoDGrid ? filter[@. sqrt((g.kr*g.dx/π)^2 + (g.l*g.dy/π)^2) < 0.65] : filter[@. sqrt((g.kr*g.dx/π)^2 + (g.l*g.dy/π)^2 + (g.m*g.dz/π)^2) < 0.65] )
  fullfilter = G<:OneDGrid ? filter[@. g.kr*g.dx/π > 0.999] : ( G<:TwoDGrid ? filter[@. sqrt((g.kr*g.dx/π)^2 + (g.l*g.dy/π)^2) > 0.999] : filter[@. sqrt((g.kr*g.dx/π)^2 + (g.l*g.dy/π)^2 + (g.m*g.dz/π)^2) > 0.999] )
  nofilter== zeros(dev, T, size(nofilter)) .+ 1 && isapprox(fullfilter, zeros(dev, T, size(fullfilter)), atol=1e-12)
end

function test_plan_flows_fftrfft(::CPU; T=Float64)
  A = ArrayType(CPU())
  return ( typeof(FourierFlows.plan_flows_fft(A(rand(Complex{T}, (4,))))) <: FFTW.cFFTWPlan{Complex{T},-1,false,1} &&
  typeof(FourierFlows.plan_flows_fft(A(rand(Complex{T}, (4, 6))))) <: FFTW.cFFTWPlan{Complex{T},-1,false,2} &&
  typeof(FourierFlows.plan_flows_fft(A(rand(Complex{T}, (4, 6, 8))))) <: FFTW.cFFTWPlan{Complex{T},-1,false,3} &&
  FourierFlows.plan_flows_fft(A(rand(Complex{T}, (4, 6, 8))), [1, 2]).region == [1, 2] &&
  typeof(FourierFlows.plan_flows_rfft(A(rand(T, (4,))))) <: FFTW.rFFTWPlan{T,-1,false,1} &&
  typeof(FourierFlows.plan_flows_rfft(A(rand(T, (4, 6))))) <: FFTW.rFFTWPlan{T,-1,false,2} &&
  typeof(FourierFlows.plan_flows_rfft(A(rand(T, (4, 6, 8))))) <: FFTW.rFFTWPlan{T,-1,false,3} &&
  FourierFlows.plan_flows_rfft(A(rand(T, (4, 6, 8))), [1, 2]).region == [1, 2] )
end

function test_plan_flows_fftrfft(::GPU; T=Float64)
  A = ArrayType(GPU())
  return ( typeof(FourierFlows.plan_flows_fft(A(rand(Complex{T}, (4,))))) == CuArrays.CUFFT.cCuFFTPlan{Complex{T},-1,false,1} &&
  typeof(FourierFlows.plan_flows_fft(A(rand(Complex{T}, (4, 6))))) == CuArrays.CUFFT.cCuFFTPlan{Complex{T},-1,false,2} &&
  typeof(FourierFlows.plan_flows_fft(A(rand(Complex{T}, (4, 6, 8))))) == CuArrays.CUFFT.cCuFFTPlan{Complex{T},-1,false,3} &&
  FourierFlows.plan_flows_fft(A(rand(Complex{T}, (4, 6, 8))), [1, 2]).region == [1, 2] &&
  typeof(FourierFlows.plan_flows_rfft(A(rand(T, (4,))))) == CuArrays.CUFFT.rCuFFTPlan{T,-1,false,1} &&
  typeof(FourierFlows.plan_flows_rfft(A(rand(T, (4, 6))))) == CuArrays.CUFFT.rCuFFTPlan{T,-1,false,2} &&
  typeof(FourierFlows.plan_flows_rfft(A(rand(T, (4, 6, 8))))) == CuArrays.CUFFT.rCuFFTPlan{T,-1,false,3} &&
  FourierFlows.plan_flows_rfft(A(rand(T, (4, 6, 8))), [1, 2]).region == [1, 2] )
end