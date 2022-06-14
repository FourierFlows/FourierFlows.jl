testnx(g, nx) = isapprox(g.nx, nx)
testny(g, ny) = isapprox(g.ny, ny)
testnz(g, nz) = isapprox(g.nz, nz)

# Physical grid tests
function testdx(dev, g::Union{OneDGrid{T}, TwoDGrid{T}, ThreeDGrid{T}}) where T
  dxgrid = @. g.x[2:end] - g.x[1:end-1]
  dxones = ArrayType(dev)(g.dx*ones(T, size(dxgrid)))
  
  return isapprox(dxgrid, dxones)
end

function testdy(dev, g::Union{TwoDGrid{T}, ThreeDGrid{T}}) where T
  dygrid = @. g.y[2:end] - g.y[1:end-1]
  dyones = ArrayType(dev)(g.dy*ones(T, size(dygrid)))
  
  return isapprox(dygrid, dyones)
end

function testdz(dev, g::ThreeDGrid{T}) where T
  dzgrid = @. g.z[2:end] - g.z[1:end-1]
  dzones = ArrayType(dev)(g.dz*ones(T, size(dzgrid)))
  
  return isapprox(dzgrid, dzones)
end

testx(g) = isapprox(g.x[end]-g.x[1], g.Lx-g.dx)
testy(g) = isapprox(g.y[end]-g.y[1], g.Ly-g.dy)
testz(g) = isapprox(g.z[end]-g.z[1], g.Lz-g.dz)

testdk(g::Union{OneDGrid, TwoDGrid, ThreeDGrid}) = CUDA.@allowscalar isapprox(g.k[2], 2π/g.Lx)
testdl(g::Union{TwoDGrid, ThreeDGrid}) = CUDA.@allowscalar isapprox(g.l[2], 2π/g.Ly)
testdm(g::ThreeDGrid) = CUDA.@allowscalar isapprox(g.m[2], 2π/g.Lz)

# Test proper arrangement of fft wavenumbers. k = [0:nx/2, -nx/2+1:-1].
flippednegatives(k, mid) = -reverse(k[mid+2:end])
function testwavenumberalignment(k, nx)
  mid = Int(nx/2)
  positives = k[2:mid]
  negatives = flippednegatives(k, mid)
  
  return isapprox(reshape(positives, mid-1), reshape(negatives, mid-1))
end

testk(g) = testwavenumberalignment(g.k, g.nx)
testl(g::Union{TwoDGrid, ThreeDGrid}) = testwavenumberalignment(g.l, g.ny)
testm(g::ThreeDGrid) = testwavenumberalignment(g.m, g.nz)
testkr(g) = CUDA.@allowscalar isapprox(cat(g.k[1:g.nkr-1], abs(g.k[g.nkr]), dims=1), g.kr)

function testgridpoints(dev::Device, g::OneDGrid{T}) where T
  X = gridpoints(g)
  dXgrid = @. X[2:end, :] - X[1:end-1, :]
  dXones = ArrayType(dev)(g.dx*ones(T, size(dXgrid)))

  return isapprox(dXgrid, dXones)
end

function testgridpoints(dev::Device, g::TwoDGrid{T}) where T
  X, Y = gridpoints(g)
  dXgrid = @. X[2:end, :] - X[1:end-1, :]
  dYgrid = @. Y[:, 2:end] - Y[:, 1:end-1]
  dXones = ArrayType(dev)(g.dx*ones(T, size(dXgrid)))
  dYones = ArrayType(dev)(g.dy*ones(T, size(dYgrid)))
  
  return isapprox(dXgrid, dXones) && isapprox(dYgrid, dYones)
end

function testgridpoints(dev::Device, g::ThreeDGrid{T}) where T
  X, Y, Z = gridpoints(g)
  dXgrid = @. X[2:end, :, :] - X[1:end-1, :, :]
  dYgrid = @. Y[:, 2:end, :] - Y[:, 1:end-1, :]
  dZgrid = @. Z[:, :, 2:end] - Z[:, :, 1:end-1]
  dXones = ArrayType(dev)(g.dx*ones(T, size(dXgrid)))
  dYones = ArrayType(dev)(g.dy*ones(T, size(dYgrid)))
  dZones = ArrayType(dev)(g.dz*ones(T, size(dZgrid)))
  
  return isapprox(dXgrid, dXones) && isapprox(dYgrid, dYones) && isapprox(dZgrid, dZones)
end

function testdealias(grid::OneDGrid)
  fh = ones(Complex{eltype(grid)}, size(grid.kr))
  dealias!(fh, grid)
  
  kmax = round(maximum(grid.kr)*2/3)
  
  for i₁ = 1:grid.nkr
    if CUDA.@allowscalar grid.kr[i₁] < kmax
      fh[i₁] = 0
    end
  end

  return isapprox(sum(abs.(fh)), 0)
end

function testdealias(grid::TwoDGrid)
  fh = ones(Complex{eltype(grid)}, size(grid.Krsq))
  dealias!(fh, grid)
  
  kmax = round(maximum(grid.kr)*2/3)
  lmax = round(maximum(abs.(grid.l))*2/3)

  for i₂ = 1:grid.nl, i₁ = 1:grid.nkr
    if (
        (CUDA.@allowscalar grid.kr[i₁] < kmax ) &&
        (CUDA.@allowscalar grid.l[i₂] < lmax && CUDA.@allowscalar grid.l[i₂] ≥ -lmax)
        )
      fh[i₁, i₂] = 0
    end
  end
  
  return isapprox(sum(abs.(fh)), 0)
end

function testdealias(grid::ThreeDGrid)
  fh = ones(Complex{eltype(grid)}, size(grid.Krsq))
  dealias!(fh, grid)
  
  kmax = round(maximum(grid.kr)*2/3)
  lmax = round(maximum(abs.(grid.l))*2/3)
  mmax = round(maximum(abs.(grid.m))*2/3)

  for i₃ = 1:grid.nm, i₂ = 1:grid.nl, i₁ = 1:grid.nkr
    if (
        (CUDA.@allowscalar grid.kr[i₁] < kmax ) &&
        (CUDA.@allowscalar grid.l[i₂] < lmax && CUDA.@allowscalar grid.l[i₂] ≥ -lmax) &&
        (CUDA.@allowscalar grid.m[i₃] < mmax && CUDA.@allowscalar grid.m[i₃] ≥ -mmax)
        )
      fh[i₁, i₂, i₃] = 0
    end
  end
  
  return isapprox(sum(abs.(fh)), 0)
end

function testnodealias(grid::OneDGrid)
  fh = ones(Complex{eltype(grid)}, size(grid.kr))
  return dealias!(fh, grid) == nothing
end

function testnodealias(grid::Union{TwoDGrid, ThreeDGrid})
  fh = ones(Complex{eltype(grid)}, size(grid.Krsq))
  return dealias!(fh, grid) == nothing
end

function testtypedonedgrid(dev::Device, nx, Lx; T=Float64)
  grid = OneDGrid(nx, Lx, T=T)
  
  return typeof(grid.dx)==T && typeof(grid.x[1])==T && typeof(grid.Lx)==T && eltype(grid) == T
end

function testtypedtwodgrid(dev::Device, nx, Lx, ny=Lx, Ly=Lx; T=Float64)
  grid = TwoDGrid(nx, Lx, ny, Ly, T=T)
  
  return typeof(grid.dx)==T && typeof(grid.dy)==T && typeof(grid.x[1])==T && typeof(grid.y[1])==T && typeof(grid.Lx)==T && typeof(grid.Ly)==T && eltype(grid) == T
end

function testtypedthreedgrid(dev::Device, nx, Lx, ny=nx, Ly=Lx, nz=nx, Lz=Lx; T=Float64)
  grid = ThreeDGrid(nx, Lx, ny, Ly, nz, Lz, T=T)
  
  return typeof(grid.dx)==T && typeof(grid.dy)==T && typeof(grid.dz)==T && typeof(grid.x[1])==T && typeof(grid.y[1])==T && typeof(grid.z[1])==T && typeof(grid.Lx)==T && typeof(grid.Ly)==T && typeof(grid.Lz)==T && eltype(grid) == T
end

function testmakefilter(dev::Device, g::AbstractGrid)
  T = eltype(g)
  filter = FourierFlows.makefilter(g)
  G = typeof(g)
  nofilter = G<:OneDGrid ? filter[@. g.kr*g.dx/π < 0.65] : ( G<:TwoDGrid ? filter[@. sqrt((g.kr*g.dx/π)^2 + (g.l*g.dy/π)^2) < 0.65] : filter[@. sqrt((g.kr*g.dx/π)^2 + (g.l*g.dy/π)^2 + (g.m*g.dz/π)^2) < 0.65] )
  fullfilter = G<:OneDGrid ? filter[@. g.kr*g.dx/π > 0.999] : ( G<:TwoDGrid ? filter[@. sqrt((g.kr*g.dx/π)^2 + (g.l*g.dy/π)^2) > 0.999] : filter[@. sqrt((g.kr*g.dx/π)^2 + (g.l*g.dy/π)^2 + (g.m*g.dz/π)^2) > 0.999] )
  
  return nofilter== zeros(dev, T, size(nofilter)) .+ 1 && isapprox(fullfilter, zeros(dev, T, size(fullfilter)), atol=1e-12)
end

function test_plan_flows_fftrfft(::CPU; T=Float64)
  A = ArrayType(CPU())
  return (typeof(FourierFlows.plan_flows_fft(A(rand(Complex{T}, (4,))))) <: FFTW.cFFTWPlan{Complex{T},-1,false,1} &&
  typeof(FourierFlows.plan_flows_fft(A(rand(Complex{T}, (4, 6))))) <: FFTW.cFFTWPlan{Complex{T},-1,false,2} &&
  typeof(FourierFlows.plan_flows_fft(A(rand(Complex{T}, (4, 6, 8))))) <: FFTW.cFFTWPlan{Complex{T},-1,false,3} &&
  FourierFlows.plan_flows_fft(A(rand(Complex{T}, (4, 6, 8))), [1, 2]).region == [1, 2] &&
  typeof(FourierFlows.plan_flows_rfft(A(rand(T, (4,))))) <: FFTW.rFFTWPlan{T,-1,false,1} &&
  typeof(FourierFlows.plan_flows_rfft(A(rand(T, (4, 6))))) <: FFTW.rFFTWPlan{T,-1,false,2} &&
  typeof(FourierFlows.plan_flows_rfft(A(rand(T, (4, 6, 8))))) <: FFTW.rFFTWPlan{T,-1,false,3} &&
  
  return FourierFlows.plan_flows_rfft(A(rand(T, (4, 6, 8))), [1, 2]).region == [1, 2])
end

function test_plan_flows_fftrfft(::GPU; T=Float64)
  A = ArrayType(GPU())
  return (typeof(FourierFlows.plan_flows_fft(A(rand(Complex{T}, (4,))))) == CUDA.CUFFT.cCuFFTPlan{Complex{T},-1,false,1} &&
  typeof(FourierFlows.plan_flows_fft(A(rand(Complex{T}, (4, 6))))) == CUDA.CUFFT.cCuFFTPlan{Complex{T},-1,false,2} &&
  typeof(FourierFlows.plan_flows_fft(A(rand(Complex{T}, (4, 6, 8))))) == CUDA.CUFFT.cCuFFTPlan{Complex{T},-1,false,3} &&
  FourierFlows.plan_flows_fft(A(rand(Complex{T}, (4, 6, 8))), [1, 2]).region == [1, 2] &&
  typeof(FourierFlows.plan_flows_rfft(A(rand(T, (4,))))) == CUDA.CUFFT.rCuFFTPlan{T,-1,false,1} &&
  typeof(FourierFlows.plan_flows_rfft(A(rand(T, (4, 6))))) == CUDA.CUFFT.rCuFFTPlan{T,-1,false,2} &&
  typeof(FourierFlows.plan_flows_rfft(A(rand(T, (4, 6, 8))))) == CUDA.CUFFT.rCuFFTPlan{T,-1,false,3} &&
  
  return FourierFlows.plan_flows_rfft(A(rand(T, (4, 6, 8))), [1, 2]).region == [1, 2])
end

function test_aliased_fraction(dev, aliased_fraction)
  nx, Lx = 16, 2π
  ny, Ly = 32, 2π
  nz, Lz = 34, 2π
  
  g₁ = OneDGrid(nx, Lx; aliased_fraction = aliased_fraction)
  g₂ = TwoDGrid(nx, Lx, ny, Ly; aliased_fraction = aliased_fraction)
  g₃ = ThreeDGrid(nx, Lx, ny, Ly, nz, Lz; aliased_fraction = aliased_fraction)

  lower_end(n) = floor(Int, (1 - aliased_fraction)/2 * n) + 1
  upper_end(n) = ceil(Int, (1 + aliased_fraction)/2 * n)
  upper_end_r(n) = Int(n/2)+1
  
  kralias = aliased_fraction==0 ? nothing : lower_end(nx):upper_end_r(nx)
  kalias = aliased_fraction==0 ? nothing : lower_end(nx):upper_end(nx)
  lalias = aliased_fraction==0 ? nothing : lower_end(ny):upper_end(ny)
  malias = aliased_fraction==0 ? nothing : lower_end(nz):upper_end(nz)
    
  return (g₁.aliased_fraction == aliased_fraction && g₂.aliased_fraction == aliased_fraction && 
  g₃.aliased_fraction == aliased_fraction && g₁.kralias == kralias && g₁.kalias == kalias && 
  g₂.kralias == kralias && g₂.kalias == kalias && g₂.lalias == lalias && g₃.kralias == kralias &&
  g₃.kalias == kalias && g₃.lalias == lalias && g₃.malias == malias)
end
