testnx(g, nx) = isapprox(g.nx, nx)
testny(g, ny) = isapprox(g.ny, ny)
testnz(g, nz) = isapprox(g.nz, nz)

# Physical grid tests
function testdx(g::Union{OneDGrid{T,A}, TwoDGrid{T,A}, ThreeDGrid{T,A}}) where {T, A}
  dxgrid = @. g.x[2:end] - g.x[1:end-1]
  dxones = A(g.dx*ones(T, size(dxgrid)))
  isapprox(dxgrid, dxones)
end

function testdy(g::Union{TwoDGrid{T,A}, ThreeDGrid{T,A}}) where {T, A}
  dygrid = @. g.y[2:end] - g.y[1:end-1]
  dyones = A(g.dy*ones(T, size(dygrid)))
  isapprox(dygrid, dyones)
end

function testdz(g::ThreeDGrid{T,A}) where {T, A}
  dzgrid = @. g.z[2:end] - g.z[1:end-1]
  dzones = A(g.dz*ones(T, size(dzgrid)))
  isapprox(dzgrid, dzones)
end

testx(g) = isapprox(g.x[end]-g.x[1], g.Lx-g.dx)
testy(g) = isapprox(g.y[end]-g.y[1], g.Ly-g.dy)
testz(g) = isapprox(g.z[end]-g.z[1], g.Lz-g.dz)

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

function testgridpoints(g::TwoDGrid{T}) where T
  X, Y = gridpoints(g)
  dXgrid = @. X[2:end, :] - X[1:end-1, :]
  dYgrid = @. Y[:, 2:end] - Y[:, 1:end-1]
  dXones = g.dx*ones(T, size(dXgrid))
  dYones = g.dy*ones(T, size(dYgrid))
  isapprox(dXgrid, dXones) && isapprox(dYgrid, dYones)
end

function testgridpoints(g::ThreeDGrid{T}) where T
  X, Y, Z = gridpoints(g)
  dXgrid = @. X[2:end, :, :] - X[1:end-1, :, :]
  dYgrid = @. Y[:, 2:end, :] - Y[:, 1:end-1, :]
  dZgrid = @. Z[:, :, 2:end] - Z[:, :, 1:end-1]
  dXones = g.dx*ones(T, size(dXgrid))
  dYones = g.dy*ones(T, size(dYgrid))
  dZones = g.dz*ones(T, size(dZgrid))
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
  typeof(gr.dx)==T && typeof(gr.x[1])==T && typeof(gr.Lx)==T 
end

function testtypedtwodgrid(dev::Device, nx, Lx, ny=Lx, Ly=Lx; T=Float64)
  gr = TwoDGrid(nx, Lx, ny, Ly, T=T)
  typeof(gr.dx)==T && typeof(gr.dy)==T && typeof(gr.x[1])==T && typeof(gr.y[1])==T && typeof(gr.Lx)==T && typeof(gr.Ly)==T 
end

function testtypedthreedgrid(dev::Device, nx, Lx, ny=nx, Ly=Lx, nz=nx, Lz=Lx; T=Float64)
  gr = ThreeDGrid(nx, Lx, ny, Ly, nz, Lz, T=T)
  typeof(gr.dx)==T && typeof(gr.dy)==T && typeof(gr.dz)==T && typeof(gr.x[1])==T && typeof(gr.y[1])==T && typeof(gr.z[1])==T && typeof(gr.Lx)==T && typeof(gr.Ly)==T && typeof(gr.Lz)==T 
end

function testmakefilter(dev::Device, g::AbstractGrid{T}) where T
  filter = FourierFlows.makefilter(g)
  G = typeof(g)
  nofilter = G<:OneDGrid ? filter[@. g.kr*g.dx/π < 0.65] : ( G<:TwoDGrid ? filter[@. sqrt((g.kr*g.dx/π)^2 + (g.l*g.dy/π)^2) < 0.65] : filter[@. sqrt((g.kr*g.dx/π)^2 + (g.l*g.dy/π)^2 + (g.m*g.dz/π)^2) < 0.65] )
  fullfilter = G<:OneDGrid ? filter[@. g.kr*g.dx/π > 0.999] : ( G<:TwoDGrid ? filter[@. sqrt((g.kr*g.dx/π)^2 + (g.l*g.dy/π)^2) > 0.999] : filter[@. sqrt((g.kr*g.dx/π)^2 + (g.l*g.dy/π)^2 + (g.m*g.dz/π)^2) > 0.999] )
  nofilter== zeros(dev, T, size(nofilter)) .+ 1 && isapprox(fullfilter, zeros(dev, T, size(fullfilter)), atol=1e-12)
end
