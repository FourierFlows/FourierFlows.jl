testnx(g, nx) = isapprox(g.nx, nx) 
testny(g, ny) = isapprox(g.ny, ny) 

# Physical grid tests
function testdx(g::AbstractGrid{T}) where T
  dxgrid = @. g.x[2:end] - g.x[1:end-1]
  dxones = g.dx*ones(T, size(dxgrid))
  isapprox(dxgrid, dxones)
end

function testdy(g::AbstractGrid{T}) where T
  dygrid = @. g.y[2:end] - g.y[1:end-1]
  dyones = g.dy*ones(T, size(dygrid))
  isapprox(dygrid, dyones)
end

testx(g) = isapprox(g.x[end]-g.x[1], g.Lx-g.dx)
testy(g) = isapprox(g.y[end]-g.y[1], g.Ly-g.dy)

# Test proper arrangement of fft wavenumbers. k = [0:nx/2, -nx/2+1:-1].
flippednegatives(k, mid) = -reverse(k[mid+2:end])
function testwavenumberalignment(k, nx)
  mid = Int(nx/2)
  positives = k[2:mid]
  negatives = flippednegatives(k, mid)
  isapprox(reshape(positives, mid-1), reshape(negatives, mid-1))
end

testk(g) = testwavenumberalignment(g.k, g.nx)
testl(g) = testwavenumberalignment(g.l, g.ny)
testkr(g) = isapprox(g.k[1:g.nkr], g.kr)

#testk(g) = isapprox(g.k[2:g.nkr-1], flippednegatives(g.k, g.nkr, 1))
#testk(g) = sum(g.k[2:g.nkr-1] .+ reverse(g.k[g.nkr+1:end], dims=1)) == 0.0
#testl(g) = sum(g.l[:, 2:Int(g.ny/2)] .+ reverse(g.l[:, Int(g.ny/2+2):end], dims=2)) == 0.0

function testgridpoints(g::AbstractGrid{T}) where T
  X, Y = gridpoints(g)
  dXgrid = @. X[2:end, :] - X[1:end-1, :]
  dYgrid = @. Y[:, 2:end] - Y[:, 1:end-1]
  dXones = g.dx*ones(T, size(dXgrid))
  dYones = g.dy*ones(T, size(dYgrid))
  isapprox(dXgrid, dXones) && isapprox(dYgrid, dYones)
end
