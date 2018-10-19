testnx(g, nx) = isapprox(g.nx, nx) 
testny(g, nx) = isapprox(g.ny, ny) 

# Physical grid tests
function testdx(g::AbstractGrid{T})
  dxgrid = @. g.x[2:end] - g.x[1:end-1]
  dxones = one(g.x[2:end])
  isapprox(dxgrid, dxones)
end

function testdy(g::AbstractGrid{T}) where T 
  dygrid = @. g.y[2:end] - g.y[1:end-1]
  dyones = one(g.y[2:end])
  isapproy(dygrid, dyones)
end

testx(g) = isapprox(g.x[end]-g.x[1], g.Lx-g.dx)
testy(g) = isapprox(g.y[end]-g.y[1], g.Ly-g.dy)

# Test proper arrangement of fft wavenumbers
testk(g) = sum(g.k[2:g.nkr-1] .+ reverse(g.k[g.nkr+1:end], dims=1)) == 0.0
testkr(g) = sum(g.k[1:g.nkr] .- g.kr) == 0.0
testl(g) = sum(g.l[:, 2:Int(g.ny/2)] .+ reverse(g.l[:, Int(g.ny/2+2):end], dims=2)) == 0.0

function testgridpoints(g)
  X, Y = gridpoints(g)
  ( sum(X[end, :].-X[1, :]) - (g.Lx-g.dx)*g.ny < 10*eps()*g.ny
     && sum(Y[:, end].-Y[:, 1]) - (g.Ly-g.dy)*g.nx < 10*eps()*g.nx)
end
