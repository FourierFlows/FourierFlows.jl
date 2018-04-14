# Physical grid tests
testdx(g) = abs(sum(g.x[2:end].-g.x[1:end-1]) - (g.nx-1)*g.dx) < 1e-15*(g.nx-1)
testdy(g) = abs(sum(g.y[2:end].-g.y[1:end-1]) - (g.ny-1)*g.dy) < 1e-15*(g.nx-1)

testx(g) = g.x[end]-g.x[1] == g.Lx-g.dx
testy(g) = g.y[end]-g.y[1] == g.Ly-g.dy

# Test proper arrangement of fft wavenumbers
testk(g) = sum(g.k[2:g.nkr-1] .+ flipdim(g.k[g.nkr+1:end], 1)) == 0.0
testkr(g) = sum(g.k[1:g.nkr] .- g.kr) == 0.0
testl(g) = sum(g.l[:, 2:Int(g.ny/2)] .+ flipdim(g.l[:, Int(g.ny/2+2):end], 2)) == 0.0

testX(g) = sum(g.X[end, :].-g.X[1, :]) - (g.Lx-g.dx)*g.ny < 1e-15*g.ny
testY(g) = sum(g.Y[:, end].-g.Y[:, 1]) - (g.Ly-g.dy)*g.nx < 1e-15*g.nx
testK(g) = sum(g.K[2:g.nkr-1, :] .+ flipdim(g.K[g.nkr+1:end, :], 1)) == 0.0
testL(g) = sum(g.L[:, 2:Int(g.ny/2)] .+ flipdim(g.L[:, Int(g.ny/2+2):end], 2)) == 0.0
testKr(g) = sum(g.K[1:g.nkr, :] .- g.Kr) == 0.0
testLr(g) = sum(g.L[1:g.nkr, :] .- g.Lr) == 0.0

g2 = CuTwoDGrid(32, 2π, 24, 4π)

@test testdx(g2)
@test testdy(g2)
@test testx(g2)
@test testk(g2)
@test testkr(g2)

# Cannot get tests of arrays with dimension (1, n) to work.
#@test testy(g2)
#@test testl(g2)

# Tests of 2D grid objects do not work with CuArrays.
#@test testX(g2)
#@test testY(g2)
#@test testK(g2)
#@test testL(g2)
#@test testKr(g2)
#@test testLr(g2)
