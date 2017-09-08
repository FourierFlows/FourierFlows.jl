include("../../src/physics/barotropicqg.jl")

using BarotropicQG, PyPlot
using BarotropicQGProblems

nx = 256
dt = 1e-3

g, p, v, eq = nondimensional_southern_ocean(nx, 0.1, 1.0, 0.0, 10.0)
ts = AB3TimeStepper(dt, eq.LC)

function test_plot(v)
  @printf("U = %f\n", v.U)
  #clf; imshow(v.zet, cmap="RdBu_r")
  clf; imshow(v.sp);
  pause(0.01)
end

nloops = 100
nsteps = 10

fig = figure(1)
test_plot(v)

for i = 1:nloops
  stepforward!(v, nsteps, ts, eq, p, g)
  updatevars!(v, g) 
  test_plot(v)
end
