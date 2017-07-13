include("../../src/physics/barotropicqg.jl")

using BarotropicQG, PyPlot
using BarotropicQGProblems

nx  = 128
dt  = 2e-2
nu  = 1e-6
nun = 4

g, p, v, eq = twodturb(nx, nu, nun; beta=10.0)
ts = RK4TimeStepper(dt, eq.LC)

function test_plot(v)
  clf() 
  pcolormesh(g.x, g.y, v.zet);
  pause(0.01)
end

nloops = 100
nsteps = 400

fig = figure(1)
test_plot(v)

for i = 1:nloops
  @time stepforward!(v, nsteps, ts, eq, p, g)
  updatevars!(v, p, g) 
  test_plot(v)
end

