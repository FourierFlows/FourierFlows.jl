include("../../src/fourierflows.jl")

using FourierFlows,
      PyPlot

import FourierFlows.BarotropicQG

nx  = 128
dt  = 1e-1
nu  = 1e-6
nun = 4

#g, p, v, eq = BarotropicQG.Setups.twodturb_withmountains(
#g, p, v, eq = BarotropicQG.Setups.gaussian_topo_turb(
#  nx; nu=nu, nun=nun, topotype="mountain", hbump=0.4)
g, p, v, eq = BarotropicQG.Setups.forced_gaussian_topo_turb(
  nx; nu=nu, nun=nun, topotype="mountain", hbump=0.4)
#ts = ForwardEulerTimeStepper(dt, eq.LC)
ts = ETDRK4TimeStepper(dt, eq.LC)


function test_plot(g, v)
  # Make a plot that compared two-dimensional turbulence solved by
  # the TwoDTurb and BarotropicQG modules.
  pcolormesh(g.x, g.y, v.zeta);
  axis("square")
  pause(0.01)
end

nloops = 100
nsteps = 200

fig = figure()
test_plot(g, v)

for i = 1:nloops
  @time stepforward!(v, ts, eq, p, g; nsteps=nsteps)
  BarotropicQG.updatevars!(v, p, g)
  test_plot(g, v)
end
