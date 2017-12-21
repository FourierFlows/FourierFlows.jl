include("../../src/FourierFlows.jl")

using FourierFlows,
      PyPlot

import FourierFlows.BarotropicQG
import FourierFlows.TwoDTurb

nx  = 128
dt  = 5e-2
nu  = 1e-6
nun = 4

g2, p2, v2, eq2 = TwoDTurb.Setups.simplenondim(nx; nu=nu, nun=nun)
ts2 = ETDRK4TimeStepper(dt, eq2.LC)

gQ, pQ, vQ, eqQ = BarotropicQG.Setups.twodturb(nx; nu=nu, nun=nun)
tsQ = ETDRK4TimeStepper(dt, eqQ.LC)

BarotropicQG.set_zeta!(vQ, pQ, gQ, v2.q)

function test_plot(axs, g, v2, vQ)
  # Make a plot that compared two-dimensional turbulence solved by
  # the TwoDTurb and BarotropicQG modules.

  axes(axs[1])
  pcolormesh(g.x, g.y, v2.q);
  #axis("square")

  axes(axs[2])
  pcolormesh(g.x, g.y, vQ.q);
  #axis("square")

  pause(0.01)
end


nloops = 100
nsteps = 100

fig, axs = subplots(ncols=2, nrows=1, sharex=true)
test_plot(axs, gQ, v2, vQ)

for i = 1:nloops
  @time stepforward!(vQ, nsteps, tsQ, eqQ, pQ, gQ)
  BarotropicQG.updatevars!(vQ, pQ, gQ)

  @time stepforward!(v2, nsteps, ts2, eq2, p2, g2)
  TwoDTurb.updatevars!(v2, p2, g2)

  test_plot(axs, gQ, v2, vQ)
end
