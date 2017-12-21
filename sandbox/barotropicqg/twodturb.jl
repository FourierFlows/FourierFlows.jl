include("../../src/FourierFlows.jl")

using FourierFlows,
      PyPlot

import FourierFlows.TwoDTurb
import FourierFlows.BarotropicQG, BarotropicQGProblems

nx  = 128
dt  = 1e-2
nu  = 1e-4
nun = 4

g, p1, v1, eq1 = TwoDTurbProblems.simplenondim(nx; nu=nu, nun=nun)
g, p2, v2, eq2 = BarotropicQGProblems.twodturb(nx; nu=nu, nun=nun)

ts1 = TimeSteppers.RK4TimeStepper(dt, eq1.LC)
ts2 = TimeSteppers.RK4TimeStepper(dt, eq2.LC)

function test_plot(axs, g, v1, v2)
  # Make a plot that compared two-dimensional turbulence solved by
  # the TwoDTurb and BarotropicQG modules.

  #clf()

  axes(axs[1]); cla()
  pcolormesh(g.x, g.y, v1.q);
  axis("square")

  axes(axs[2])
  pcolormesh(g.x, g.y, v2.q);
  axis("square")

  pause(0.01)
end


nloops = 100
nsteps = 100

fig, axs = subplots(ncols=2, nrows=1)
test_plot(axs, g, v1, v2)

for i = 1:nloops
  #@time TimeSteppers.stepforward!(v1, nsteps, ts1, eq1, p1, g)
  @time BarotropicQG.stepforward!(v2, nsteps, ts2, eq2, p2, g)

  updatevars!(v1, p1, g)
  updatevars!(v2, p2, g)

  test_plot(axs, g, v1, v2)
end
