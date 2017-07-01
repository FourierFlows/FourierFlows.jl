include("../src/hydrostaticWaveEqn.jl")

import TurbOnlySolver
import Solver
import Problems

using Utils
using PyPlot

function quickplot(v, g)
  fig, axs = subplots(nrows=1, ncols=2)

  axes(axs[1])
  imshow(abs.(v.qh))

  axes(axs[2])
  imshow(v.q)

  pause(10)

  clf()
end

dt = 0.1
nx = 256
nsteps = 400
nloops = 10

# Initialize problem
g, p, v, qts, ats = Problems.test_turb_problem(nx, dt)

TurbOnlySolver.stepforward!(1, qts, v, p, g)
TurbOnlySolver.updatevars!(v, p, g)

maxsp = maximum(sqrt.(v.U.^2.0+v.V.^2.0))
@printf("CFL number: %f\n", maxsp*qts.dt/g.dx)

for i = 1:nloops

  @printf("(nx = %4d) %12s stepforward %d times: ", nx, "turb only", nsteps)
  @time TurbOnlySolver.stepforward!(nsteps, qts, v, p, g)

  TurbOnlySolver.updatevars!(v, p, g)

  fig = figure(1); clf
  imshow(v.q)
  pause(0.1)

end
