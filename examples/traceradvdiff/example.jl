using FourierFlows, PyPlot, JLD2

import FourierFlows.TracerAdvDiff

# Numerical parameters and time-stepping parameters
nx  = 128      # 2D resolution = nx^2
stepper = "FilteredRK4"   # timestepper
dt  = 0.1     # timestep
nsteps = 50 # total number of time-steps
nsubs  = 1   # number of time-steps for plotting
               # (nsteps must be multiple of nsubs)

# Physical parameters
Lx  = 2π       # domain size

uin(x, y) = -0.1*sin.(1*x).*cos.(1*y)
vin(x, y) =  0.1*cos.(1*x).*sin.(1*y)


prob = TracerAdvDiff.ConstDiffSteadyFlowProblem(;
  grid = nothing,
    nx = nx,
    Lx = Lx,
    ny = nothing,
    Ly = nothing,
     κ = 1.0,
     η = nothing,
     u = uin,
     v = vin,
    dt = dt,
  stepper = stepper
  )

s, v, p, g, eq, ts = prob.state, prob.vars, prob.params, prob.grid, prob.eqn, prob.ts;

TracerAdvDiff.set_c!(prob, 0.01*ones(size(g.X))+0.0*sin.(2*g.X).*sin.(2*g.Y))


function plot_output(prob, fig, axs; drawcolorbar=false)

  # Plot the PV field and the evolution of energy and enstrophy.

  s, v, p, g = prob.state, prob.vars, prob.params, prob.grid
  t = round(prob.state.t, 2)
  TracerAdvDiff.updatevars!(prob)
  cla()
  pcolormesh(g.X, g.Y, v.c)
  # pcolormesh(fftshift(g.Kr), fftshift(g.Lr), fftshift(abs.(prob.state.sol)))
  title("\$c(x,y,\\mu t= $t )\$")
  if drawcolorbar;colorbar();end
  pause(0.001)
end

fig, axs = subplots(ncols=1, nrows=1, figsize=(8, 8))
plot_output(prob, fig, axs; drawcolorbar=true)

while prob.step < nsteps
  stepforward!(prob, nsubs)
  plot_output(prob, fig, axs; drawcolorbar=false)
end
