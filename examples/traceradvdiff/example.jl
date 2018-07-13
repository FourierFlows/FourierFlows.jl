using FourierFlows, PyPlot, JLD2

import FourierFlows.TracerAdvDiff

# Numerical parameters and time-stepping parameters
nx  = 128      # 2D resolution = nx^2
stepper = "RK4"   # timestepper
dt  = 0.005     # timestep
nsteps = 100 # total number of time-steps
nsubs  = 50   # number of time-steps for plotting
               # (nsteps must be multiple of nsubs)

# Physical parameters
Lx  = 2π       # domain size


gr = TwoDGrid(nx, Lx)
X, Y = gr.X, gr.Y

psiampl = 0.2
m, n = 1, 1
psiin = psiampl*cos.(m*X).*cos.(n*Y)
uin(x, y) = +psiampl*n*cos.(m*x).*sin.(n*y)
vin(x, y) = -psiampl*m*sin.(m*x).*cos.(n*y)

# uvel, vvel = 0.2, 0.1
# psiin = -uvel*Y + vvel*X
# uin(x, y) = uvel + 0*x
# vin(x, y) = vvel + 0*x


prob = TracerAdvDiff.ConstDiffSteadyFlowProblem(;
  grid = nothing,
    nx = nx,
    Lx = Lx,
    ny = nothing,
    Ly = nothing,
     κ = 0.00,
     η = nothing,
     u = uin,
     v = vin,
    dt = dt,
  stepper = stepper
  )

s, v, p, g, eq, ts = prob.state, prob.vars, prob.params, prob.grid, prob.eqn, prob.ts;


σ = 0.1
c0func(x, y) = 0.1*exp.(-(x.^2+y.^2)/(2σ^2))
# + 0*ones(size(g.X)) + 0*sin.(1/2*g.X).*sin.(1/2*g.Y)

c0 = c0func.(g.X, g.Y)
tf = nsteps*dt
cf = c0func.(g.X-uvel*tf, g.Y-vvel*tf)


TracerAdvDiff.set_c!(prob, c0)


function plot_output(prob, fig, axs; drawcolorbar=false)

  # Plot the PV field and the evolution of energy and enstrophy.

  s, v, p, g = prob.state, prob.vars, prob.params, prob.grid
  t = round(prob.state.t, 2)
  TracerAdvDiff.updatevars!(prob)
  cla()
  pcolormesh(g.X, g.Y, v.c)
  if drawcolorbar;colorbar();end
  contour(g.X, g.Y, psiin, 15, colors="k", linewidths=0.7)
  title("\$c(x,y,\\mu t= $t )\$")
  pause(0.001)
end

fig, axs = subplots(ncols=1, nrows=1, figsize=(8, 8))
plot_output(prob, fig, axs; drawcolorbar=true)

while prob.step < nsteps
  stepforward!(prob, nsubs)
  plot_output(prob, fig, axs; drawcolorbar=false)
end
