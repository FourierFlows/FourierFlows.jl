using 
  FourierFlows, 
  FourierFlows.TracerAdvDiff,
  PyPlot

# Numerical parameters and time-stepping parameters
nx  = 128         # 2D resolution = nx^2
stepper = "RK4"   # timestepper
dt  = 0.02        # timestep
nsubs = 200       # number of time-steps for plotting
nsteps = 4nsubs   # total number of time-steps (must be multiple of nsubs)

# Physical parameters
Lx = 2π       # domain size
kap = 0.002   # diffusivity

# Flow field
Psi = 0.2
m, n = 1, 1
uvel(x, y) = +Psi * n * cos(m*x) * sin(n*y)
vvel(x, y) = -Psi * m * sin(m*x) * cos(n*y)

# Initial condition c0 = c(x, y, t=0)
C, dc = 0.1, 0.1
x0, y0 = 1.2, 0
c0(x, y) = C*exp( -((x-x0)^2+(y-y0)^2) / (2dc^2))

# Generate problem
prob = TracerAdvDiff.ConstDiffProblem(nx=nx, Lx=Lx, kap=kap, u=uvel, v=vvel, dt=dt, stepper=stepper, steadyflow=true) 
TracerAdvDiff.set_c!(prob, c0)

# Calculate streamfunction of flow field for plotting
X, Y = prob.grid.X, prob.grid.Y
psi = @. Psi * cos(m*X) * cos(n*Y)

"Plot the flow streamlines and tracer concentration."
function plotoutput(prob, ax; drawcolorbar=false)

  TracerAdvDiff.updatevars!(prob)
  t = round(prob.state.t, digits=2)

  sca(ax)
  cla()
  pcolormesh(X, Y, prob.vars.c)

  if drawcolorbar; colorbar(); end

  contour(X, Y, psi, 15, colors="k", linewidths=0.7)
  axis("equal")
  axis("square")
  title("shading: \$c(x, y, t= $t )\$, contours: \$\\psi(x, y)\$")

  pause(0.001)

  nothing
end

fig, ax = subplots(ncols=1, nrows=1, figsize=(8, 8))
plotoutput(prob, ax; drawcolorbar=true)

# Step problem forward
while prob.step < nsteps
  stepforward!(prob, nsubs)
  plotoutput(prob, ax)
end
