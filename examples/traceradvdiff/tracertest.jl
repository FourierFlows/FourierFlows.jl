include("../../src/fourierflows.jl")

using FourierFlows

import FourierFlows, FourierFlows.TracerAdvDiff

import FourierFlows: M0, Cx1, Cy1, Cy2, myn, Myn
import FourierFlows: Diagnostic

include("./tracerutils.jl")
       



nx, kap, Lx, Ly = 128, 1e-4, 2*pi, 1.0      # Domain
U, ep, k1       = 1.0, 0.4, 2*pi/Lx         # Specify barotropic flow

h(x)  = Ly*(1.0 - ep*sin(k1*x))             # Topography
hx(x) = -Ly*k1*ep*cos(k1*x)                 # Topographic x-gradient

ub(x, y, t) = U/h(x)                        # Barotropic x-velocity
wb(x, y, t) = y*ub(x, y, t)*hx(x)/h(x)      # Topography-following z-velocity

CFL, del, umax = 0.2, Lx/nx, U/(Ly*(1-ep))  # Time stepping parameters
dt, tfinal     = CFL*del/umax, 4.0*Ly*Lx/U  # Time-step and final time
nsteps         = ceil(Int, tfinal/dt)       # Number of steps
substeps       = ceil(Int, 0.01*nsteps)     # Number of snapshots and diags

name = @sprintf("./plots/wavychannel_ep%d_nx%04d", Int(ep*10), nx)



# Initial condition
xi, yi, dc = 0.0, -0.5*Ly, 0.05 
ci(x, y) = exp( -((x-xi)^2.0+(y-yi)^2.0)/(2.0*dc^2) ) / (2.0*pi*dc^2.0)




# Initialize problem
g = FourierFlows.TwoDGrid(nx, Lx, nx, Ly+ep; y0=-Ly-ep)
prob = TracerProblem(g, kap, ub, wb, CFL*del/umax)
TracerAdvDiff.set_c!(prob, ci)




# Diagnostics
calc_xcen(prob)   = Cx1(prob.vars.c, prob.grid)
calc_ycen(prob)   = Cy1(prob.vars.c, prob.grid)
calc_yvar(prob)   = Cy2(prob.vars.c, prob.grid)
calc_ym2(prob)    = myn(prob.vars.c, prob.grid, 2)
calc_intym2(prob) = Myn(prob.vars.c, prob.grid, 2)

xcen   = Diagnostic(calc_xcen,   prob)
ycen   = Diagnostic(calc_ycen,   prob)
yvar   = Diagnostic(calc_yvar,   prob)
ym2    = Diagnostic(calc_ym2,    prob)
intym2 = Diagnostic(calc_intym2, prob)

diags = [xcen, ycen, yvar, ym2, intym2]




# Plotting
mask = Array{Bool}(g.nx, g.ny)
for j=1:g.ny, i=1:g.nx
  mask[i, j] = g.y[j] < -h(g.x[i]) ? true : false
end


t_theory = linspace(0.0, tfinal, 100)
function makeplot(axs, prob)
  vs, pr, g = FourierFlows.unpack(prob)
  cplot = NullableArray(vs.c, mask)

  flatvar = 2.0*kap*t_theory + yvar.data[1]

  axs[1][:cla]()
  axs[2][:cla]()

  axes(axs[1])
  pcolormesh(g.X, g.Y, PyObject(cplot))
  xlabel(L"$x$"); ylabel(L"$z$")

  axes(axs[2])
  plot(t_theory, flatvar;    color="C0", linestyle="-")
  plot(yvar.time, yvar.data; color="C1", linestyle="--")
  xlabel(L"$t$")
  ylabel(L"C_2 = \int (y-\bar y)^2 c \mathrm{d}y \mathrm{d}x")
    
  xlim(0.0, tfinal)

  # Formatting
  axs[1][:xaxis][:set_ticks_position]("top")
  axs[1][:xaxis][:set_label_position]("top")

  pause(0.01)
  nothing
end

fig, axs = subplots(nrows=2)




# Integrate
while prob.t < tfinal

  stepforward!(prob; nsteps=substeps)
  TracerAdvDiff.updatevars!(prob)
  increment!(diags)

  flatvar = 2.0*kap*yvar.time[end] + yvar.data[1]

  @printf(
  "step: %04d, xc: %.2f, yc: %.2f, sigma (flat): %.2e, sigma (wavy): %.2e\n", 
  prob.step, xcen.value, ycen.value, flatvar, yvar.value)

  makeplot(axs, prob)
  savefig(@sprintf("%s_%06d.png", name, prob.step), dpi=240)

end
