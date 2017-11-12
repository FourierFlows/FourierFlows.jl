include("../../src/fourierflows.jl")

using  FourierFlows, PyPlot, JLD2
import FourierFlows.TwoModeBoussinesq

import FourierFlows.TwoModeBoussinesq: 
          updatevars!, mode1speed, set_zeta!, set_planewave!, 
          mode0energy, mode1ke, mode1pe, mode1energy 


# Turbulent background
@load "twodturb_nx0512_nnu4.jld2" Z
nx, ny = size(Z)
Lx = 2π*1600e3 

# Timestepping parameters
nsteps, substeps = 4000, 40
# Wave parameters
ep    = 1e-1     # Wave amplitude ep = U*k/sig
nkw   = 16     
alpha = 1.0 
# Computed parameters
f, N, kw   = 1e-4,                 5e-3,            2π*nkw/Lx
m, sigma   = N*kw/(f*sqrt(alpha)), f*sqrt(alpha+1)
uw0, twave = ep*sigma/kw,          2π/sigma                   # Wave period
# Time-stepping and damping
dt = twave/50
nnu0, nnu1 = 4, 8
nu0 = 5e-2/(dt*(0.65*pi*nx/Lx)^nnu0)
nu1 = 1e-2/(dt*(0.65*pi*nx/Lx)^nnu1)



# Initialize problem
prob = TwoModeBoussinesq.InitialValueProblem(nx=nx, Lx=Lx, 
  nnu0=nnu0, nnu1=nnu1, f=f, N=N, m=m, dt=dt)

# diag = Diagnostic(calc, prob; nsteps=...)
etot = Diagnostic(totalenergy, prob; nsteps=nsteps)
e0   = Diagnostic(mode0energy, prob; nsteps=nsteps)
e1   = Diagnostic(mode1energy, prob; nsteps=nsteps)
ke1  = Diagnostic(mode1ke,     prob; nsteps=nsteps)
pe1  = Diagnostic(mode1pe,     prob; nsteps=nsteps)
diags = [etot, e0, e1, ke1, pe1]


# Jet + plane wave initial condition
TwoModeBoussinesq.set_zeta!(prob, Z)
TwoModeBoussinesq.set_planewave!(prob, uw0, nkw)
update!(diags)


#message(prob) = @sprintf("\$t = %.1f\$ wave periods", prob.t/twave)


x, y = prob.grid.X, prob.grid.Y
fig, axs = subplots(ncols=2, nrows=1, sharex=true, sharey=true,
  figsize=(10, 5))


message(prob, diags) = @sprintf("step: %04d, e: %.2f, e0: %.2f, e1: %.2f",
  prob.step, diags.etot.value/diags.etot.data[1], 
  diags.e0.value/diags.e0.data[1], diags.e1.value/diags.e1.data[1])

# Step forward
while prob.step < nsteps

  stepforward!(prob, diags, nsteps=substeps)
  updatevars!(prob)

  println(message(prob, diags))

  axes(axs[1])
  pcolormesh(x, y, prob.vars.Z, cmap="RdBu_r")

  axes(axs[2])
  pcolormesh(x, y, mode1speed(prob), cmap="YlGnBu_r")

  axs[1][:xaxis][:set_visible](false)
  axs[1][:yaxis][:set_visible](false)
  axs[2][:xaxis][:set_visible](false)

  tight_layout()
  pause(0.01)

  savename = @sprintf("test_%06d.png", prob.step)
  savefig(savename, dpi=240)

end
