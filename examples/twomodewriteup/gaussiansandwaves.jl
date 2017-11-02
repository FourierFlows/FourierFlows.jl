include("../../src/fourierflows.jl")

using FourierFlows,
      PyPlot

import FourierFlows.TwoModeBoussinesq
import FourierFlows.TwoModeBoussinesq: apv, mode1speed, mode1w,
  totalenergy, mode0energy, mode1energy, CFL

# Resolution and domain
nx, Lx = 256, 2π*1600e3

# Initial condition params
nkw   = 8         # Wave number   
alpha = 0.01      # Frequency parameter
ep    = 1e-2      # Wave amplitude
Ro    = 2e-1      # Eddy Rossby number
Reddy = Lx/20     # Eddy radius

# Simulation and computed params
f, N = 1e-4, 5e-3
sigma, kw = f*sqrt(1+alpha), 2π*nkw/Lx
m = N*kw/(f*sqrt(alpha))              # Vertical scale

# Numerical params
twave = 2π/sigma                      # Wave period
dt = 2e-2 * twave                     # Time-step
nnu0, nnu1 = 8, 8                     # Hyperviscous order
nu0 = 1e-2/(dt*(0.65π*nx/Lx)^nnu0)    # Hyperviscosity
nu1 = 1e-3/(dt*(0.65π*nx/Lx)^nnu1)    # Hyperviscosity

nsteps, nsubs = round(Int, 100twave/dt), round(Int, 2twave/dt)

prob = TwoModeBoussinesq.InitialValueProblem(
  nx=nx, Lx=Lx, nu0=nu0, nnu0=nnu0, nu1=nu1, nnu1=nnu1, f=f, N=N, m=m, dt=dt)

# Initial condition
x, y = prob.grid.X, prob.grid.Y
Z0 = f*Ro * exp.(-(x.^2+y.^2)/2Reddy^2)
TwoModeBoussinesq.set_zeta!(prob, Z0)

# ep = U/(Reddy*sigma) or U*kw/sigma
uw = minimum([ep*Reddy*sigma, ep*sigma/kw])
TwoModeBoussinesq.set_planewave!(prob, uw, nkw)


etot = Diagnostic(totalenergy, prob; nsteps=nsteps)
e0   = Diagnostic(mode0energy, prob; nsteps=nsteps)
e1   = Diagnostic(mode1energy, prob; nsteps=nsteps) 
diags = [etot, e0, e1]


# Plotting
fig, axs = subplots(ncols=2, nrows=1, figsize=(10, 5), 
  sharex=true, sharey=true)
xr, yr = x/Reddy, y/Reddy
Z00 = maximum(abs.(prob.vars.Z))
w00 = maximum(abs.(mode1w(prob)))


""" Plot the mode-0 available potential vorticity and vertical velocity. """
function makeplot!(axs, prob; eddylim=8, message=nothing, save=false)

  q = apv(prob)
  w = mode1w(prob)

  axs[1][:cla]()
  axs[2][:cla]()

  axes(axs[1])
  pcolormesh(xr, yr, q, cmap="RdBu_r",
    vmin=-Z00, vmax=Z00)

  axes(axs[2])
  pcolormesh(xr, yr, w, cmap="RdBu_r",
    vmin=-4w00, vmax=4w00)

  axs[1][:set_xlim](-eddylim, eddylim)
  axs[1][:set_ylim](-eddylim, eddylim)
  axs[2][:set_xlim](-eddylim, eddylim)
  axs[2][:set_ylim](-eddylim, eddylim)

  axs[1][:tick_params](axis="both", which="both", length=0)
  axs[2][:tick_params](axis="both", which="both", length=0)

  axs[1][:set_xlabel](L"x/R")
  axs[1][:set_ylabel](L"x/R")
  axs[2][:set_xlabel](L"x/R")

  if message != nothing
    text(0.00, 1.05, message, transform=axs[1][:transAxes], fontsize=14)
  end

  tight_layout(rect=(0.05, 0.05, 0.95, 0.95))
  pause(0.1)

  if save
    savefig(@sprintf("./plots/eddywave_%04d.png", prob.step), dpi=240)
  end

  nothing
end
 

# Run
makeplot!(axs, prob)
startwalltime = time()
while prob.step < nsteps

  stepforward!(prob, diags; nsteps=nsubs)
  TwoModeBoussinesq.updatevars!(prob)

  @printf(
    "wall: %.2f min, step: %04d, CFL: %.2f, e0: %.3f, e1: %.3f, etot: %.3f\n",
    (time()-startwalltime)/60, prob.step, CFL(prob, prob.ts.dt),
    e0.value/e0.data[1], e1.value/e1.data[1], etot.value/etot.data[1]
  )

  message = @sprintf(
    "\$t=%02d\$ wave periods, \$E_0=%.3f\$, \$E_1=%.3f\$, \$E_{\\mathrm{tot}}=%.6f\$",
    round(Int, prob.t/twave), 
    e0.value/e0.data[1], e1.value/e1.data[1], etot.value/etot.data[1])

  makeplot!(axs, prob; message=message, save=true)

end
