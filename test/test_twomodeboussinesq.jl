include("../src/fourierflows.jl")
include("./twomodeutils.jl")

using FourierFlows,
      PyPlot

import FourierFlows.TwoModeBoussinesq

rms(a) = sqrt(mean(a.^2))

function niwqgplot(axs, vs, pr, g, q0, Uw, R, tnd) 

  TwoModeBoussinesq.updatevars!(vs, pr, g)

  domfrac = 5
  xl, xr = -g.Lx/domfrac, g.Lx/domfrac
  yl, yr = -g.Ly/domfrac, g.Ly/domfrac

  KE = 0.5*(abs2.(vs.u) + abs2.(vs.v))

  qw = calc_qw(vs.u, vs.v, pr.f, pr, g)

  # All lengths non-dimensionalized by R
  axes(axs[1])
  pcolormesh(g.X/R, g.Y/R, vs.q*tnd, cmap="RdBu_r",
    vmin=-q0, vmax=q0)
  xlim(xl/R, xr/R); ylim(yl/R, yr/R)
  xlabel(L"x/R"); ylabel(L"y/R")
  title(L"q")

  axes(axs[2])
  pcolormesh(g.X/R, g.Y/R, KE, cmap="YlGnBu_r",
    vmin=0, vmax=0.5*Uw^2)
  xlim(xl/R, xr/R); ylim(yl/R, yr/R)
  xlabel(L"x/R")
  title(L"\sqrt{u^2+v^2}")

  E10 = 0.5*(Uw/2)^2*g.Lx^2
  E1 = 0.5*sum(abs2.(vs.uh)) + 0.5*sum(abs2.(vs.vh))

  @printf("rms Ro: %.2e, max speed: %.3f, t: %.3f, E1: %.6e\n",
     rms(vs.q)/pr.f, maximum(sqrt.(KE)), vs.t/tnd, E1/E10)

  pause(0.01)

  nothing
end

# Physical parameters
Lx  = 2*pi*200e3              # Domain extent
f0  = 1e-4                    # Inertial or Coriolis frequency
nu0 = 5e7                     # Potential vorticity hyperdiffusivity
nu1 = 5e7                     # Wave hyperviscosity
nnu = 4                       # Order of hyperviscosity and hyperdiffusivity
N0  = 5e-3                    # Buoyancy frequency
m   = 2*pi/325                # Vertical wavenumber

# Initial condition
Uw  = 1e-1                    # Wave speed
Ue  = 5e-2                    # Eddy speed
R   = Lx/15                   # Eddy radius
ke  = 2*pi/R                  # Inverse eddy scale
te  = 1/(Ue*ke)               # Eddy turn-over time
tf  = 2*pi/f0                 # Inertial period

# Numerical params
nx  = 256                     # Resolution
dt  = 1e-2 * tf               # Time-step
nsteps = Int(ceil(30*te/dt))  # Total number of time-steps
nsubsteps = ceil(Int, tf/dt)  # Number of steps between plots
nplots = ceil(Int, nsteps/nsubsteps)   # Number of plots


g  = TwoDGrid(nx, Lx)
pr = TwoModeBoussinesq.Params(nu0, nnu, nu1, nnu, f0, N0, m, -Ue, 0.0)
vs = TwoModeBoussinesq.Vars(g)
eq = TwoModeBoussinesq.Equation(pr, g)
ts = ETDRK4TimeStepper(dt, eq.LCc, eq.LCr)

# Initial condition
q0 = FourierFlows.lambdipole(Ue, R, g; center=(0.0, 0.0))
u0 = 0.5*Uw/sqrt(2) * ones(Complex{Float64}, g.nx, g.ny)
v0 = 0.5*Uw/sqrt(2) * ones(Complex{Float64}, g.nx, g.ny)

TwoModeBoussinesq.set_q!(vs, pr, g, q0)
TwoModeBoussinesq.set_uvp!(vs, pr, g, u0, v0, 0.0*u0)

q00 = maximum(q0)

# Plot
fig, axs = subplots(nrows=1, ncols=2, sharex=true, sharey=true,
  figsize=(12, 5))

niwqgplot(axs, vs, pr, g, q00*tf, 2*Uw, R, tf)

for i = 1:nplots

  @time stepforward!(vs, nsubsteps, ts, eq, pr, g)
  niwqgplot(axs, vs, pr, g, q00*tf, 2*Uw, R, tf)

end
