include("../src/fourierflows.jl")

using FourierFlows,
      PyPlot

import FourierFlows.NIWQG

rms(a) = sqrt(mean(a.^2))

function niwqgplot(axs, vs, pr, g, q0, Uw, R, tnd) 

  domfrac = 5
  xl, xr = -g.Lx/domfrac, g.Lx/domfrac
  yl, yr = -g.Ly/domfrac, g.Ly/domfrac


  # All lengths non-dimensionalized by R
  axes(axs[1])
  pcolormesh(g.X/R, g.Y/R, vs.q*tnd, cmap="RdBu_r",
    vmin=-q0, vmax=q0)
  xlim(xl/R, xr/R); ylim(yl/R, yr/R)
  xlabel(L"x/R"); ylabel(L"y/R")
  title(L"q")

  axes(axs[2])
  pcolormesh(g.X/R, g.Y/R, 0.5*abs2.(vs.phi), cmap="YlGnBu_r",
    vmin=0, vmax=0.5*Uw^2)
  xlim(xl/R, xr/R); ylim(yl/R, yr/R)
  xlabel(L"x/R")
  title(L"\sqrt{u^2+v^2}")


  @printf("rms Ro: %.2e, max speed: %.3f, t: %.3f\n",
     rms(vs.q)/pr.f, maximum(abs.(vs.phi)), vs.t/tnd)

  pause(0.01)

  nothing
end


# Physical parameters
Lx  = 2*pi*200e3              # Domain extent
f0  = 1e-4                    # Inertial or Coriolis frequency
kap = 1e8                     # Potential vorticity hyperdiffusivity
nu  = 1e8                     # Wave hyperviscosity
nnu = nkap = 4                # Order of hyperviscosity and hyperdiffusivity
N0, m = 5e-3, 2*pi/325
eta = N0^2.0/(f0*m^2.0)       # Wave dispersivity

# Initial condition
Uw  = 1e-1                    # Wave speed
Ue  = 5e-2                    # Eddy speed
R   = Lx/15                   # Eddy radius
ke  = 2*pi/R                  # Inverse eddy scale
te  = 1/(Ue*ke)               # Eddy turn-over time

# Numerical params
nx  = 256                     # Resolution
dt  = 0.01 * te               # Time-step
nsteps = Int(ceil(30*te/dt))  # Total number of time-steps
nsubsteps = ceil(Int, te/dt)  # Number of steps between plots
nplots = ceil(Int, nsteps/nsubsteps)   # Number of plots

@printf("Rossby number: %.3f", 1/(te*f0))


pr = NIWQG.Params(kap, nkap, nu, nnu, eta, f0, -Ue, 0.0)
g  = TwoDGrid(nx, Lx)
vs = NIWQG.Vars(g)
eq = NIWQG.Equation(pr, g)
ts = ETDRK4TimeStepper(dt, eq.LCc, eq.LCr)


# Initial condition
q0   = FourierFlows.lambdipole(Ue, R, g; center=(0.0, 0.0))
phi0 = (1.0+im)/sqrt(2)*Uw * ones(Complex{Float64}, g.nx, g.ny)

NIWQG.set_q!(vs, pr, g, q0)
NIWQG.set_phi!(vs, pr, g, phi0)

q00 = maximum(q0)

# Plot
fig, axs = subplots(nrows=1, ncols=2, sharex=true, sharey=true,
  figsize=(12, 5))

niwqgplot(axs, vs, pr, g, q00*te, 2*Uw, R, te)

for i = 1:nplots

  @time stepforward!(vs, ts, eq, pr, g, nsteps=nsubsteps)
  NIWQG.updatevars!(vs, pr, g)
  niwqgplot(axs, vs, pr, g, q00*te, 2*Uw, R, te)

end
