include("/Users/glwagner/Numerics/FourierFlows/src/fourierflows.jl")

using FourierFlows,
      PyPlot

import FourierFlows.TwoModeBoussinesq

include("./twomodeutils.jl")

nkw  = 16
nx   = 256                                  # Resolution
Ro   = 2e-1                                 # Eddy Rossby number

# Physical parameters
Lx    = 2*pi*1600e3                         # Domain extent
f0    = 1e-4                                # Inertial or Coriolis frequency
N0    = 5e-3                                # Buoyancy frequency
alph0 = 3                                   # Frequency parameter
nkw0  = 16                                  # Non-dimensional wavenumber
kw0   = 2*pi*nkw0/Lx                        # Wavenumber
m     = N0*kw0/(f0*sqrt(alph0))             # Vertical scale

# Initial condition
kw    = 2*pi*nkw/Lx                         # Non-dimensional wavenumber
alph  = (N0*kw/(f0*m))^2    
sig   = f0*sqrt(1+alph)                     # Wave frequency
R     = Lx/20                               # Eddy radius
tsig  = 2*pi/sig                            # Inertial period

# Numerical params
dt    = 2.0e-2 * tsig                       # Time-step
nnu   = 8                                   # Hyperviscous order
nu0 = nu1 = 1e-1/(dt*(0.65*pi*nx/Lx)^nnu)   # Hyperviscosity

# Initialize problem
g  = TwoDGrid(nx, Lx)
pr = TwoModeBoussinesq.Params(nu0, nnu, nu1, nnu, f0, N0, m)
vs = TwoModeBoussinesq.Vars(g)
eq = TwoModeBoussinesq.Equation(pr, g)
ts = ETDRK4TimeStepper(dt, eq.LCc, eq.LCr)


message(vs, pr, g) = @sprintf("\$t = %.1f\$ wave periods", vs.t/tsig)

function eddywave(uw, nsteps, name)

  vs.t = 0.0
  ep = uw*kw/sig

  @printf("
    *** %s *** 
    Ro: %.2f, alph: %.3f, sig/f: %.1f, ep: %.2f, uw: %.2f m/s, nkw: %d\n\n",
    name, Ro, alph, sig/f0, ep, uw, nkw
  )

  # Make initial condition
  Z0 = Ro*f0 * exp.(-(g.X.^2+g.Y.^2)/(2*R^2))
  TwoModeBoussinesq.set_zeta!(vs, pr, g, Z0)
  TwoModeBoussinesq.set_planewave!(vs, pr, g, uw, nkw)

  # Some plot properties
  R00 = 1.0*maximum(abs.(rossbynum(vs, pr, g)))
  S00 = 0.5*maximum(meanspeed(vs, pr, g))
  u00 = 4.0*uw

  basicplot = FourierFlows.ThreeComponentPlot(
    g, vs, pr, 
  # Component            Color limits      Color name
    rossbyq,             [-R00, R00],      "RdBu_r",
    wavespeed,           [0.0,  u00],      "YlGnBu_r",
    waveinducedflow,     [0.0, 1.0*S00],   "YlGnBu_r",
    message, @sprintf("./plots/%s", name)
  )

  FourierFlows.makeplot!(basicplot, save=true, show=true)

  # Initial energy
  E0i, E1i = TwoModeBoussinesq.calc_energies(vs, pr, g)
  Ei = E0i + E1i

  nsubs  = ceil(Int, tsig/dt)        # Number of steps between plots
  nplots = ceil(Int, nsteps/nsubs)   # Number of plots

  # Run
  startwalltime = time()
  for i = 1:nplots

    stepforward!(vs, nsubs, ts, eq, pr, g)

    TwoModeBoussinesq.updatevars!(vs, pr, g)

    q      = TwoModeBoussinesq.calc_apv(vs, pr, g)
    sp     = TwoModeBoussinesq.calc_wavespeed(vs) 
    E0, E1 = TwoModeBoussinesq.calc_energies(vs, pr, g)
    E      = E0 + E1

    @printf("
      step: %04d, t: %.3f, wall time: %.3f,
      CFL: %.3f, max Z/f: %.2e, max q/f: %.2e, max speed: %.2e, 
      E: %.6f, E0: %.6f, E1: %.6f, E0frac: %.3f, E1frac: %.3f\n\n", 
      ts.r.step, vs.t/tsig, time()-startwalltime,
      maximum([abs.(2*vs.u); abs2.(2*vs.v); vs.U; vs.V])*ts.r.dt/g.dx, 
      maximum(vs.Z)/pr.f, maximum(q)/pr.f, maximum(sp), 
      E/Ei, E0/E0i, E1/E1i, E0/Ei, E1/Ei,
    )

    FourierFlows.makeplot!(basicplot; save=true, show=true)

  end

end


nsteps = 40 * ceil(Int, tsig/dt)
epz    = [1e-2, 5e-2]

for iep = 1:length(epz)

  ep = epz[iep]
  uw = ep*sig/kw
  name = @sprintf("eddy_ep%02d", Int(ep*100))

  eddywave(uw, nsteps, name)

end
