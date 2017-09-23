include("/Users/glwagner/Numerics/FourierFlows/src/fourierflows.jl")

using FourierFlows,
      PyPlot

import FourierFlows.TwoModeBoussinesq

include("./twomodeutils.jl")




# Physical parameters
Lx   = 2*pi*100e3                          # Domain extent
f0   = 1e-4                                # Inertial or Coriolis frequency
N0   = 5e-3                                # Buoyancy frequency
sig  = 2*f0                                # Wave frequency
nkw  = 16                                  # Non-dimensional wavenumber
kw   = 2*pi*nkw/Lx                         # Wavenumber
m    = N0*kw/sqrt(sig^2-f0^2)              # Vertical scale

# Initial condition
Ro   = 1e-1                                # Eddy Rossby number
ep   = 2e-2                                # Wave nonlinearity
uw   = ep*f0*Lx/(2*pi*nkw)                 # Wave velocity
R    = Lx/20                               # Eddy radius
tsig = 2*pi/f0                             # Inertial period

# Numerical params
nx   = 256                                 # Resolution
dt   = 2e-2 * tsig                         # Time-step
nnu  = 8                                   # Hyperviscous order
nu0  = 1e-1/(dt*(0.65*pi*nx/Lx)^nnu)       # Zeroth mode hyperviscosity
nu1  = 1e-1/(dt*(0.65*pi*nx/Lx)^nnu)       # First mode hyperviscosity
nsteps = 10000                             # Total number of time-steps
nsubs  = 100 #ceil(Int, 4*tf/dt)           # Number of steps between plots
nplots = ceil(Int, nsteps/nsubs)           # Number of plots




# Initialize problem
g  = TwoDGrid(nx, Lx)
pr = TwoModeBoussinesq.Params(nu0, nnu, nu1, nnu, f0, N0, m)
vs = TwoModeBoussinesq.Vars(g)
eq = TwoModeBoussinesq.Equation(pr, g)
ts = ETDRK4TimeStepper(dt, eq.LCc, eq.LCr)


# Make initial condition
Z0 = Ro*f0 * exp.(-(g.X.^2+g.Y.^2)/(2*R^2))

TwoModeBoussinesq.set_zeta!(vs, pr, g, Z0)
TwoModeBoussinesq.set_planewave!(vs, pr, g, uw, nkw)


# Plot initialization
fig, axs = subplots(nrows=2, ncols=2, sharex=true, sharey=true,
  figsize=(6, 6))

# Potential plot functions
rossbynum(vs, pr, g)    = TwoModeBoussinesq.calc_apv(vs, pr, g) / pr.f
waveu(vs, pr, g)        = real.(vs.u + conj.(vs.u))
wavev(vs, pr, g)        = real.(vs.v + conj.(vs.v))
wavew(vs, pr, g)        = real.(vs.w + conj.(vs.w))
wavespeed(vs, pr, g)    = sqrt(waveu(vs, pr, g).^2.0 + wavev(vs, pr, g).^2.0)
wavepressure(vs, pr, g) = real.(vs.p + conj.(vs.p)) 
wavebuoyancy(vs, pr, g) = real.(im*pr.m*vs.p - im*pr.m*conj.(vs.p))
meanspeed(vs, pr, g)    = sqrt.(vs.U.^2.0 + vs.V.^2.0)

function waveinducedflow(vs, pr, g)
  uw, vw = calc_uw(sig, vs, pr, g)
  return sqrt.(uw.^2.0 + vw.^2.0)
end

function apvinducedflow(vs, pr, g)
  q = TwoModeBoussinesq.calc_apv(vs, pr, g)
  psiqh = -g.invKKrsq.*rfft(q)
  uq = irfft(-im*g.Lr.*psiqh, g.nx)
  vq = irfft( im*g.Kr.*psiqh, g.nx)
  return sqrt.(uq.^2.0+vq.^2.0)
end

Sp0 = meanspeed(vs, pr, g)

function apvinducedflow_diff(vs, pr, g)
  q = TwoModeBoussinesq.calc_apv(vs, pr, g)
  psiqh = -g.invKKrsq.*rfft(q)
  uq = irfft(-im*g.Lr.*psiqh, g.nx)
  vq = irfft( im*g.Kr.*psiqh, g.nx)
  return sqrt.(uq.^2.0+vq.^2.0) - Sp0
end




# Save some initial properties
E0i, E1i = TwoModeBoussinesq.calc_energies(vs, pr, g)
Ei = E0i + E1i

# Some plot properties
d = 1.0 # domain fraction for plot
R00 = 1.5*maximum(abs.(rossbynum(vs, pr, g)))
p00 = 2.0*maximum(wavepressure(vs, pr, g))
b00 = 2.0*maximum(wavebuoyancy(vs, pr, g))
S00 = 1.0*maximum(meanspeed(vs, pr, g))




pl = FourComponentPlot(
  g, vs, pr, 
  rossbynum,           L"q/f",             [-R00, R00],      "RdBu_r",
  wavebuoyancy,        L"b",               [-b00, b00],      "RdBu_r",
 #apvinducedflow,      L"|\nabla \psi^q|", [0.0,  S00],      "YlGnBu_r",
 #apvinducedflow_diff, L"|\nabla \psi^q|", [-S00, S00]/10,   "RdBu_r",
  apvinducedflow,      L"|\nabla \psi^q|", [0.0,  S00],      "YlGnBu_r",
  waveinducedflow,     L"|\nabla \psi^w|", [0.0,  1e-2*S00], "YlGnBu_r",
 #wavepressure, L"p",              [-p00, p00], "RdBu_r",
 #wavespeed,    L"\sqrt{u^2+v^2}", [0.0, 2*uw], "YlGnBu_r",
  [-0.5*d*g.Lx, 0.5*d*g.Lx], 
  [-0.5*d*g.Ly, 0.5*d*g.Ly], 
  R, 
  L"x/R", L"y/R"
)

makeplot!(axs, pl)




# Run
for i = 1:nplots

  @time stepforward!(vs, nsubs, ts, eq, pr, g)

  TwoModeBoussinesq.updatevars!(vs, pr, g)

  spw = waveinducedflow(vs, pr, g)
  println("\n max wave induced speed: ", maximum(spw))

  q  = TwoModeBoussinesq.calc_apv(vs, pr, g)
  sp = TwoModeBoussinesq.calc_speed1(vs) 
  E0, E1 = TwoModeBoussinesq.calc_energies(vs, pr, g)
  E = E0 + E1

  @printf("
    step: %04d, t: %.3f, CFL: %.3f,
    max Ro: %.2e, max speed: %.2e, 
    E: %.6f, E0: %.6f, E1: %.6f, E0frac: %.3f, E1frac: %.3f\n\n", 
    ts.r.step, vs.t/tsig, 
    maximum([abs.(2*vs.u); abs2.(2*vs.v); vs.U; vs.V])*ts.r.dt/g.dx, 
    maximum(q)/pr.f, maximum(sp), E/Ei, E0/E0i, E1/E1i, E0/Ei, E1/Ei,
  )

  makeplot!(axs, pl)

end
