include("/Users/glwagner/Numerics/FourierFlows/src/fourierflows.jl")

using FourierFlows,
      PyPlot,
      JLD

import FourierFlows.TwoModeBoussinesq

include("./twomodeutils.jl")
include("./maketurb.jl")

alph = 1.0
nkw  = 16
nx   = 384                                  # Resolution
ep   = 5e-2                                 # Wave nonlinearity
Ro   = 1e-1                                 # Eddy Rossby number
name = @sprintf("wavesturb_a10%02d", nkw)

# Plot parameters
plotpath = "./plots"
iplot = 0
dpi = 240

# Physical parameters
Lx    = 2*pi*1600e3                         # Domain extent
f0    = 1e-4                                # Inertial or Coriolis frequency
N0    = 5e-3                                # Buoyancy frequency
kw    = 2*pi*nkw/Lx                         # Wavenumber

# Initial condition
sig   = f0*sqrt(1+alph)                     # Wave frequency
m     = N0*kw/(f0*sqrt(alph))               # Vertical scale
uw    = ep*sig/kw                           # Wave velocity
R     = Lx/20                               # Eddy radius
tsig  = 2*pi/sig                            # Inertial period

# Numerical params
dt    = 2.0e-2 * tsig                       # Time-step
nnu   = 8                                   # Hyperviscous order
nu0   = 1e-1/(dt*(0.65*pi*nx/Lx)^nnu)       # Zeroth mode hyperviscosity
nu1   = 1e-1/(dt*(0.65*pi*nx/Lx)^nnu)       # First mode hyperviscosity
nsteps = 200 * ceil(Int, tsig/dt)           # Total number of time-steps
nsubs  = ceil(Int, tsig/dt)                # Number of steps between plots
nplots = ceil(Int, nsteps/nsubs)           # Number of plots

turbname = "turb_384_1.jld"
#Z0 = matureturb(nx, Lx, 5*dt, nu0, nnu, f0; kpeak=384, Ro=Ro, Ro0=1.3*Ro)
#save(turbname, "Z", Z0)
Z0 = load("./turb_384_0.jld", "Z")

fig, axs = subplots()
imshow(Z0)
savefig(@sprintf("%s.png", turbname), facecolor="w", dpi=240)


@printf("
  *** %s *** 
  max Ro: %.2f, alph: %.3f, sig/f: %.1f, ep: %.2f, uw: %.2f m/s, nkw: %d\n\n",
  name, maximum(Z0/f0), alph, sig/f0, ep, uw, nkw
)


# Initialize problem
g  = TwoDGrid(nx, Lx)
pr = TwoModeBoussinesq.Params(nu0, nnu, nu1, nnu, f0, N0, m)
vs = TwoModeBoussinesq.Vars(g)
eq = TwoModeBoussinesq.Equation(pr, g)
ts = ETDRK4TimeStepper(dt, eq.LCc, eq.LCr)

# Make initial condition
TwoModeBoussinesq.set_zeta!(vs, pr, g, Z0)
TwoModeBoussinesq.set_planewave!(vs, pr, g, uw, nkw)

# Plot initialization
fig3, axs3 = subplots(nrows=1, ncols=3, sharex=true, sharey=true,
  figsize=(12, 4))

fig2, axs2 = subplots(nrows=1, ncols=2, sharex=true, sharey=true,
  figsize=(10, 5))

# Potential plot functions
rossbyq(vs, pr, g)      = TwoModeBoussinesq.calc_apv(vs, pr, g) / pr.f
rossbynum(vs, pr, g)    = vs.Z / pr.f
waveu(vs, pr, g)        = real.(vs.u + conj.(vs.u))
wavev(vs, pr, g)        = real.(vs.v + conj.(vs.v))
wavew(vs, pr, g)        = real.(vs.w + conj.(vs.w))
wavespeed(vs, pr, g)    = sqrt.(waveu(vs, pr, g).^2.0 + wavev(vs, pr, g).^2.0)
wavepressure(vs, pr, g) = real.(vs.p + conj.(vs.p)) 
wavebuoyancy(vs, pr, g) = real.(im*pr.m*vs.p - im*pr.m*conj.(vs.p))
meanspeed(vs, pr, g)    = sqrt.(vs.U.^2.0 + vs.V.^2.0)
message(vs, pr, g) = @sprintf("\$t = %.1f\$ wave periods", vs.t/tsig)

function waveinducedflow(vs, pr, g)
  uw, vw = calc_uw(sig, vs, pr, g)
  return sqrt.(uw.^2.0 + vw.^2.0)
end

function waveinducedu(vs, pr, g)
  uw, vw = calc_uw(sig, vs, pr, g)
  return uw
end

function waveinducedv(vs, pr, g)
  uw, vw = calc_uw(sig, vs, pr, g)
  return vw
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


# Some plot properties
d = 1.0 # domain fraction for plot
R00 = 0.8*maximum(abs.(rossbynum(vs, pr, g)))
p00 = 2.0*maximum(wavepressure(vs, pr, g))
b00 = 1.5*maximum(wavebuoyancy(vs, pr, g))
S00 = 1.0*maximum(meanspeed(vs, pr, g))
u00 = 4.0*uw


basicplot = TwoComponentPlot(
  g, vs, pr, 
# Component            Title               Color limits      Color name
  rossbynum,           "",                 [-R00, R00],      "RdBu_r",
  wavespeed,           "",                 [0.0,  u00],      "YlGnBu_r",
  [-0.5*d*g.Lx, 0.5*d*g.Lx], [-0.5*d*g.Ly, 0.5*d*g.Ly], R, L"x/R", L"y/R",
  message, @sprintf("./plots/%s_basic", name), dpi
)

waveinducedplot = ThreeComponentPlot(
  g, vs, pr, 
# Component            Title               Color limits      Color name
  wavespeed,           "",                 [0.0,  u00],      "YlGnBu_r",
  waveinducedu,        "",                 [-S00, S00],      "RdBu_r",
  waveinducedv,        "",                 [-S00, S00],      "RdBu_r",
  [-0.5*d*g.Lx, 0.5*d*g.Lx], [-0.5*d*g.Ly, 0.5*d*g.Ly], R, L"x/R", L"y/R",
  message, @sprintf("./plots/%s_waveflow", name), dpi
)

apvplot = TwoComponentPlot(
  g, vs, pr, 
# Component            Title               Color limits      Color name
  rossbyq,             "",                 [-R00, R00],      "RdBu_r",
  apvinducedflow,      "",                 [0.0, 0.5*S00],   "YlGnBu_r",
  [-0.5*d*g.Lx, 0.5*d*g.Lx], [-0.5*d*g.Ly, 0.5*d*g.Ly], R, L"x/R", L"y/R",
  message, @sprintf("./plots/%s_apvflow", name), dpi
)

flowsplot = ThreeComponentPlot(
  g, vs, pr, 
# Component            Title               Color limits      Color name
  rossbyq,             "",                 [-R00, R00],      "RdBu_r",
  wavespeed,           "",                 [0.0,  u00],      "YlGnBu_r",
  waveinducedflow,     "",                 [0.0, 0.5*S00],   "YlGnBu_r",
  [-0.5*d*g.Lx, 0.5*d*g.Lx], [-0.5*d*g.Ly, 0.5*d*g.Ly], R, L"x/R", L"y/R",
  message, @sprintf("./plots/%s_inducedv", name), dpi
)

figure(fig2[:number])
makeplot!(axs2, basicplot)
saveplot!(basicplot)

figure(fig3[:number])
makeplot!(axs3, flowsplot)
saveplot!(flowsplot)

# Initial energy
E0i, E1i = TwoModeBoussinesq.calc_energies(vs, pr, g)
Ei = E0i + E1i


# Run
startwalltime = time()
for i = 1:nplots

  stepforward!(vs, nsubs, ts, eq, pr, g)

  TwoModeBoussinesq.updatevars!(vs, pr, g)

  q      = TwoModeBoussinesq.calc_apv(vs, pr, g)
  sp     = TwoModeBoussinesq.calc_speed1(vs) 
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

  figure(fig2[:number])
  makeplot!(axs2, basicplot)
  saveplot!(basicplot)

  figure(fig3[:number])
  makeplot!(axs3, flowsplot)
  saveplot!(flowsplot)

end
