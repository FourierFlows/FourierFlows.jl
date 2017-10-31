include("../../src/fourierflows.jl")

using FourierFlows,
      PyPlot,
      JLD

import FourierFlows,
       FourierFlows.TwoDTurb,
       FourierFlows.TwoModeBoussinesq

include("./twomodeutils.jl")



# Parameters
ep     = 0.1
nkw    = 0                                    # Near-inertial wave
name   = "testcompare"
f0     = 1e-4                                 # Planetary vorticity
N0     = 5e-3                                 # Buoyancy frequency
m      = 2*pi / 1600

nplots = 1000
nsubs  = 100

nx     = 2048
Lx     = 2*pi*200e3                           # Domain extent
dt     = 2e-2 * 2*pi/f0
nnu0   = 4
nnu1   = 8
nu0    = 1e-1/(dt*(0.65*pi*nx/Lx)^nnu0)  # Hyperviscosity
nu1    = 1e-1/(dt*(0.65*pi*nx/Lx)^nnu1)  # Hyperviscosity

sig    = sqrt(f0^2.0 + N0^2.0*(nkw*2*pi/Lx)^2.0 / m^2.0)
alph   = sig^2.0/f0^2.0 - 1.0
tsig   = 2*pi / sig






# Grid
g = TwoDGrid(nx, Lx)

# Two-D turbulence
vs1 = TwoDTurb.Vars(g)
pr1 = TwoDTurb.Params(nu0, nnu0)
eq1 = TwoDTurb.Equation(pr1, g)
ts1 = ETDRK4TimeStepper(dt, eq1.LC)

# Two-mode Boussinesq
vs2 = TwoModeBoussinesq.Vars(g)
pr2 = TwoModeBoussinesq.Params(nu0, nnu0, nu1, nnu1, f0, N0, m)
eq2 = TwoModeBoussinesq.Equation(pr2, g)
ts2 = ETDRK4TimeStepper(dt, eq2.LCc, eq2.LCr)





# Load initial vorticity field
savename = @sprintf("twodturb_nx%04d_nnu%d.jld", nx, nnu0)

@printf("Loading initial vorticity... ")
t0 = time()
Z0 = load(savename, "Z")
@printf("... done. (t = %.4f s)\n", time()-t0)

# Set for 2D turbulence
TwoDTurb.set_q!(vs1, pr1, g, Z0)
maxsp = maximum([vs1.U; vs1.V])

# Eddy length scale:
ke = (
  FourierFlows.parsevalsum2(sqrt.(g.KKrsq).*g.KKrsq.*abs2.(vs1.psih), g) / 
    FourierFlows.parsevalsum2(g.KKrsq.*abs2.(vs1.psih), g)
)

# Initial wave velocity:
uw = 2.0 * 0.1*maxsp #2.0*f0*ep/ke
TwoModeBoussinesq.set_zeta!(vs2, pr2, g, Z0)
TwoModeBoussinesq.set_planewave!(vs2, pr2, g, uw, nkw)





# Plotting
message(vs, pr, g)   = @sprintf("\$t = %.1f\$ wave periods", vs.t/tsig)
rossby2d(vs, pr, g)  = vs.q / f0
rossbytwomode(vs, pr, g) = TwoModeBoussinesq.calc_apv(vs, pr, g) / pr.f
meanspeed(vs, pr, g) = sqrt.(vs.U.^2.0 + vs.V.^2.0)
waveu(vs, pr, g)     = real.(vs.u + conj.(vs.u))
wavev(vs, pr, g)     = real.(vs.v + conj.(vs.v))
wavespeed(vs, pr, g) = sqrt.(waveu(vs, pr, g).^2.0 + wavev(vs, pr, g).^2.0)

function waveinducedflow(vs, pr, g)
  uw, vw = calc_uw(sig, vs, pr, g)
  sqrt.(uw.^2.0 + vw.^2.0)
end

R00 = 0.8*maximum(abs.(rossbytwomode(vs2, pr2, g)))
S00 = 1.0*maximum(meanspeed(vs2, pr2, g))
u00 = 4.0*uw

nowavesplot = FourierFlows.OneComponentPlot(g, vs1, pr1, rossby2d, [-R00, R00], 
  "RdBu_r", message, @sprintf("./plots/%s_2d", name), 240
)

wavesplot = FourierFlows.ThreeComponentPlot(g, vs2, pr2, 
# Component            Color limits      Color name
  rossbytwomode,       [-R00, R00],      "RdBu_r",
  wavespeed,           [0.0,  u00],      "YlGnBu_r",
  waveinducedflow,     [0.0, 0.1*S00],   "YlGnBu_r",
  message, @sprintf("./plots/%s_twomode", name), 240
)

FourierFlows.makeplot!(nowavesplot, save=true)
FourierFlows.makeplot!(wavesplot, save=true)

# Initial energy
E0i, E1i = TwoModeBoussinesq.calc_energies(vs2, pr2, g)
Ei = E0i + E1i





# Run
@printf("
  *** %s *** 
  max Ro: %.2f, alph: %.3f, sig/f: %.1f, ep: %.2f, uw: %.2f m/s, nkw: %d\n\n",
  name, maximum(Z0/f0), alph, sig/f0, ep, uw, nkw
)

startwalltime = time()
for i = 1:nplots

  stepforward!(vs1, ts1, eq1, pr1, g; nsteps=nsubs)
  stepforward!(vs2, ts2, eq2, pr2, g; nsteps=nsubs)

  TwoDTurb.updatevars!(vs1, pr1, g)
  TwoModeBoussinesq.updatevars!(vs2, pr2, g)

  q      = TwoModeBoussinesq.calc_apv(vs2, pr2, g)
  sp     = TwoModeBoussinesq.calc_wavespeed(vs2) 
  E0, E1 = TwoModeBoussinesq.calc_energies(vs2, pr2, g)
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

  FourierFlows.makeplot!(nowavesplot, save=true)
  FourierFlows.makeplot!(wavesplot, save=true)

end
