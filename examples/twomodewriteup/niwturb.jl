include("../../src/fourierflows.jl")

using PyPlot,
      JLD

import FourierFlows,
       FourierFlows.TwoModeBoussinesq

include("./twomodeutils.jl")



# Parameters
ep     = 0.02
nkw    = 0                                    # Near-inertial wave
f0, N0 = 1e-4, 5e-3                           # Planetary vorticity
m0     = 2*pi/2400
m      = 3*m0
name   = @sprintf("niwturb_ep%03d_m%.3f", Int(ep*1000), Int(m/m0))

nx   = 1024
Lx   = 2*pi*200e3                           # Domain extent

sig  = sqrt(f0^2.0 + N0^2.0*(nkw*2*pi/Lx)^2.0 / m^2.0)
alph = sig^2.0/f0^2.0 - 1.0
tsig = 2.0*pi / sig

dt   = 1e-2 * 2*pi/sig
nnu0 = 4
nnu1 = 4
nu0  = 1e-1/(dt*(0.65*pi*nx/Lx)^nnu0)  # Hyperviscosity
nu1  = 1e-1/(dt*(0.65*pi*nx/Lx)^nnu1)  # Hyperviscosity

nsteps = 10000
nsubs  = Int(ceil(tsig/dt))




# Two-mode Boussinesq
g  = FourierFlows.TwoDGrid(nx, Lx)
vs = TwoModeBoussinesq.Vars(g)
pr = TwoModeBoussinesq.Params(nu0, nnu0, nu1, nnu1, f0, N0, m)
eq = TwoModeBoussinesq.Equation(pr, g)
ts = ETDRK4TimeStepper(dt, eq.LCc, eq.LCr)

# Load initial vorticity field
savename = @sprintf("twodturb_nx%04d_nnu%d.jld", nx, nnu0)
Z0 = load(savename, "Z")

# Set for 2D turbulence
TwoModeBoussinesq.set_zeta!(vs, pr, g, Z0)
maxsp = maximum([vs.U; vs.V])
maxRo = maximum(vs.Z/pr.f)

# Eddy length scale:
ke = (
  FourierFlows.parsevalsum2(sqrt.(g.KKrsq).*g.KKrsq.*abs2.(vs.psih), g) / 
    FourierFlows.parsevalsum2(g.KKrsq.*abs2.(vs.psih), g)
)


# Initial wave velocity:
uw = 2.0 * ep/maxRo*maxsp #2.0*f0*ep/ke
u = uw*ones(g.nx, g.ny)
TwoModeBoussinesq.set_uvp!(vs, pr, g, u, 0.0*u, 0.0*u)
#TwoModeBoussinesq.set_planewave!(vs, pr, g, uw, nkw)

E0i, E1i = TwoModeBoussinesq.calc_energies(vs, pr, g)
Ei = E0i + E1i




# Plotting
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

simplemessage(vs, pr, g)   = @sprintf("\$t = %.1f\$ wave periods", vs.t/tsig)
function energymessage(vs, pr, g)

  # Calc energy
  E0, E1 = TwoModeBoussinesq.calc_energies(vs, pr, g)
  E = E0 + E1

  msg = (
    @sprintf("\$t = %04d\$ wave periods, max Z/f: %.2f, ",
      round(Int, vs.t/tsig), maximum(vs.Z/pr.f))* 
    @sprintf("total E: %4.2f, E0: %4.2f, E1: %4.2f",
      E/Ei, E0/E0i, E1/E1i)
  )

  return msg
end
  
R00 = 0.8*maximum(abs.(rossbytwomode(vs, pr, g)))
S00 = 1.0*maximum(meanspeed(vs, pr, g))
u00 = 4.0*uw

wavesplot = FourierFlows.ThreeComponentPlot(g, vs, pr, 
# Component            Color limits      Color name
  rossbytwomode,       [-R00, R00],      "RdBu_r",
  wavespeed,           [0.0,  u00],      "YlGnBu_r",
  waveinducedflow,     [0.0, 1.0*S00],   "YlGnBu_r",
  energymessage, @sprintf("./plots/%s_twomode", name), 240
)

FourierFlows.makeplot!(wavesplot, save=true)






# Run
@printf("
  *** %s *** 
  max Ro: %.2f, alph: %.3f, sig/f: %.1f, ep: %.2f, uw: %.2f m/s, nkw: %d\n",
  name, maximum(Z0/f0), alph, sig/f0, ep, uw, nkw
)

startwalltime = time()
for step = 1:Int(ceil(nsteps/nsubs))

  stepforward!(vs, ts, eq, pr, g; nsteps=nsubs)
  TwoModeBoussinesq.updatevars!(vs, pr, g)

  q  = TwoModeBoussinesq.calc_apv(vs, pr, g)
  sp = TwoModeBoussinesq.calc_wavespeed(vs) 

  E0, E1 = TwoModeBoussinesq.calc_energies(vs, pr, g)
  E = E0 + E1

  @printf("
    step: %04d, t: %.3f, wall time: %.3f min, CFL: %.3f, max Z/f: %.2e, 
    E: %.6f, E0: %.6f, E1: %.6f, E0frac: %.3f, E1frac: %.3f\n\n", 
    ts.step, vs.t/tsig, (time()-startwalltime)/60,
    maximum([abs.(2*vs.u); abs2.(2*vs.v); vs.U; vs.V])*ts.dt/g.dx, 
    maximum(vs.Z)/pr.f, E/Ei, E0/E0i, E1/E1i, E0/Ei, E1/Ei,
  )

  FourierFlows.makeplot!(wavesplot, save=true, show=true)

end
