include("../../src/fourierflows.jl")

using PyPlot,
      JLD

import FourierFlows,
       FourierFlows.TwoDTurb

include("./twomodeutils.jl")



# Parameters
name   = "testcompare"
ep     = 0.1
nkw    = 0                                    # Near-inertial wave
f0, N0 = 1e-4, 5e-3                           # Planetary vorticity
m      = 2*pi / 1600

nx   = 1024
Lx   = 2*pi*200e3                           # Domain extent
dt   = 1e-3 * 2*pi/f0
nnu0 = 4
nnu1 = 8
nu0  = 1e-1/(dt*(0.65*pi*nx/Lx)^nnu0)  # Hyperviscosity
nu1  = 1e-1/(dt*(0.65*pi*nx/Lx)^nnu1)  # Hyperviscosity

sig  = sqrt(f0^2.0 + N0^2.0*(nkw*2*pi/Lx)^2.0 / m^2.0)
alph = sig^2.0/f0^2.0 - 1.0
tsig = 2.0*pi / sig

nsteps = 1000
nsubs  = 1 #Int(ceil(2*pi/f0 / dt))




# Two-mode Boussinesq
g  = FourierFlows.TwoDGrid(nx, Lx)
vs = TwoDTurb.Vars(g)
pr = TwoDTurb.Params(nu0, nnu0)
eq = TwoDTurb.Equation(pr, g)
ts = ETDRK4TimeStepper(dt, eq.LC)

# Load initial vorticity field
savename = @sprintf("twodturb_nx%04d_nnu%d.jld", nx, nnu0)
Z0 = load(savename, "Z")

# Set for 2D turbulence
TwoDTurb.set_q!(vs, pr, g, Z0)
maxsp = maximum([vs.U; vs.V])

# Eddy length scale:
ke = (
  FourierFlows.parsint(sqrt.(g.KKrsq).*g.KKrsq.*abs2.(vs.psih), g) / 
    FourierFlows.parsint(g.KKrsq.*abs2.(vs.psih), g)
)

fig, axs = subplots(ncols=2, figsize=(12, 4))
axs[1][:imshow](vs.q)
axs[2][:imshow](sqrt.(vs.U.^2.0 + vs.V.^2.0))
title(@sprintf("Eddy scale 1/kL: %.2f", 1/(ke*Lx)))
show()




# Plotting
message(vs, pr, g)   = @sprintf("\$t = %.1f\$ wave periods", vs.t/tsig)
rossby2d(vs, pr, g)  = vs.q / f0
meanspeed(vs, pr, g) = sqrt.(vs.U.^2.0 + vs.V.^2.0)

R00 = 0.8*maximum(abs.(rossby2d(vs, pr, g)))
S00 = 1.0*maximum(meanspeed(vs, pr, g))

meanplot = FourierFlows.TwoComponentPlot(g, vs, pr, 
# Component            Color limits      Color name
  rossby2d,         [-R00, R00],      "RdBu_r",
  meanspeed,        [0.0,  S00],      "YlGnBu_r",
  message, @sprintf("./plots/%s_2d", name), 240
)

FourierFlows.makeplot!(meanplot, save=true)




# Run
startwalltime = time()
for step = 1:Int(ceil(nsteps/nsubs))

  stepforward!(vs, ts, eq, pr, g; nsteps=nsubs)
  TwoDTurb.updatevars!(vs, pr, g)

  sp = meanspeed(vs, pr, g)

  @printf("
    step: %04d, t: %.3f, wall time: %.3f, CFL: %.3f, max Z/f: %.2e\n\n",
    ts.step, vs.t/tsig, time()-startwalltime,
    maximum([vs.U; vs.V])*ts.dt/g.dx, maximum(vs.q)/f0
  )

  FourierFlows.makeplot!(meanplot, save=true, show=true)

end
