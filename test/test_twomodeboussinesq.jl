include("../../src/fourierflows.jl")

using FourierFlows,
      PyPlot

import FourierFlows.TwoModeBoussinesq

nx  = 128
dt  = 1e-1
Lx  = 2.0*pi
nu  = 1e-6
nun = 4

nsteps = 2000


g  = TwoDGrid(nx, Lx)
p  = TwoModeBoussinesq.Params(nu, nun)
v  = TwoModeBoussinesq.Vars(g)
eq = TwoModeBoussinesq.Equation(p, g)
ts = ETDRK4TimeStepper(dt, eq.LCr, eq.LCc)

# Initial condition
U0, Lu = 0.001, Lx/20
q0 = rand(nx, nx)
#u0 = U0*exp.(-(g.X.^2+g.Y.^2)/(2*Lu^2))
u0 = U0*ones(Complex{Float64}, nx, nx)
TwoModeBoussinesq.set_q!(v, p, g, q0)
#TwoModeBoussinesq.set_uvp!(v, p, g, u0, 0.0*u0, 0.0*u0)

nplots = 100
nsubsteps = Int(nsteps/nplots)

# Plot
fig, axs = subplots(nrows=1, ncols=2, sharex=true, sharey=true,
  figsize=(12, 5))

for i = 1:nplots

  @time stepforward!(v, nsubsteps, ts, eq, p, g)
  TwoModeBoussinesq.updatevars!(v, p, g)

  @printf "max(u): %.3f, t: %.3f\n" maximum(abs.(v.u)) v.t 

  axes(axs[1])
  pcolormesh(g.x, g.y, v.q)
  axis("tight")
  title(L"q")

  axes(axs[2])
  #pcolormesh(g.x, g.y, sqrt.(real.(v.u).^2 + real.(v.v).^2))
  pcolormesh(g.x, g.y, real.(v.u + conj.(v.u)), vmin=-2*U0, vmax=2*U0, 
    cmap="RdBu_r")
  axis("tight")
  title(L"u")

  #axes(axs[2, 1])
  #pcolormesh(g.x, g.y, real.(v.v))
  #axis("tight")
  #title(L"v")

  #axes(axs[2, 2])
  #pcolormesh(g.x, g.y, real.(v.p))
  #axis("tight")
  #title(L"p")

  pause(0.1)

  #savefig("fourtimesteppers.png"; dpi=240, facecolor="w")

end
