include("../../src/fourierflows.jl")

using FourierFlows, PyPlot
import FourierFlows.TwoModeBoussinesq
import FourierFlows.TwoModeBoussinesq: totalenergy, mode0energy, mode1energy


# Parameters
  n = 256          # Resolution
  L = 2π           # Domain size
nν0 = 4            # Hyperviscosity for 0th mode
nν1 = 4
 ν0 = 1e-2
 ν1 = 1e-2
  f = 1.0
  N = 20
  α = 1.0
nkw = 16

 Re = L/15
 Ro = 0.1
 uw = 1

  σ = f*sqrt(1+α)                                                               
 kw = 2π*nkw/L                                                                 
  m = N*kw/(f*sqrt(α))                                                          
 tσ = 2π/σ
 dt = 1e-2*tσ
ntσ = 100

nsteps = Int(ntσ*tσ/dt)


# Setup
prob = TwoModeBoussinesq.InitialValueProblem(nx=n, Lx=L, ν0=ν0, nν0=nν0, 
  ν1=ν1, nν1=nν1, f=f, N=N, m=m, dt=dt)

x, y = prob.grid.X, prob.grid.Y
Z = f*Ro * exp.(-(x.^2+y.^2)/2Re^2)
TwoModeBoussinesq.set_Z!(prob, Z)

Umax = maximum(sqrt.(prob.vars.U.^2+prob.vars.V.^2))
TwoModeBoussinesq.set_planewave!(prob, uw*Umax, nkw)

etot = Diagnostic(totalenergy, prob; nsteps=nsteps)
  e0 = Diagnostic(mode0energy, prob; nsteps=nsteps)
  e1 = Diagnostic(mode1energy, prob; nsteps=nsteps)
diags = [etot, e0, e1]


# Run and plot
stepforward!(prob, diags; nsteps=nsteps)
TwoModeBoussinesq.updatevars!(prob)

fig, axs = subplots(ncols=2, nrows=1, figsize=(10, 5))

axes(axs[1])
pcolormesh(x, y, real.(prob.vars.u))

axes(axs[2])
plot(etot.time, etot.data/etot.data[1])
plot(e0.time, e0.data/e0.data[1])
plot(e1.time, e1.data/e1.data[1])

show()
