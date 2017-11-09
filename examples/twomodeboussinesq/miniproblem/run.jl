include("../../../src/fourierflows.jl")

using FourierFlows, PyPlot
import FourierFlows.TwoModeBoussinesq
import FourierFlows.TwoModeBoussinesq: totalenergy, mode0energy, mode1energy


# -- Parameters --
   nkw = 8    
     n = 256
     L = 2π*100e3*nkw
     f = 1e-4
     N = 5e-3
     α = 1.0             # Frequency parameter
     ε = 5e-2            # Wave amplitude
    Ro = 5e-2            # Eddy Rossby number
 Reddy = L/10            # Eddy radius
   nν0 = 4
   nν1 = 4
    ν0 = 1e12
    ν1 = 1e6
dtfrac = 1e-2
  name = "mini"
nsteps = 1000
 nsubs = 100

σ = f*sqrt(1+α)
twave = 2π/σ
kw = 2π*nkw/L
m = N*kw/(f*sqrt(α))
uw = minimum([ε*Reddy*σ, ε*σ/kw])
dt = dtfrac * twave

# Setup
#prob = TwoModeBoussinesq.InitialValueProblem(
prob = TwoModeBoussinesq.PrognosticAPVInitialValueProblem(
  nx=n, Lx=L, ν0=ν0, nν0=nν0, ν1=ν1, nν1=nν1, f=f, N=N, m=m, dt=dt)

# Initial condition
x, y = prob.grid.X, prob.grid.Y
Q0 = f*Ro * exp.(-(x.^2+y.^2)/(2*Reddy^2))
TwoModeBoussinesq.set_Q!(prob, Q0)

TwoModeBoussinesq.set_planewave!(prob, uw, nkw)

etot = Diagnostic(totalenergy, prob; nsteps=nsteps)
e0   = Diagnostic(mode0energy, prob; nsteps=nsteps)
e1   = Diagnostic(mode1energy, prob; nsteps=nsteps) 
diags = [etot, e0, e1]

fig, axs = subplots(ncols=2)

# Run
startwalltime = time()
while prob.step < nsteps

  stepforward!(prob, diags; nsteps=nsubs)
  TwoModeBoussinesq.updatevars!(prob)

  walltime = (time()-startwalltime)/60

  log1 = @sprintf("step: %04d, t: %d, ", prob.step, prob.t/twave)
  log2 = @sprintf("ΔE: %.3f, Δe: %.3f, Δ(E+e): %.6f, τ: %.2f min",
    e0.value/e0.data[1], e1.value/e1.data[1], etot.value/etot.data[1],
    walltime)

  println(log1*log2)

  axes(axs[1]); imshow(prob.vars.Q)
  axes(axs[2]); imshow(real.(prob.vars.u))
  savefig(@sprintf("mini_%04d.png", prob.step), dpi=240)

end
