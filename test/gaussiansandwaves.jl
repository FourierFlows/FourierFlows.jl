include("../src/fourierflows.jl")

using  FourierFlows, FourierFlows.NIWQG, PyPlot
import FourierFlows.NIWQG


# Rossby number and "epsilon" (ep = Uwave * kwave / sigma)
# Use these to control amplitude of Gaussian eddy and plane wave
Ro = 0.05 #0.0 #0.05
ep = 0.1

n, L = 256, 2π*200e3
f, N = 1e-4, 5e-3
dt   = 2e-2 * 2π/f

nsubs  = 25 #round(Int, 2π/(f*dt))
nsteps = 10000

nnu, nkap = 4, 4
nu  = 1e-1/(dt*(0.65π*n/L)^nnu)
kap = 1e-1/(dt*(0.65π*n/L)^nkap)

println("nu: ", nu, " kap: ", kap)

nkw   = 4
kw    = nkw*2π/L
sigma = 1.2f
alpha = sigma^2/f^2 - 1
kapw  = kw/sqrt(alpha)
eta   = kapw^2/f

uw = ep*sigma/kw
u00 = uw * sigma/(f*sqrt(2)*sqrt(2+alpha))
v00 = u00 * f/sigma

Reddy = L/20

prob = NIWQG.InitialValueProblem(n, L, dt, kap, nkap, nu, nnu, eta, f)

x, y = prob.grid.X, prob.grid.Y
q0 = f*Ro * exp.(-(x.^2+y.^2)/(2*Reddy^2))
NIWQG.set_q!(prob, q0)

phi0 = uw*exp.(im*kw*x)
NIWQG.set_phi!(prob, phi0)


# Diagnostics
action   = Diagnostic(NIWQG.niwke,    prob; nsteps=nsteps)
qgenergy = Diagnostic(NIWQG.qgenergy, prob; nsteps=nsteps)
wavepe   = Diagnostic(NIWQG.niwpe,    prob; nsteps=nsteps)
diags = [action, qgenergy, wavepe]

coupledei = qgenergy.data[1] + wavepe.data[1]


# Plotting
fig, axs = subplots(ncols=2, nrows=1, sharex=true, sharey=true,
  figsize=(10, 5))


# Loop
@printf("\nTime-stepping. Ro: %.3f, ep: %.3f\n", Ro, ep)
while prob.step < nsteps
  FourierFlows.stepforward!(prob, diags; nsteps=nsubs)
  NIWQG.updatevars!(prob)
  increment!(diags)

  @printf("step: %04d, action: %.3f, qgenergy: %.3f, wavepe: %.3f\n",
    prob.step, action.value/action.data[1], qgenergy.value/coupledei,
    wavepe.value/coupledei)

  xe, ye = x/Reddy, y/Reddy

  axs[1][:cla]()
  axs[2][:cla]()

  axes(axs[1])
  pcolormesh(xe, ye, prob.vars.q, cmap="RdBu_r")

  axes(axs[2])
  pcolormesh(xe, ye, real.(prob.vars.phi), cmap="RdBu_r",
    vmin=-2*uw, vmax=2*uw)

  Rlim = 6
  axs[1][:set_xlim](-Rlim, Rlim)
  axs[1][:set_ylim](-Rlim, Rlim)
  axs[2][:set_xlim](-Rlim, Rlim)
  axs[2][:set_ylim](-Rlim, Rlim)

  tight_layout()

  pause(0.01)
end
