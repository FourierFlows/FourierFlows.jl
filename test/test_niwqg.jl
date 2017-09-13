include("../src/fourierflows.jl")

using FourierFlows,
      PyPlot

import FourierFlows.NIWQG

function rms(q)
  return sqrt(mean(q.^2))
end


function niwqgplot(axs, vs, pr, g, Uw, ke, tnd) 

  axes(axs[1])
  pcolormesh(1e-3*g.X, 1e-3*g.Y, vs.q*tnd, cmap="RdBu_r")
  axis("tight")
  title(L"q")

  axes(axs[2])
  pcolormesh(1e-3*g.X, 1e-3*g.Y, abs.(vs.phi),
    cmap="YlGnBu_r")
  axis("tight")
  title(L"\sqrt{u^2+v^2}")

  @printf("rms Ro: %.2e, max speed: %.3f, t: %.3f\n",
     rms(vs.q)/pr.f, maximum(pr.kw^2*abs.(vs.phi)), vs.t/tnd)

  nothing
end


# Parameters
Lx  = 2*pi*200e3
f0  = 1e-4
nu = kap = 5e7
nnu = nkap = 4

# Dispersivity
N0, m = 5e-3, 2*pi/325
eta = 0.5/f0*(N0/m)^2.0
kapw = N0/(m*f0)

nx  = 256
dt  = 1e-1 * 2*pi/f0
nsteps = 20000

# Initial condition
Uw  = 5e-1          # Wave speed
Ue  = 5e-2          # Eddy speed
R   = 84e3          # Eddy radius
te  = R/(2*pi*Ue)   # Eddy turn-over time


g  = TwoDGrid(nx, Lx)
#p  = NIWQG.Params(kap, nkap, nu, nnu, eta, f0)
pr = NIWQG.MeanFlowParams(kap, nkap, nu, nnu, eta, f0, Ue, 0.0)
vs = NIWQG.Vars(g)
eq = NIWQG.Equation(p, g)
ts = ETDRK4TimeStepper(dt, eq.LCr, eq.LCc)


# Initial condition
q0   = FourierFlows.lambdipole(Ue, R, g)
phi0 = Uw/kapw^2 * ones(Complex{Float64}, g.nx, g.ny)

# Normalize
#Ro = 0.1
#q0   = rand(g.nx, g.ny)
#q0 = Ro*p.f*q0 ./ rms(q0)

NIWQG.set_q!(vs, pr, g, q0)
NIWQG.set_phi!(vs, pr, g, phi0)

# Plot
fig, axs = subplots(nrows=1, ncols=2, sharex=true, sharey=true,
  figsize=(12, 5))

niwqgplot(axs, vs, pr, g, 2*Uw, te)
pause(0.01) 
#show()

nplots = 100
nsubsteps = Int(nsteps/nplots)

for i = 1:nplots

  @time stepforward!(vs, nsubsteps, ts, eq, pr, g)

  NIWQG.updatevars!(vs, pr, g)

  niwqgplot(axs, vs, pr, g, 2*Uw, te)
  pause(0.1)


end
