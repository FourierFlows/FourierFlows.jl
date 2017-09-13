include("../src/fourierflows.jl")

using FourierFlows,
      PyPlot

import FourierFlows.NIWQG

function rms(q)
  return sqrt(mean(q.^2))
end


function niwqgplot(axs, vs, pr, g, Uw, R, tnd) 

  domfrac = 20
  xc, yc = mean(g.x), mean(g.y)
  xl, xr = xc - g.Lx/domfrac, xc + g.Lx/domfrac
  yl, yr = yc - g.Ly/domfrac, yc + g.Ly/domfrac

  # All lengths non-dimensionalized by R
  axes(axs[1])
  pcolormesh(g.X/R, g.Y/R, vs.q*tnd, cmap="RdBu_r")
  xlim(xl/R, xr/R); ylim(yl/R, yr/R)
  axis("tight")
  title(L"q")

  axes(axs[2])
  pcolormesh(g.X/R, g.Y/R, abs.(vs.phi), cmap="YlGnBu_r")
  xlim(xl/R, xr/R); ylim(yl/R, yr/R)
  axis("tight")
  title(L"\sqrt{u^2+v^2}")

  @printf("rms Ro: %.2e, max speed: %.3f, t: %.3f\n",
     rms(vs.q)/pr.f, maximum(abs.(vs.phi)), vs.t/tnd)

  nothing
end


# Parameters
Lx  = 2*pi*200e3
f0  = 1e-4
kap = 5e7
nu = 5e8
nnu = nkap = 4

# Dispersivity
N0, m = 5e-3, 2*pi/325
eta = N0^2.0/(f0*m^2.0)

nx  = 512
dt  = 1e-3 * 2*pi/f0
nsteps = 2000

# Initial condition
Uw  = 5e-2          # Wave speed
Ue  = 5e-2          # Eddy speed
R   = 84e3          # Eddy radius
ke  = 2*pi/R        # Inverse eddy scale
te  = R/(2*pi*Ue)   # Eddy turn-over time


g  = TwoDGrid(nx, Lx)
#p  = NIWQG.Params(kap, nkap, nu, nnu, eta, f0)
pr = NIWQG.MeanFlowParams(kap, nkap, nu, nnu, eta, f0, Ue, 0.0)
vs = NIWQG.Vars(g)
eq = NIWQG.Equation(pr, g)
ts = ETDRK4TimeStepper(dt, eq.LCr, eq.LCc)


# Initial condition
q0   = FourierFlows.lambdipole(Ue, R, g)
phi0 = (1.0+im)/sqrt(2)*Uw * ones(Complex{Float64}, g.nx, g.ny)

# Normalize
#Ro = 0.1
#q0   = rand(g.nx, g.ny)
#q0 = Ro*p.f*q0 ./ rms(q0)

NIWQG.set_q!(vs, pr, g, q0)
NIWQG.set_phi!(vs, pr, g, phi0)

# Plot
fig, axs = subplots(nrows=1, ncols=2, sharex=true, sharey=true,
  figsize=(12, 5))

niwqgplot(axs, vs, pr, g, 2*Uw, R, te)
pause(0.01) 
#show()

nplots = 100
nsubsteps = Int(nsteps/nplots)

for i = 1:nplots

  @time stepforward!(vs, nsubsteps, ts, eq, pr, g)

  NIWQG.updatevars!(vs, pr, g)

  niwqgplot(axs, vs, pr, g, 2*Uw, R, te)
  pause(0.1)


end
