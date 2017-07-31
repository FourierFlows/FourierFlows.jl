include("../../src/fourierflows.jl")

using FourierFlows,
      PyPlot,
      FourierFlows.BarotropicQG

nx  = 128
dt  = 1e-1
nu  = 0e-6
nun = 4

beta = 0.2
Lx = 2.0*pi
mu = 0.0

g  = Grid(nx, Lx)
p  = FreeDecayParams(g, beta, 0.0, mu, nu, nun)
v  = FreeDecayVars(g)
eq = Equation(p, g)

# Rossby wave initial condition

zeta0 = 1e-2
kwave = 3.0
lwave = 2.0
omega = -beta*kwave/(kwave^2.0 + lwave^2.0)

set_zeta!(v, p, g, zeta0*cos.(kwave*g.X).*cos.(lwave*g.Y))

ts = ETDRK4TimeStepper(dt, eq.LC)
# ts = ForwardEulerTimeStepper(dt, eq.LC)

fig, axs = subplots(nrows=1, ncols=2, sharey=true, sharex=true)

function test_plot(g, v, zetaan, fig, axs; iter=1)
  # Make a plot that compared two-dimensional turbulence solved by
  # the TwoDTurb and BarotropicQG modules.

  im1 = axs[1][:pcolormesh](g.X, g.Y, v.zeta)
  im2 = axs[2][:pcolormesh](g.X, g.Y, abs.(v.zeta-zetaan))
  # axis("square")
  if iter == 1
    colorbar(im1, ax=axs[1])
    colorbar(im2, ax=axs[2])
  end
  pause(0.01)
end



nloops = 100
nsteps = 200

#test_plot(g, v, zetaan, fig, axs)

for i = 1:nloops
  @time stepforward!(v, ts, eq, p, g, nsteps)
  BarotropicQG.updatevars!(v, p, g)
  zetaan = zeta0*cos.(kwave*(g.X-omega/kwave*v.t)).*cos.(lwave*g.Y)
  test_plot(g, v, zetaan, fig, axs; iter=i)
end
