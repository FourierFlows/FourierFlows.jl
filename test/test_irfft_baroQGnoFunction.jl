include("../src/fourierflows.jl")

using FourierFlows,
    #   PyPlot,
      FourierFlows.BarotropicQG

nx  = 128
dt  = 1e-17
nu  = 1e-6
nun = 4

beta = 0.5
Lx = 2.0*pi
mu = 1e-3

g  = Grid(nx, Lx)
p  = FreeDecayParams(g, beta, 0.0, mu, nu, nun)
v  = FreeDecayVars(g)
eq = Equation(p, g)

# a random initial condition
zeta0 = randn(nx, nx)
zeta0 = zeta0 - mean(zeta0)

set_zeta!(v, p, g, zeta0)

ts = ETDRK4TimeStepper(dt, eq.LC)
# ts = ForwardEulerTimeStepper(dt, eq.LC)

# fig, axs = subplots(nrows=1, ncols=2, sharey=true, sharex=true)
#
# function test_plot(g, v, zeta0, fig, axs; iter=1)
#   # Make a plot that compared two-dimensional turbulence solved by
#   # the TwoDTurb and BarotropicQG modules.
#
#   im1 = axs[1][:pcolormesh](g.X, g.Y, v.zeta)
#   im2 = axs[2][:pcolormesh](g.X, g.Y, abs.(v.zeta-zeta0))
#   # axis("square")
#   if iter == 1
#     colorbar(im1, ax=axs[1])
#     colorbar(im2, ax=axs[2])
#   end
#   pause(0.01)
# end


nloops = 1
nsteps = 1

@time stepforward!(v, ts, eq, p, g, nsteps=1)
BarotropicQG.updatevars!(v, p, g)
# test_plot(g, v, zeta0, fig, axs; iter=1)

test = norm(zeta0-v.zeta)/norm(zeta0) < 1.0e-14
println("test_baroQG_irfft:    ", test)
