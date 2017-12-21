include("../src/FourierFlows.jl")

using FourierFlows,
      FourierFlows.BarotropicQG

module baroQGtests

using FourierFlows,
      FourierFlows.BarotropicQG

""" Test that the time-stepper is not doing anything wild by evolving a random
initial condition for dt=1e-16 looking at relative error of the norm. """
function test_baroQG_timestep(grid, params, vars, eq, dt)

    # a random initial condition
    zeta0 = randn(grid.nx, grid.nx)
    zeta0 = zeta0 - mean(zeta0)

    set_zeta!(vars, params, grid, zeta0)

    ts = ETDRK4TimeStepper(dt, eq.LC)
    # ts = ForwardEulerTimeStepper(dt, eq.LC)

    nloops = 1
    nsteps = 1

    @time stepforward!(vars, ts, eq, params, grid, nsteps=1)
    BarotropicQG.updatevars!(vars, params, grid)
    println("error = ", norm(zeta0-vars.zeta)/norm(zeta0))

    norm(zeta0-vars.zeta)/norm(zeta0) < 1.0e-12
end


end # end IRFFTbaroQGtests


import baroQGtests: test_baroQG_timestep

# Run tests
nx  = 128
dt  = 1e-16
nu  = 1e-6
nun = 4

beta = 0.5
Lx = 2.0*pi
mu = 1e-3

g  = Grid(nx, Lx)
p  = FreeDecayParams(g, beta, 0.0, mu, nu, nun)
v  = FreeDecayVars(g)
eq = Equation(p, g)

println("test_baroQG_irfft:    ", test_baroQG_timestep(g, p, v, eq, dt))
