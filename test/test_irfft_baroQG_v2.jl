include("../src/fourierflows.jl")

using FourierFlows,
      FourierFlows.BarotropicQG


module IRFFTbaroQGtests

""" Test whether the irfft's are implemented correctly by evolving a random
initial condition for dt=1e-15 looking at relative error of the norm. """
function test_baroQG_irfft(grid, params, vars, eq)

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
    println(norm(zeta0-v.zeta)/norm(zeta0))

    norm(zeta0-v.zeta)/norm(zeta0) < 1.0e-14
end


end # end IRFFTbaroQGtests


import IRFFTbaroQGtests: test_baroQG_irfft


# Run tests
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

println("test_baroQG_irfft:    ", test_baroQG_irfft(g, p, v, eq))
