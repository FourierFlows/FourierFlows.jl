import FourierFlows.BarotropicQG

# -----------------------------------------------------------------------------
# BAROQG's TEST FUNCTIONS

""" Test that the time-stepper is not doing anything wild by evolving a random
initial condition for dt=1e-16 looking at relative error of the norm. """
function test_baroQG_timestep(grid, params, vars, eq, dt;
            stepper="ForwardEuler")

    filter = ones(v.sol)
    if stepper == "FilteredForwardEuler" || stepper == "FilteredETDRK4"
        # create a filtered ts (simpler is Euler)
        ts = FilteredForwardEulerTimeStepper(dt, g, v;
          filterorder=4.0, innerfilterK=0.65, outerfilterK=1)
        # and use its filter to apply it to the initial condition
        filter = ts.filter
    end

    if stepper == "ForwardEuler"
        ts = ForwardEulerTimeStepper(dt, eq.LC)
    elseif stepper == "FilteredForwardEuler"
        ts = FilteredForwardEulerTimeStepper(dt, g, v;
          filterorder=4.0, innerfilterK=0.65, outerfilterK=1)
    elseif stepper == "AB3"
        ts = AB3TimeStepper(dt, eq.LC)
    elseif stepper == "RK4"
        ts = RK4TimeStepper(dt, eq.LC)
    elseif stepper == "ETDRK4"
        ts = ETDRK4TimeStepper(dt, eq.LC)
    elseif stepper == "FilteredETDRK4"
        ts = FilteredETDRK4TimeStepper(dt, eq.LC, g)
    end


    # the initial conditions (IC)
    ζ0  = randn(g.nx, g.ny)
    ζ0 = ζ0-mean(ζ0)
    ζ0h_copy = rfft(ζ0).*filter
    # the following removes power from IC from wavenumbers larger than 0.3*kmax
    # (this is needed for the filtered time-steppers to give machine precision
    # during the tests)
    ζ0h_copy[sqrt.( (g.Kr*g.dx/π).^2  + (g.Lr*g.dy/π).^2) .> 0.3]=0
    ζ0 = irfft(ζ0h_copy, g.nx)
    ζ0h = rfft(ζ0)


    BarotropicQG.set_zeta!(vars, params, grid, ζ0)

    nsteps = 5

    stepforward!(vars, ts, eq, params, grid, nsteps=nsteps)
    BarotropicQG.updatevars!(vars, params, grid)

    isapprox(ζ0, vars.zeta, rtol=nsteps*1e-14)
end



# -----------------------------------------------------------------------------
# Running the tests

nx  = 128
dt  = 1e-16
nu  = 1e-6
nun = 4

beta = 0.5
Lx = 2.0*pi
mu = 1e-3

g  = BarotropicQG.Grid(nx, Lx)
p  = BarotropicQG.FreeDecayParams(g, beta, 0.0, mu, nu, nun)
v  = BarotropicQG.FreeDecayVars(g)
eq = BarotropicQG.Equation(p, g)


@test test_baroQG_timestep(g, p, v, eq, dt; stepper = "ForwardEuler")
@test test_baroQG_timestep(g, p, v, eq, dt; stepper = "FilteredForwardEuler")
@test test_baroQG_timestep(g, p, v, eq, dt; stepper = "AB3")
@test test_baroQG_timestep(g, p, v, eq, dt; stepper = "RK4")
@test test_baroQG_timestep(g, p, v, eq, dt; stepper = "ETDRK4")
@test test_baroQG_timestep(g, p, v, eq, dt; stepper = "FilteredETDRK4")
