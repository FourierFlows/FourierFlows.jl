import FourierFlows.BarotropicQG

# -----------------------------------------------------------------------------
# BAROQG's TEST FUNCTIONS

""" Test that the time-stepper is not doing anything wild by evolving a random
initial condition for dt=1e-16 looking at relative error of the norm. """
function test_baroQG_timestep(grid, params, vars, eq, dt;
            stepper="ForwardEuler")

    filter = ones(v.zetah)
    if stepper == "FilteredForwardEuler" || stepper == "FilteredETDRK4" || stepper == "FilteredRK4"
        # create a filtered ts (simpler is Euler)
        ts = FourierFlows.autoconstructtimestepper("FilteredForwardEuler",
                                                    dt, eq.LC, g)
        # and use its filter to apply it to the initial condition
        filter = ts.filter
    end

    ts = FourierFlows.autoconstructtimestepper(stepper, dt, eq.LC, g)
    prob = FourierFlows.Problem(g, v, p, eq, ts)

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

    BarotropicQG.set_zeta!(prob, ζ0)

    nsteps = 5

    stepforward!(prob, nsteps)
    BarotropicQG.updatevars!(prob)

    isapprox(ζ0, prob.vars.zeta, rtol=nsteps*1e-14)
end



# -----------------------------------------------------------------------------
# Running the tests

nx  = 64
dt  = 1e-16
ν  = 1e-6
νn = 2

f0 = 1.0
 β = 0.5
Lx = 2.0*pi
 μ = 1e-3

η(x,y) = zeros(nx, nx)
FU(t)  = 0

g  = BarotropicQG.Grid(nx, Lx)
p  = BarotropicQG.Params(g, f0, β, FU, η, μ, ν, νn)
v  = BarotropicQG.Vars(g)
eq = BarotropicQG.Equation(p, g)

@test test_baroQG_timestep(g, p, v, eq, dt; stepper = "ForwardEuler")
@test test_baroQG_timestep(g, p, v, eq, dt; stepper = "FilteredForwardEuler")
@test test_baroQG_timestep(g, p, v, eq, dt; stepper = "AB3")
@test test_baroQG_timestep(g, p, v, eq, dt; stepper = "RK4")
@test test_baroQG_timestep(g, p, v, eq, dt; stepper = "FilteredRK4")
@test test_baroQG_timestep(g, p, v, eq, dt; stepper = "ETDRK4")
@test test_baroQG_timestep(g, p, v, eq, dt; stepper = "FilteredETDRK4")
