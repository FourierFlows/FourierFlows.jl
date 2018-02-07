import FourierFlows.BarotropicQG

# -----------------------------------------------------------------------------
# BAROQG's TEST FUNCTIONS

""" Test that the time-stepper is not doing anything wild by evolving a random
initial condition for dt=1e-16 looking at relative error of the norm. """
function test_baroQG_RossbyWave(stepper, dt, nsteps, g, p, v, eq)

    ts = FourierFlows.autoconstructtimestepper(stepper, dt, eq.LC, g)
    prob = FourierFlows.Problem(g, v, p, eq, ts)

    s, v, p, g = prob.state, prob.vars, prob.params, prob.grid

    # the Rossby wave initial condition
     ampl = 1e-2
    kwave = 3.0*2π/g.Lx
    lwave = 2.0*2π/g.Ly
        ω = -p.beta*kwave/(kwave^2.0 + lwave^2.0)
       ζ0 = ampl*cos.(kwave*g.X).*cos.(lwave*g.Y)
      ζ0h = rfft(ζ0)

    BarotropicQG.set_zeta!(prob, ζ0)

    stepforward!(prob, nsteps)
    BarotropicQG.updatevars!(prob)

    ζ_theory = ampl*cos.(kwave*(g.X-ω/kwave*s.t)).*cos.(lwave*g.Y)
    # println(norm(ζ_theory - v.zeta)/norm(ζ_theory))
    isapprox(ζ_theory, v.zeta, rtol=g.nx*g.ny*nsteps*1e-12)
end



# -----------------------------------------------------------------------------
# Running the tests

nx  = 64
ν  = 0.0
νn = 2

f0 = 1.0
 β = 2.0
Lx = 2π
 μ = 0.0

η(x,y) = zeros(nx, nx)
FU(t)  = 0

g  = BarotropicQG.Grid(nx, Lx)
p  = BarotropicQG.Params(g, f0, β, FU, η, μ, ν, νn)
v  = BarotropicQG.Vars(g)
eq = BarotropicQG.Equation(p, g)

dt, nsteps  = 1e-3, 200
@test test_baroQG_RossbyWave("ETDRK4", dt, nsteps, g, p, v, eq)

dt, nsteps  = 1e-3, 200
@test test_baroQG_RossbyWave("RK4", dt, nsteps, g, p, v, eq)

dt, nsteps  =1e-3, 200
@test test_baroQG_RossbyWave("FilteredRK4", dt, nsteps, g, p, v, eq)

# dt, nsteps  = 1e-4, 2000
# @test test_baroQG_RossbyWave("AB3", dt, nsteps, g, p, v, eq)

dt, nsteps  = 1e-4, 2000
@test test_baroQG_RossbyWave("ForwardEuler", dt, nsteps, g, p, v, eq)

dt, nsteps  = 1e-4, 2000
@test test_baroQG_RossbyWave("FilteredForwardEuler", dt, nsteps, g, p, v, eq)
