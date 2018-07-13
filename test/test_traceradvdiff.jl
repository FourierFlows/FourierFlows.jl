import FourierFlows.TracerAdvDiff

# -----------------------------------------------------------------------------
# TRACERADVDIFF's TEST FUNCTIONS

"""
    test_constatvel(; kwargs...)

Evolvesa a Rossby wave and compares with the analytic solution.
"""
function test_constatvel(stepper, dt, nsteps)

  nx  = 128
  Lx  = 2π

  uvel, vvel = 0.2, 0.1
  uin(x, y) = uvel + 0*x
  vin(x, y) = vvel + 0*x

  prob = TracerAdvDiff.ConstDiffSteadyFlowProblem(;
    nx = nx, Lx = Lx, κ = 0.00, u = uin, v = vin, dt = dt, stepper = stepper)

  s, v, p, g = prob.state, prob.vars, prob.params, prob.grid

  σ = 0.1
  c0func(x, y) = 0.1*exp.(-(x.^2+y.^2)/(2σ^2))
  # + 0*ones(size(g.X)) + 0*sin.(1/2*g.X).*sin.(1/2*g.Y)

  c0 = c0func.(g.X, g.Y)
  tfinal = nsteps*dt
  cfinal = c0func.(g.X-uvel*tfinal, g.Y-vvel*tfinal)

  TracerAdvDiff.set_c!(prob, c0)

  stepforward!(prob, nsteps)

  TracerAdvDiff.updatevars!(prob)

  isapprox(cfinal, v.c, rtol=g.nx*g.ny*nsteps*1e-12)
end


# -----------------------------------------------------------------------------
# Running the tests

dt, nsteps  = 1e-2, 100
@test test_constatvel("RK4", dt, nsteps)
