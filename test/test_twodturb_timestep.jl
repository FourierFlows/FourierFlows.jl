import FourierFlows.TwoDTurb

# -----------------------------------------------------------------------------
# TWODTURB's TEST FUNCTIONS

function makebasicturbproblem(n, L, ν, nν)
  g  = TwoDTurb.TwoDGrid(n, L)
  p  = TwoDTurb.Params(ν, nν)
  v  = TwoDTurb.Vars(g)
  eq = TwoDTurb.Equation(p, g)

  g, p, v, eq
end

function teststepforward(g, p, v, eq; dt=1e-16, nsteps=10,
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
  qi = randn(g.nx, g.ny)
  qi = qi-mean(qi)
  qih_copy = rfft(qi).*filter
  # the following removes power from IC from wavenumbers larger than 0.3*kmax
  # (this is needed for the filtered time-steppers to give machine precision
  # during the tests)
  qih_copy[sqrt.( (g.Kr*g.dx/π).^2  + (g.Lr*g.dy/π).^2) .> 0.3]=0
  qi = irfft(qih_copy, g.nx)
  qih = rfft(qi)

  prob = Problem(g, v, p, eq, ts)
  TwoDTurb.set_q!(prob, qi)

  absq₀ = sum(abs.(prob.vars.sol))

  stepforward!(prob; nsteps=nsteps)

  absq₁ = sum(abs.(prob.vars.sol))

  isapprox(absq₀, absq₁, rtol=nsteps*1e-14)
end

function teststepforward(n::Int, L, ν, nν::Int; stepper="ForwardEuler")
  g, p, v, eq = makebasicturbproblem(n, L, ν, nν)
  teststepforward(g, p, v, eq; stepper=stepper)
end



# -----------------------------------------------------------------------------
# Running the tests

@test teststepforward(128, 2π, 1e-2, 2; stepper="ForwardEuler")
@test teststepforward(128, 2π, 1e-2, 2; stepper="FilteredForwardEuler")
@test teststepforward(128, 2π, 1e-2, 2; stepper="AB3")
@test teststepforward(128, 2π, 1e-2, 2; stepper="RK4")
@test teststepforward(128, 2π, 1e-2, 2; stepper="ETDRK4")
@test teststepforward(128, 2π, 1e-2, 2; stepper="FilteredETDRK4")
