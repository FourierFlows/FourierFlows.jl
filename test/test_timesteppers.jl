using Base.Test, FourierFlows
import FourierFlows.TwoDTurb

steppers = [
  "ForwardEuler",
  "FilteredForwardEuler",
  "AB3",
  "RK4",
  "ETDRK4",
  "FilteredETDRK4",
]

steppersteps = [
  1,
  1,
  3,
  1,
  1,
  1,
]

filteredsteppers = [
  "FilteredForwardEuler",
  "FilteredETDRK4",
]

"""
Build a twodturb problem, and use it to test time-stepping methods.
"""
function testtwodturbstepforward(n=64, L=2π, ν=0.0, nν=2; 
                                 dt=1e-16, nsteps=1, stepper="ForwardEuler")

  g  = TwoDTurb.Grid(n, L)
  p  = TwoDTurb.Params(ν, nν)
  v  = TwoDTurb.Vars(g)
  eq = TwoDTurb.Equation(p, g)
  ts = FourierFlows.autoconstructtimestepper(stepper, dt, sol, g)
  
  # Initial condition (IC)
  qi = rand(g.nx, g.ny)
  qih = rfft(qi)

  # Remove high wavenumber power from IC
  cutoff = 0.3
  highwavenumbers = (g.kr/g.k[2]).^2 .+ (g.l/g.l[2]).^2 .> cutoff # max is 0.5
  qih[highwavenumbers] =  0
  qi = irfft(qih, g.nx)

  prob = Problem(g, v, p, eq, ts)
  TwoDTurb.set_q!(prob, qi)

  absq₀ = abs.(prob.state.sol)
  stepforward!(prob, nsteps)
  absq₁ = abs.(prob.state.sol)

  isapprox(sum(absq₀-absq₁), 0.0, rtol=n^2*nsteps*1e-15)
end


# Run the tests
n = 64 

for (i, stepper) in enumerate(steppers)
  @test testtwodturbstepforward(n; stepper=stepper, nsteps=steppersteps[i])
end
