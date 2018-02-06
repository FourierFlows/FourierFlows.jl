import FourierFlows.TwoDTurb


"""
Build a twodturb initial value problem and use it to test time-stepping methods.
We integrate a random IC from 0 to t=tf.
The amplitude of the initial condition is kept low (e.g. multiplied by 1e-5)
so that nonlinear terms do not come into play. This way we can compare the final
state qh(t) with qh(0)*exp(-ν k^2 t).
"""
function testtwodturbstepforward(n=64, L=2π, μ=0.1, nμ=0, ν=0.0, nν=2;
                                 nsteps=100, stepper="ForwardEuler")

  dt = tf/steps

  prob = TwoDTurb.InitialValueProblem(nx=n, Lx=L, ny=n, Ly=L, ν=ν, nν=nν, dt=dt,
                                        stepper=stepper)
  g = prob.grid

  # Initial condition (IC)
  qih = (rand(g.nkr, g.nl) + im*rand(g.nkr, g.nl))*1e-5

  # Remove high wavenumber power from IC
  cutoff = 0.3
  highwavenumbers = sqrt.( (g.Kr*g.dx/π).^2 + (g.Lr*g.dy/π).^2 ) .> cutoff
  qih[highwavenumbers] =  0
  qih[1,1] = 0
  qi = irfft(qih, g.nx)
  qih = rfft(qi)

  q0h = deepcopy(qih)

  TwoDTurb.set_q!(prob, qi)

  stepforward!(prob, nsteps)

  qfh = deepcopy(q0h).*exp.(-ν*g.KKrsq*prob.state.step*dt)

  # println(sum(abs.(qfh-prob.state.sol)) / sum(abs.(qfh)))
  isapprox(qfh, prob.state.sol, rtol=n^2*nsteps*1e-13)
end



# Dictionary of (stepper, nsteps) pairs to test. Each stepper is tested by
# stepping forward nstep times.
steppersteps = Dict([
  ("ForwardEuler", 1000),
  ("FilteredForwardEuler", 1000),
  ("RK4", 100),
  ("FilteredRK4", 100),
  ("ETDRK4", 100),
  ("FilteredETDRK4", 100),
])

# Run the tests
n = 64
tf = 0.1


for (stepper, steps) in steppersteps
  @test testtwodturbstepforward(n; stepper=stepper, nsteps=steps)
end
