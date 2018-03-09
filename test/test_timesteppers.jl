import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy

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


"""
Build a twodturb initial value problem and use it to test time-stepping methods.
We integrate a random IC from 0 to t=tf.
The amplitude of the initial condition is kept low (e.g. multiplied by 1e-5)
so that nonlinear terms do not come into play. This way we can compare the final
state qh(t=tf) with qh(t=0)*exp(-ν k^nν tf).
We choose to linear drag (nν=0) so that we can test the energy since in that
case E(t=tf) = E(t=0)*exp(-2ν tf).
"""
function testtwodturbstepforward(n=64, L=2π, ν=1e-2, nν=0;
                                 nsteps=100, stepper="ForwardEuler")

  dt = tf/nsteps

  prob = TwoDTurb.InitialValueProblem(nx=n, Lx=L, ny=n, Ly=L, ν=ν, nν=nν, dt=dt,
                                        stepper=stepper)
  g, p, v, s = prob.grid, prob.params, prob.vars, prob.state

  E = Diagnostic(energy, prob, nsteps=nsteps)
  diags = [E]

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
  E0 = energy(s, v, g)

  stepforward!(prob, diags, nsteps)

  qfh = deepcopy(q0h).*exp.(-p.ν*g.KKrsq.^p.nν.*s.step*dt)
   Ef = E0.*exp.(-2*p.ν*s.step*dt)

  # println(sum(abs.(qfh-s.sol)) / sum(abs.(qfh)))
  isapprox(qfh, prob.state.sol, rtol=n^2*nsteps*1e-13) &  isapprox(Ef, E[E.count], rtol=n^2*nsteps*1e-13)
end




# Run the tests
n = 64
tf = 0.2


for (stepper, steps) in steppersteps
  @test testtwodturbstepforward(n; stepper=stepper, nsteps=steps)
end
