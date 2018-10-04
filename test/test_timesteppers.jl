using 
  FourierFlows.TracerAdvDiff,
  FFTW

using LinearAlgebra: mul!, ldiv!

# Dictionary of (stepper, nsteps) pairs to test. Each stepper is tested by
# stepping forward nstep times.
steppersteps = Dict([
  ("ForwardEuler", 2000),
  ("FilteredForwardEuler", 2000),
  ("AB3", 400),
  ("FilteredAB3", 400),
  ("RK4", 40),
  ("FilteredRK4", 40),
  ("ETDRK4", 40),
  ("FilteredETDRK4", 40)
])

dualsteppersteps = Dict([
  ("DualRK4", 40),
  ("DualFilteredRK4", 40),
  ("DualETDRK4", 40),
  ("DualFilteredETDRK4", 40),
])

function testtraceradvdifftestforward(n=64, L=2π, kap=1e-2; nsteps=100, stepper="ForwardEuler")
  Lc = L/20
  t0 = Lc^2/2kap
  t1 = 2t0

  csol(x, t) = @. exp(-x^2/(4*kap*t)) / sqrt(4π*kap*t)

  dt = t0/nsteps
  prob = TracerAdvDiff.ConstDiffProblem(nx=n, Lx=L, kap=kap, dt=dt, stepper=stepper, steadyflow=true)
  g, p, v, s = unpack(prob)

  c0 = csol(g.X, t0)
  c1 = csol(g.X, t1) # analytical solution

  set_c!(prob, c0)
  stepforward!(prob, nsteps)
  updatevars!(prob)

  isapprox(c1, prob.vars.c, rtol=n^2*nsteps*1e-13)
end

for (stepper, steps) in steppersteps
  @test testtraceradvdiffstepforward(n; stepper=stepper, nsteps=steps)
end

#=
"""
    testtwodturbstepforward(n, L, nu, nnu; kwargs...)

Build a twodturb initial value problem and use it to test time-stepping methods.
We integrate a random IC from 0 to t=tf. The amplitude of the initial condition
is kept low (e.g. multiplied by 1e-5) so that nonlinear terms do not come into play.
This way we can compare the final state qh(t=tf) with qh(t=0)*exp(-nu k^nnu tf).
We choose linear drag (nnu=0) so that we can test the energy since in that
case E(t=tf) = E(t=0)*exp(-2nu tf).
"""
function testtwodturbstepforward(n=64, L=2π, nu=1e-2, nnu=0; nsteps=100, stepper="ForwardEuler")
  dt = tf/nsteps
  prob = TwoDTurb.InitialValueProblem(nx=n, Lx=L, ny=n, Ly=L, nu=nu, nnu=nnu, dt=dt, stepper=stepper)
  g, p, v, s = prob.grid, prob.params, prob.vars, prob.state
  E = Diagnostic(energy, prob, nsteps=nsteps)
  diags = [E]

  qih = 1e-5*(rand(g.nkr, g.nl) + im*rand(g.nkr, g.nl)) # initial condition
  qih[1, 1] = 0

  # Remove high wavenumber power from IC
  cutoff = 0.3
  highwavenumbers = @. sqrt((g.Kr*g.dx/π)^2 + (g.Lr*g.dy/π)^2 ) > cutoff
  @. qih[highwavenumbers] = 0

  qi = irfft(qih, g.nx)
  qih = rfft(qi)
  q0h = deepcopy(qih)

  TwoDTurb.set_q!(prob, qi)
  E0 = energy(s, v, g)
  stepforward!(prob, diags, nsteps)

  qfh = @. q0h*exp(-p.nu*g.KKrsq^p.nnu*s.step*dt)
  Ef = @. E0*exp(-2*p.nu*s.step*dt)
  isapprox(qfh, prob.state.sol, rtol=n^2*nsteps*1e-13) &  isapprox(Ef, E[E.count], rtol=n^2*nsteps*1e-13)
end

"""
    testverticallyfourierstepforward(n, L, nu, nnu; kwargs...)

Build a Vertically Fourier Boussinesq initial value problem and use it to test
time-stepping methods.
"""
function testverticallyfourierstepforward(n=64, L=2π, nu=1e-2, nnu=0; nsteps=100, stepper="DualRK4")
  dt = tf/nsteps
  prob = VerticallyFourierBoussinesq.Problem(
    nx=n, Lx=L, ny=n, Ly=L, nu0=nu, nnu0=nnu, nu1=nu, nnu1=nnu, dt=dt, stepper=stepper)
  g, p, v, s = prob.grid, prob.params, prob.vars, prob.state

  Zih = 1e-5*(rand(g.nkr, g.nl) + im*rand(g.nkr, g.nl)) # initial condition
  Zih[1, 1] = 0

  uih = 1e-5*(rand(g.nk, g.nl) + im*rand(g.nk, g.nl)) # initial condition
  vih = 1e-5*(rand(g.nk, g.nl) + im*rand(g.nk, g.nl)) # initial condition
  pih = 1e-5*(rand(g.nk, g.nl) + im*rand(g.nk, g.nl)) # initial condition
  uih[1, 1] = 0
  vih[1, 1] = 0
  pih[1, 1] = 0

  # Remove high wavenumber power from IC
  cutoff = 0.3
  highwavenumbersr = @. sqrt((g.Kr*g.dx/π)^2 + (g.Lr*g.dy/π)^2 ) > cutoff
  highwavenumbersc = @. sqrt((g.K*g.dx/π)^2 + (g.L*g.dy/π)^2 ) > cutoff
  Zih[highwavenumbersr] .= 0
  uih[highwavenumbersc] .= 0
  vih[highwavenumbersc] .= 0
  pih[highwavenumbersc] .= 0

  Zi = irfft(Zih, g.nx)
  Zih = rfft(Zi)
  Z0h = deepcopy(Zih)

  ui = ifft(uih)
  uih = fft(ui)
  u0h = deepcopy(uih)

  vi = ifft(vih)
  vih = fft(vi)
  v0h = deepcopy(vih)

  pi = ifft(vih)
  pih = fft(vi)
  p0h = deepcopy(vih)

  VerticallyFourierBoussinesq.set_Z!(prob, Zi)
  stepforward!(prob, nsteps)

  Zfh = @. Z0h*exp(-p.nu0*g.KKrsq^p.nnu0*s.step*dt)
  isapprox(Zfh, prob.state.solr, rtol=n^2*nsteps*1e-13)
end

# Run the tests
n = 64
tf = 0.2

for (stepper, steps) in steppersteps
  @test testtwodturbstepforward(n; stepper=stepper, nsteps=steps)
end

for (stepper, steps) in dualsteppersteps
  @test testverticallyfourierstepforward(n; stepper=stepper, nsteps=steps)
end
=#
