function test_diagnosticsteps(dev::Device=CPU(); nsteps=100, freq=1, ndata=ceil(Int, (nsteps+1)/freq))
  prob = Problem(nx=6, Lx=2π, dev=dev)
  getone(prob) = 1
  d = Diagnostic(getone, prob, nsteps=nsteps, freq=freq, ndata=ndata)
  stepforward!(prob, d, nsteps)
  expectedsteps = cat([0], freq:freq:nsteps, dims=1)
  d.steps[1:d.i] == expectedsteps
end

function test_scalardiagnostics(dev::Device=CPU(); nx=6, Lx=2π, kappa=1e-2, nsteps=100, freq=1, ndata=ceil(Int, (nsteps+1)/freq))
  k1 = 2π/Lx
   τ = 1/(kappa*k1^2) # time-scale for diffusive decay
  dt = 1e-4 * τ # dynamics are resolved.

  prob = Problem(nx=nx, Lx=Lx, kappa=kappa, dt=dt, stepper="ETDRK4", dev=dev)
  g = prob.grid

  ct(t) = exp(-kappa*k1^2*t)

  t1 = dt*nsteps
  c0 = @. sin(k1*g.x)
  c1 = @. ct(t1) * c0 # analytical solution

  set_c!(prob, c0)

  variance(prob) = parsevalsum2(prob.sol, prob.grid) / Lx
  d = Diagnostic(variance, prob, nsteps=nsteps, freq=freq, ndata=ndata)

  stepforward!(prob, d, nsteps)

  da = @. 0.5*ct(d[:t])^2
  isapprox(d[:data], da)
end

function test_basicdiagnostics(dev::Device=CPU(); nx=6, Lx=2π, kappa=1e-2)
  k1 = 2π/Lx
   τ = 1/(kappa*k1^2) # time-scale for diffusive decay
  dt = 1e-2 * τ # dynamics are resolved.

  prob = Problem(nx=nx, Lx=Lx, kappa=kappa, dt=dt, stepper="ETDRK4", dev=dev)
  g = prob.grid

  c0(x) = sin(k1*x)
  set_c!(prob, c0)

  getsol(prob) = prob.sol
  soldiag = Diagnostic(getsol, prob, nsteps=1)

  stepforward!(prob, soldiag, 1)

  isapprox(soldiag.data[soldiag.i], prob.sol)
end
