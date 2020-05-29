function test_diagnosticsteps(dev::Device=CPU(); nsteps=100, freq=1, ndata=ceil(Int, (nsteps+1)/freq))
  prob = Problem(nx=6, Lx=2π, dev=dev)
  getone(prob) = 1
  d = Diagnostic(getone, prob, nsteps=nsteps, freq=freq, ndata=ndata)
  stepforward!(prob, d, nsteps)
  expectedsteps = cat([0], freq:freq:nsteps, dims=1)
  
  return d.steps[1:d.i] == expectedsteps
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
  return isapprox(d[:data], da)
end

function test_basicdiagnostics(dev::Device=CPU(); nx=6, Lx=2π, kappa=1e-2)
  k1 = 2π/Lx
   τ = 1/(kappa*k1^2) # time-scale for diffusive decay
  dt = 1e-2 * τ # dynamics are resolved.

  prob = Problem(nx=nx, Lx=Lx, kappa=kappa, dt=dt, stepper="ETDRK4", dev=dev)
  g = prob.grid
  
  c0 = @. sin(k1*g.x)
  set_c!(prob, c0)

  getsol(prob) = prob.sol
  soldiag = Diagnostic(getsol, prob, nsteps=1)

  stepforward!(prob, soldiag, 1)

  return isapprox(soldiag.data[soldiag.i], prob.sol)
end

function test_extenddiagnostic(dev::Device=CPU(); nx=6, Lx=2π, kappa=1e-2)
  k1 = 2π/Lx
   τ = 1/(kappa*k1^2) # time-scale for diffusive decay
  dt = 1e-2 * τ # dynamics are resolved.

  prob = Problem(nx=nx, Lx=Lx, kappa=kappa, dt=dt, stepper="ETDRK4", dev=dev)
  g = prob.grid
  
  c0 = @. sin(k1*g.x)
  set_c!(prob, c0)
  getsol(prob) = prob.sol
  
  soldiag1 = Diagnostic(getsol, prob, nsteps=1)
  soldiag2 = Diagnostic(getsol, prob, nsteps=4)
  
  nsteps_initially1 = length(soldiag1.t)
  nsteps_extend1 = 12
  FourierFlows.extend!(soldiag1, nsteps_extend1)
  
  nsteps_initially2 = length(soldiag2.t)
  FourierFlows.extend!(soldiag2)
  
  return length(soldiag1.t) == nsteps_initially1 + nsteps_extend1 && length(soldiag2.t) == 2*nsteps_initially2
end

function test_incrementdiagnostic(dev::Device=CPU(); nx=6, Lx=2π, kappa=1e-2)
  k1 = 2π/Lx
   τ = 1/(kappa*k1^2) # time-scale for diffusive decay
  dt = 1e-2 * τ # dynamics are resolved.

  prob = Problem(nx=nx, Lx=Lx, kappa=kappa, dt=dt, stepper="ETDRK4", dev=dev)
  g = prob.grid
  
  c0 = @. sin(k1*g.x)
  set_c!(prob, c0)
  
  dummydiag1(prob) = 2.0
  dummydiag2(prob) = 2.0im
  
  diag1 = Diagnostic(dummydiag1, prob, nsteps=10)
  diag2 = Diagnostic(dummydiag2, prob, nsteps=10)
 
  diags = [diag1, diag2]
  
  FourierFlows.increment!(diags)
  
  return diag1.data[2]==2.0 && diag2.data[2]==2.0im
end

function test_lastindex_getindex(dev::Device=CPU(); nx=6, Lx=2π, kappa=1e-2)
  k1 = 2π/Lx
   τ = 1/(kappa*k1^2) # time-scale for diffusive decay
  dt = 1e-2 * τ # dynamics are resolved.

  prob = Problem(nx=nx, Lx=Lx, kappa=kappa, dt=dt, stepper="ETDRK4", dev=dev)
  g = prob.grid
  
  c0 = @. sin(k1*g.x)
  set_c!(prob, c0)
  
  dummydiag1(prob) = 2.0
  dummydiag2(prob) = 2.0im
  
  diag1 = Diagnostic(dummydiag1, prob, nsteps=10)
  diag2 = Diagnostic(dummydiag2, prob, nsteps=10)
 
  diags = [diag1, diag2]
  
  lastindex_at_beggining = FourierFlows.lastindex.(diags)
  
  ntimes=9
  for j=1:ntimes
    FourierFlows.increment!(diags)
  end
  
  lastindex_at_end = FourierFlows.lastindex.(diags)
  
  return lastindex_at_beggining == 1.0*ones(2) && lastindex_at_end == (ntimes+1)*ones(2) && getindex(diag2, 2:5) == 2.0im*ones(4)
end
