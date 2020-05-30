function test_diagnosticsteps(dev::Device=CPU(); nsteps=100, freq=1, ndata=ceil(Int, (nsteps+1)/freq))
  prob = Problem(nx=6, Lx=2π, dev=dev)
  getone(prob) = 1
  diagnostic = Diagnostic(getone, prob, nsteps=nsteps, freq=freq, ndata=ndata)
  stepforward!(prob, diagnostic, nsteps)
  expectedsteps = cat([0], freq:freq:nsteps, dims=1)
  
  return diagnostic.steps[1:diagnostic.i] == expectedsteps
end

function test_scalardiagnostics(dev::Device=CPU(); nx=6, Lx=2π, kappa=1e-2, nsteps=100, freq=1, ndata=ceil(Int, (nsteps+1)/freq))
  k1 = 2π/Lx
  dt = 1e-2 / (kappa*k1^2) # time-scale for diffusive decay dynamics are resolved
  prob = Problem(nx=nx, Lx=Lx, kappa=kappa, dt=dt, stepper="RK4", dev=dev)
  c0 = @. sin(k1 * prob.grid.x)

  ct(t) = exp(-kappa*k1^2*t)

  t1 = dt*nsteps
  c0 = @. sin(k1*prob.grid.x)
  c1 = @. ct(t1) * c0 # analytical solution

  set_c!(prob, c0)

  variance(prob) = parsevalsum2(prob.sol, prob.grid) / Lx
  diagnostic = Diagnostic(variance, prob, nsteps=nsteps, freq=freq, ndata=ndata)

  stepforward!(prob, diagnostic, nsteps)

  da = @. 0.5*ct(diagnostic[:t])^2
  return isapprox(diagnostic[:data], da)
end

function test_basicdiagnostics(dev::Device=CPU(); nx=6, Lx=2π, kappa=1e-2)
  k1 = 2π/Lx
  dt = 1e-2 / (kappa*k1^2) # time-scale for diffusive decay dynamics are resolved
  prob = Problem(nx=nx, Lx=Lx, kappa=kappa, dt=dt, stepper="RK4", dev=dev)
  c0 = @. cos(k1 * prob.grid.x)
  set_c!(prob, c0)

  getsol(prob) = prob.sol
  soldiag = Diagnostic(getsol, prob, nsteps=1)

  stepforward!(prob, soldiag, 1)

  return isapprox(soldiag.data[soldiag.i], prob.sol)
end

function test_extenddiagnostic(dev::Device=CPU(); nx=6, Lx=2π, kappa=1e-2)
  k1 = 2π/Lx
  dt = 1e-2 / (kappa*k1^2) # time-scale for diffusive decay dynamics are resolved
  prob = Problem(nx=nx, Lx=Lx, kappa=kappa, dt=dt, stepper="RK4", dev=dev)
  c0 = @. sin(k1 * prob.grid.x)
  set_c!(prob, c0)

  getsol(prob) = prob.sol
  
  soldiagnostic1 = Diagnostic(getsol, prob, nsteps=1)
  soldiagnostic2 = Diagnostic(getsol, prob, nsteps=4)
  
  nsteps_initially1 = length(soldiagnostic1.t)
  nsteps_extend1 = 12
  FourierFlows.extend!(soldiagnostic1, nsteps_extend1)
  
  nsteps_initially2 = length(soldiagnostic2.t)
  FourierFlows.extend!(soldiagnostic2)
  
  return length(soldiagnostic1.t) == nsteps_initially1 + nsteps_extend1 && length(soldiagnostic2.t) == 2*nsteps_initially2
end

function test_incrementdiagnostic(dev::Device=CPU(); nx=6, Lx=2π, kappa=1e-2)
  k1 = 2π/Lx
  dt = 1e-2 / (kappa*k1^2) # time-scale for diffusive decay dynamics are resolved
  prob = Problem(nx=nx, Lx=Lx, kappa=kappa, dt=dt, stepper="ETDRK4", dev=dev)
  c0 = @. sin(k1 * prob.grid.x)
  set_c!(prob, c0)
  
  dummydiagnostic1(prob) = 2.0
  dummydiagnostic2(prob) = 2.0im
  
  diagnostic1 = Diagnostic(dummydiagnostic1, prob, nsteps=10)
  diagnostic2 = Diagnostic(dummydiagnostic2, prob, nsteps=10)
 
  diagnostics = [diagnostic1, diagnostic2]
  
  FourierFlows.increment!(diagnostics)
  
  return diagnostic1.data[2]==2.0 && diagnostic2.data[2]==2.0im
end

function test_lastindex_getindex(dev::Device=CPU(); nx=6, Lx=2π, kappa=1e-2)
  k1 = 2π/Lx
  dt = 1e-2 / (kappa*k1^2) # time-scale for diffusive decay dynamics are resolved
  prob = Problem(nx=nx, Lx=Lx, kappa=kappa, dt=dt, stepper="ETDRK4", dev=dev)
  c0 = @. sin(k1 * prob.grid.x)
  set_c!(prob, c0)
  
  dummydiagnostic1(prob) = 2.0
  dummydiagnostic2(prob) = 2.0im
  
  nsteps = 10
  
  diagnostic1 = Diagnostic(dummydiagnostic1, prob, nsteps=nsteps)
  diagnostic2 = Diagnostic(dummydiagnostic2, prob, nsteps=nsteps)
 
  diagnostics = [diagnostic1, diagnostic2]
  
  lastindex_at_beggining = FourierFlows.lastindex.(diagnostics)
  
  ntimes=7
  for j=1:ntimes
    FourierFlows.increment!(diagnostics)
  end
  
  lastindex_at_end = FourierFlows.lastindex.(diagnostics)
  
  return lastindex_at_beggining == 1.0*ones(2) && lastindex_at_end == (ntimes+1)*ones(2) && getindex(diagnostic2, 2:5) == 2.0im*ones(4) && lastindex(diagnostic1) == ntimes+1 && lastindex(diagnostic2) == ntimes+1
end
