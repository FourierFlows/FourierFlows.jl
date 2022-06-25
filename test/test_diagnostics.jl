function test_diagnosticsteps(dev::Device=CPU(); nsteps=100, freq=1, ndata=ceil(Int, (nsteps+1)/freq))
  prob = Problem(dev; nx=6, Lx=2π)
  getone(prob) = 1
  diagnostic = Diagnostic(getone, prob, nsteps=nsteps, freq=freq, ndata=ndata)
  stepforward!(prob, diagnostic, nsteps)
  expectedsteps = cat([0], freq:freq:nsteps, dims=1)
  
  return diagnostic.steps[1:diagnostic.i] == expectedsteps
end

function test_scalardiagnostics(dev::Device=CPU(); nx=6, Lx=2π, κ=1e-2, nsteps=100, freq=1, ndata=ceil(Int, (nsteps+1)/freq))
  k₀ = 2π / Lx
  dt = 1e-2 / (κ * k₀^2) # time-scale for diffusive decay dynamics are resolved
  prob = Problem(dev; nx, Lx, κ, dt, stepper="RK4")
  c0 = @. sin(k₀ * prob.grid.x)

  ct(t) = exp(-κ * k₀^2 * t)

  c₀ = @. sin(k₀ * prob.grid.x)

  set_c!(prob, c₀)

  variance(prob) = parsevalsum2(prob.sol, prob.grid) / Lx
  diagnostic = Diagnostic(variance, prob, nsteps=nsteps, freq=freq, ndata=ndata)

  stepforward!(prob, diagnostic, nsteps)

  da = @. 0.5 * ct(diagnostic[:t])^2

  return isapprox(diagnostic[:data], da)
end

function test_basicdiagnostics(dev::Device=CPU(); nx=6, Lx=2π, κ=1e-2)
  k₀ = 2π / Lx
  dt = 1e-2 / (κ * k₀^2) # time-scale for diffusive decay dynamics are resolved
  prob = Problem(dev; nx, Lx, κ, dt, stepper="RK4")
  c₀ = @. cos(k₀ * prob.grid.x)
  set_c!(prob, c₀)

  getsol(prob) = prob.sol
  soldiag = Diagnostic(getsol, prob, nsteps=1)

  stepforward!(prob, soldiag, 1)

  return isapprox(soldiag.data[soldiag.i], prob.sol)
end

function test_extenddiagnostic(dev::Device=CPU(); nx=6, Lx=2π, κ=1e-2)
  k₀ = 2π / Lx
  dt = 1e-2 / (κ * k₀^2) # time-scale for diffusive decay dynamics are resolved
  prob = Problem(dev; nx, Lx, κ, dt, stepper="RK4")
  c₀ = @. sin(k₀ * prob.grid.x)
  set_c!(prob, c₀)

  getsol(prob) = prob.sol
  
  soldiagnostic₁ = Diagnostic(getsol, prob, nsteps=1)
  soldiagnostic₂ = Diagnostic(getsol, prob, nsteps=4)
  
  nsteps_initially₁ = length(soldiagnostic₁.t)
  nsteps_extend₁ = 12
  FourierFlows.extend!(soldiagnostic₁, nsteps_extend1)
  
  nsteps_initially₂ = length(soldiagnostic₂.t)
  FourierFlows.extend!(soldiagnostic₂)
  
  return length(soldiagnostic₁.t) == nsteps_initially₁ + nsteps_extend₁ &&
         length(soldiagnostic₂.t) == 2 * nsteps_initially₂
end

function test_incrementdiagnostic(dev::Device=CPU(); nx=6, Lx=2π, κ=1e-2)
  k₀ = 2π / Lx
  dt = 1e-2 / (κ * k₀^2) # time-scale for diffusive decay dynamics are resolved
  prob = Problem(dev; nx, Lx, κ, dt, stepper="ETDRK4")
  c₀ = @. sin(k₀ * prob.grid.x)
  set_c!(prob, c₀)
  
  dummydiagnostic₁(prob) = 2.0
  dummydiagnostic₂(prob) = 2.0im
  
  diagnostic₁ = Diagnostic(dummydiagnostic₁, prob, nsteps=10)
  diagnostic₂ = Diagnostic(dummydiagnostic₂, prob, nsteps=10)
 
  diagnostics = [diagnostic₁, diagnostic₂]
  
  FourierFlows.increment!(diagnostics)
  
  return diagnostic₁.data[2]==2.0 && diagnostic₂.data[2]==2.0im
end

function test_getindex(dev::Device=CPU(); nx=6, Lx=2π, κ=1e-2)
  k₀ = 2π / Lx
  dt = 1e-2 / (κ * k₀^2) # time-scale for diffusive decay dynamics are resolved
  prob = Problem(dev; nx, Lx, κ, dt, stepper="ETDRK4")
  c₀ = @. sin(k₀ * prob.grid.x)
  set_c!(prob, c₀)
  
  dummydiagnostic₁(prob) = 2.0
  dummydiagnostic₂(prob) = 2.0im
  
  nsteps = 10
  
  diagnostic₁ = Diagnostic(dummydiagnostic₁, prob, nsteps=nsteps)
  diagnostic₂ = Diagnostic(dummydiagnostic₂, prob, nsteps=nsteps)
 
  diagnostics = [diagnostic₁, diagnostic₂]
  
  ntimes=7
  for j=1:ntimes
    FourierFlows.increment!(diagnostics)
  end
  
  return getindex(diagnostic₁, 2:5) == 2 * ones(4) &&
         getindex(diagnostic₂, 2:5) == 2im * ones(4)
end
