function constantdiffusionproblem(stepper; nx=128, Lx=2π, kappa=1e-2, nsteps=1000, dev=CPU())
   τ = 1/kappa  # time-scale for diffusive decay
  dt = 1e-9 * τ # dynamics are resolved

  prob = Problem(nx=nx, Lx=Lx, kappa=kappa, dt=dt, stepper=stepper, dev=dev)
  
  # a gaussian initial condition c(x, t=0)
  c0ampl, σ = 0.01, 0.2
  c0func(x) = c0ampl * exp(-x^2/(2σ^2))
  c0 = c0func.(prob.grid.x)

  # analytic solution for for 1D heat equation with constant κ
  tfinal = nsteps*dt
  σt = sqrt(2*kappa*tfinal + σ^2)
  cfinal = @. c0ampl * σ/σt * exp(-prob.grid.x^2/(2*σt^2))

  set_c!(prob, c0)
  tcomp = @elapsed stepforward!(prob, nsteps)
  updatevars!(prob)

  prob, c0, cfinal, nsteps, tcomp
end

function varyingdiffusionproblem(stepper; nx=128, Lx=2π, kappa=1e-2, nsteps=1000, dev=CPU())
   τ = 1/kappa  # time-scale for diffusive decay
  dt = 1e-9 * τ # dynamics are resolved

  kappa = kappa*ones(nx) # this is actually a constant diffusion but defining it
                         # as an array makes stepforward! call function calcN!
                         # instead of just the linear coefficients L*sol

  prob = Problem(nx=nx, Lx=Lx, kappa=kappa, dt=dt, stepper=stepper, dev=dev)
  
  # a gaussian initial condition c(x, t=0)
  c0ampl, σ = 0.01, 0.2
  c0func(x) = c0ampl * exp(-x^2/(2σ^2))
  c0 = c0func.(prob.grid.x)

  # analytic solution for for 1D heat equation with constant κ
  tfinal = nsteps*dt
  σt = sqrt(2*kappa[1]*tfinal + σ^2)
  cfinal = @. c0ampl * σ/σt * exp(-prob.grid.x^2/(2*σt^2))

  set_c!(prob, c0)
  tcomp = @elapsed stepforward!(prob, nsteps)
  updatevars!(prob)

  prob, c0, cfinal, nsteps, tcomp
end

function constantdiffusionproblem_stoptime(stepper; nx=128, Lx=2π, kappa=1e-2, stop_time=1e-4+1e-6/π, dev=CPU())
   τ = 1/kappa  # time-scale for diffusive decay
  dt = 1e-9 * τ # dynamics are resolved

  prob = Problem(nx=nx, Lx=Lx, kappa=kappa, dt=dt, stepper=stepper, dev=dev)
  
  # a gaussian initial condition c(x, t=0)
  c0ampl, σ = 0.01, 0.2
  c0func(x) = c0ampl * exp(-x^2/(2σ^2))
  c0 = c0func.(prob.grid.x)

  # analytic solution for for 1D heat equation with constant κ
  tfinal = stop_time
  σt = sqrt(2*kappa*tfinal + σ^2)
  cfinal = @. c0ampl * σ/σt * exp(-prob.grid.x^2/(2*σt^2))

  set_c!(prob, c0)
  tcomp = @elapsed stepforward!(prob, tfinal)
  updatevars!(prob)

  prob, c0, cfinal, tfinal, tcomp
end

function constantdiffusiontest(stepper, dev::Device=CPU(); kwargs...)
  prob, c0, cfinal, nsteps, tcomp = constantdiffusionproblem(stepper; kwargs...)
  normmsg = "$stepper: relative error ="
  @printf("% 40s %.2e (%.3f s)\n", normmsg, norm(cfinal-Array(prob.vars.c))/norm(cfinal), tcomp)
  isapprox(cfinal, Array(prob.vars.c), rtol=nsteps*rtol_timesteppers)
end

function constantdiffusiontest_stoptime(stepper, dev::Device=CPU(); kwargs...)
  prob, c0, cfinal, tfinal, tcomp = constantdiffusionproblem_stoptime(stepper; kwargs...)
  normmsg = "$stepper: relative error ="
  @printf("% 40s %.2e (%.3f s)\n", normmsg, norm(cfinal-Array(prob.vars.c))/norm(cfinal), tcomp)
  isapprox(cfinal, Array(prob.vars.c), rtol=prob.clock.step*rtol_timesteppers)
end

function varyingdiffusiontest(stepper, dev::Device=CPU(); kwargs...)
  prob, c0, cfinal, nsteps, tcomp = varyingdiffusionproblem(stepper; kwargs...)
  normmsg = "$stepper: relative error ="
  @printf("% 40s %.2e (%.3f s)\n", normmsg, norm(cfinal-Array(prob.vars.c))/norm(cfinal), tcomp)
  isapprox(cfinal, Array(prob.vars.c), rtol=nsteps*rtol_timesteppers)
end
