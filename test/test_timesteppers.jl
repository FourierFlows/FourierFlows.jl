function constantdiffusionproblem(stepper; nx=128, Lx=2π, kappa=1e-2, nsteps=1000, dev=CPU())
   τ = 1/kappa  # time-scale for diffusive decay
  dt = 1e-9 * τ # dynamics are resolved

  prob = Problem(nx=nx, Lx=Lx, kappa=kappa, dt=dt, stepper=stepper, dev=dev)
  g = prob.grid

  # a gaussian initial condition c(x, t=0)
  c0ampl, σ = 0.01, 0.2
  c0func(x) = @. c0ampl*exp(-x^2/(2σ^2))
  c0 = c0func.(g.x)

  # analytic solution for for 1D heat equation with constant κ
  tfinal = nsteps*dt
  σt = sqrt(2*kappa*tfinal + σ^2)
  cfinal = @. c0ampl*σ/σt * exp(-g.x^2/(2*σt^2))

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
  g = prob.grid

  # a gaussian initial condition c(x, t=0)
  c0ampl, σ = 0.01, 0.2
  c0func(x) = @. c0ampl*exp(-x^2/(2σ^2))
  c0 = c0func.(g.x)

  # analytic solution for for 1D heat equation with constant κ
  tfinal = nsteps*dt
  σt = sqrt(2*kappa[1]*tfinal + σ^2)
  cfinal = @. c0ampl*σ/σt * exp(-g.x^2/(2*σt^2))

  set_c!(prob, c0)
  tcomp = @elapsed stepforward!(prob, nsteps)
  updatevars!(prob)

  prob, c0, cfinal, nsteps, tcomp
end


function constantdiffusiontest(stepper, dev::Device=CPU(); kwargs...)
  prob, c0, c1, nsteps, tcomp = constantdiffusionproblem(stepper; kwargs...)
  normmsg = "$stepper: relative error ="
  @printf("% 40s %.2e (%.3f s)\n", normmsg, norm(c1-prob.vars.c)/norm(c1), tcomp)
  isapprox(c1, prob.vars.c, rtol=nsteps*rtol_timesteppers)
end

function varyingdiffusiontest(stepper, dev::Device=CPU(); kwargs...)
  prob, c0, c1, nsteps, tcomp = varyingdiffusionproblem(stepper; kwargs...)
  normmsg = "$stepper: relative error ="
  @printf("% 40s %.2e (%.3f s)\n", normmsg, norm(c1-prob.vars.c)/norm(c1), tcomp)
  isapprox(c1, prob.vars.c, rtol=nsteps*rtol_timesteppers)
end
