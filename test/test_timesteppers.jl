function diffusionproblem(stepper; nx=6, Lx=2π, kap=1e-2, nsteps=1000)
  k1 = 2π/Lx
   τ = 1/(kap*k1^2) # time-scale for diffusive decay
  dt = 1e-9 * τ # dynamics are resolved.

  prob = Problem(nx=nx, Lx=Lx, kap=kap, dt=dt, stepper=stepper)
  g = prob.grid

  t1 = dt*nsteps
  c0 = @. sin(k1*g.x)
  c1 = @. exp(-kap*k1^2*t1) * c0 # analytical solution

  set_c!(prob, c0)
  tcomp = @elapsed stepforward!(prob, nsteps)
  updatevars!(prob)

  prob, c0, c1, nsteps, tcomp
end

function diffusiontest(stepper; kwargs...)
  prob, c0, c1, nsteps, tcomp = diffusionproblem(stepper; kwargs...)
  normmsg = "$stepper: relative error ="
  @printf("% 40s %.2e (%.3f s)\n", normmsg, norm(c1-prob.vars.c)/norm(c1), tcomp)
  isapprox(c1, prob.vars.c, rtol=nsteps*rtol_timesteppers)
end
