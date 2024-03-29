gausian_solution(x, t; c₀=0.01, σ=0.2, κ=1e-2) = c₀ * σ / sqrt(σ^2 + 2κ*t) * exp(-x^2 / (2(σ^2 + 2κ*t)))

construct_diffusion_problem(stepper, κ, dt; nx=128, Lx=2π, dev=CPU()) =
  Problem(dev; nx, Lx, κ, dt, stepper)

function constantdiffusion_stepforward(prob, nsteps; c₀=0.01, σ=0.2, κ=1e-2)
  # a gaussian initial condition c(x, t=0)
  c_initial = @. gausian_solution(prob.grid.x, 0; c₀, σ, κ)

  # analytic solution for for 1D heat equation with constant κ
  t_final = nsteps*prob.clock.dt
  c_final = @. gausian_solution(prob.grid.x, t_final; c₀, σ, κ)

  set_c!(prob, c_initial)
  t_compute = @elapsed stepforward!(prob, nsteps)
  updatevars!(prob)

  return c_initial, c_final, prob, t_compute
end

function constantdiffusion_step_until(prob, t_final; c₀=0.01, σ=0.2, κ=1e-2)
  # a gaussian initial condition c(x, t=0)
  c_initial = @. gausian_solution(prob.grid.x, 0; c₀, σ, κ)

  # analytic solution for for 1D heat equation with constant κ
  c_final = @. gausian_solution(prob.grid.x, t_final; c₀, σ, κ)

  set_c!(prob, c_initial)
  t_compute = @elapsed step_until!(prob, t_final)
  updatevars!(prob)

  return c_initial, c_final, prob, t_compute
end

nx = 128
 κ = 1e-2
dt = 1e-9 * 1/κ # make sure diffusive decay dynamics are resolved

function constantdiffusiontest_stepforward(stepper, dev::Device=CPU(); kwargs...)
  nsteps = 1000
    
  prob = construct_diffusion_problem(stepper, κ, dt; nx, Lx=2π, dev)
  c_initial, c_final, prob, t_compute = constantdiffusion_stepforward(prob, nsteps; c₀=0.01, σ=0.2, κ)
  normmsg = "$stepper: relative error ="
  @printf("% 40s %.2e (%.3f s)\n", normmsg, norm(c_final-Array(prob.vars.c))/norm(c_final), t_compute)

  return isapprox(c_final, Array(prob.vars.c), rtol=prob.clock.step*rtol_timesteppers)
end

function varyingdiffusiontest_stepforward(stepper, dev::Device=CPU(); kwargs...)
  nsteps = 1000
  
  prob = construct_diffusion_problem(stepper, κ * ones(nx), dt; nx, Lx=2π, dev)
  c_initial, c_final, prob, t_compute = constantdiffusion_stepforward(prob, nsteps; c₀=0.01, σ=0.2, κ)
  normmsg = "$stepper: relative error ="
  @printf("% 40s %.2e (%.3f s)\n", normmsg, norm(c_final-Array(prob.vars.c))/norm(c_final), t_compute)

  return isapprox(c_final, Array(prob.vars.c), rtol=prob.clock.step*rtol_timesteppers)
end

function constantdiffusiontest_step_until(stepper, dev::Device=CPU(); kwargs...)
  t_final = 1000*dt + 1e-6/π # make sure t_final is not integer multiple of dt
  
  prob = construct_diffusion_problem(stepper, κ, dt; nx, Lx=2π, dev)
  c_initial, c_final, prob, t_compute = constantdiffusion_step_until(prob, t_final; c₀=0.01, σ=0.2, κ)
  normmsg = "$stepper: relative error ="
  @printf("% 40s %.2e (%.3f s)\n", normmsg, norm(c_final-Array(prob.vars.c))/norm(c_final), t_compute)

  return isapprox(c_final, Array(prob.vars.c), rtol=prob.clock.step*rtol_timesteppers)
end
