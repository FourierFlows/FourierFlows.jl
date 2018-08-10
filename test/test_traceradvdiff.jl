# traceradvdiff.jl's test function

import FourierFlows.TracerAdvDiff

"""
    test_constvel(; kwargs...)

Advects a gaussian concentration c0(x, y, t) with a constant velocity flow
u(x, y, t) = uvel and v(x, y, t) = vvel and compares the final state with
cfinal = c0(x-uvel*tfinal, y-vvel*tfinal)
"""
function test_constvel(stepper, dt, nsteps)

  nx, Lx = 128, 2π
  uvel, vvel = 0.2, 0.1
  u(x, y) = uvel
  v(x, y) = vvel

  prob = TracerAdvDiff.Problem(; nx=nx, Lx=Lx, kap=0.0, u=u, v=v, dt=dt, stepper=stepper)
  s, v, p, g = prob.state, prob.vars, prob.params, prob.grid

  σ = 0.1
  c0func(x, y) = 0.1*exp.(-(x.^2+y.^2)/(2σ^2))

  c0 = c0func.(g.X, g.Y)
  tfinal = nsteps*dt
  cfinal = c0func.(g.X .- uvel*tfinal, g.Y .- vvel*tfinal)

  TracerAdvDiff.set_c!(prob, c0)

  stepforward!(prob, nsteps)
  TracerAdvDiff.updatevars!(prob)

  isapprox(cfinal, v.c, rtol=g.nx*g.ny*nsteps*1e-12)
end


"""
    test_timedependenttvel(; kwargs...)

Advects a gaussian concentration c0(x, y, t) with a time-varying velocity flow
u(x, y, t) = uvel and v(x, y, t) = vvel*sign(-t+tfinal/2) and compares the final
state with cfinal = c0(x-uvel*tfinal, y)
"""
function test_timedependentvel(stepper, dt, tfinal)

  nx, Lx = 128, 2π
  nsteps = round(Int, tfinal/dt)

  if !isapprox(tfinal, nsteps*dt, rtol=1e-12)
      error("tfinal is not multiple of dt")
  end

  uvel, vvel = 0.2, 0.1
  u(x, y, t) = uvel
  v(x, y, t) = t <= tfinal/2 ? vvel : -vvel

  prob = TracerAdvDiff.ConstDiffProblem(; nx=nx, Lx=Lx, kap=0.0, u=u, v=v, dt=dt, stepper=stepper)
  s, v, p, g = prob.state, prob.vars, prob.params, prob.grid

  σ = 0.1
  c0func(x, y) = 0.1*exp.(-(x.^2+y.^2)/(2σ^2))

  c0 = c0func.(g.X, g.Y)
  tfinal = nsteps*dt
  cfinal = c0func.(g.X .- uvel*tfinal, g.Y)

  TracerAdvDiff.set_c!(prob, c0func)

  stepforward!(prob, nsteps)
  TracerAdvDiff.updatevars!(prob)

  isapprox(cfinal, v.c, rtol=g.nx*g.ny*nsteps*1e-12)
end


"""
    test_diffusion(; kwargs...)

Diffuses a gaussian concentration c0(x, y, t) and compares the final state with
the analytic solution of the heat equation, cfinal
"""
function test_diffusion(stepper, dt, tfinal; steadyflow = true)

  nx  = 128
  Lx  = 2π
  kap = 0.01
  nsteps = round(Int, tfinal/dt)

  if !isapprox(tfinal, nsteps*dt, rtol=1e-12)
      error("tfinal is not multiple of dt")
  end

  grid = TwoDGrid(nx, Lx)
  prob = TracerAdvDiff.ConstDiffProblem(; steadyflow=steadyflow, grid=grid, nx=nx,
    Lx=Lx, kap=kap, dt=dt, stepper=stepper)
  s, v, p, g = prob.state, prob.vars, prob.params, prob.grid

  c0ampl, σ = 0.1, 0.1
  c0func(x, y) = c0ampl*exp.(-(x.^2+y.^2)/(2σ^2))

  c0 = c0func.(g.X, g.Y)
  tfinal = nsteps*dt
  σt = sqrt(2*kap*tfinal + σ^2)
  cfinal = c0ampl*σ^2/σt^2 * exp.(-(g.X.^2 + g.Y.^2)/(2*σt^2))

  TracerAdvDiff.set_c!(prob, c0)

  stepforward!(prob, nsteps)
  TracerAdvDiff.updatevars!(prob)

  isapprox(cfinal, v.c, rtol=g.nx*g.ny*nsteps*1e-12)
end


"""
    test_hyperdiffusion(; kwargs...)

Diffuses a gaussian concentration c0(x, y, t) using hyperdiffusivity and
compares the final state with the analytic solution of the heat equation, cfinal
"""
function test_hyperdiffusion(stepper, dt, tfinal; steadyflow = true)

    nx  = 128
    Lx  = 2π
    kap = 0.0   # no diffusivity
    eta = kap   # no diffusivity
   kaph = 0.01  # hyperdiffusivity coeff
  nkaph = 1     # nkaph=1 converts hyperdiffusivity to plain diffusivity
                # so we can compare with the analytic solution of heat equation

  nsteps = round(Int, tfinal/dt)

  if !isapprox(tfinal, nsteps*dt, rtol=1e-12)
      error("tfinal is not multiple of dt")
  end

   g = TwoDGrid(nx, Lx)

  u, v = 0*g.X, 0*g.X

  vs = TracerAdvDiff.Vars(g)
  pr = TracerAdvDiff.ConstDiffSteadyFlowParams(eta, kap, kaph, nkaph, u, v, g)
  eq = TracerAdvDiff.Equation(pr, g)
  ts = FourierFlows.autoconstructtimestepper(stepper, dt, eq.LC, g)
  prob = FourierFlows.Problem(g, vs, pr, eq, ts)

  c0ampl, σ = 0.1, 0.1
  c0func(x, y) = c0ampl*exp.(-(x.^2+y.^2)/(2σ^2))

  c0 = c0func.(g.X, g.Y)
  tfinal = nsteps*dt
  σt = sqrt(2*kaph*tfinal + σ^2)
  cfinal = c0ampl*σ^2/σt^2 * exp.(-(g.X.^2 + g.Y.^2)/(2*σt^2))

  TracerAdvDiff.set_c!(prob, c0)

  stepforward!(prob, nsteps)
  TracerAdvDiff.updatevars!(prob)

  isapprox(cfinal, vs.c, rtol=g.nx*g.ny*nsteps*1e-12)
end


# --
# Run tests
# --

stepper = "RK4"

dt, nsteps  = 1e-2, 40
@test test_constvel(stepper, dt, nsteps)

dt, tfinal  = 0.002, 0.1
@test test_timedependentvel(stepper, dt, tfinal)

dt, tfinal  = 0.005, 0.1
@test test_diffusion(stepper, dt, tfinal; steadyflow=true)

dt, tfinal  = 0.005, 0.1
@test test_diffusion(stepper, dt, tfinal; steadyflow=false)

dt, tfinal  = 0.005, 0.1
@test test_hyperdiffusion(stepper, dt, tfinal)
