import FourierFlows.BarotropicQG

using FFTW, Statistics, Random
import Statistics: mean


# -----------------------------------------------------------------------------
# BAROQG's TEST FUNCTIONS

"""
    test_baroQG_RossbyWave(; kwargs...)

Evolvesa a Rossby wave and compares with the analytic solution.
"""
function test_baroQG_RossbyWave(stepper, dt, nsteps)
    nx = 64
  beta = 2.0
    Lx = 2π
    mu = 0.0
    nu = 0.0

  # the following if statement is called so that all the cases of
  # Problem() fuction are tested
  if stepper=="ForwardEuler"
    eta = zeros(nx, nx)
  else
    eta(x, y) = 0*x
  end

  prob = BarotropicQG.InitialValueProblem(nx=nx, Lx=Lx, eta=eta, beta=beta, mu=mu, nu=nu, stepper=stepper, dt=dt)
  s, v, p, g = prob.state, prob.vars, prob.params, prob.grid

  # the Rossby wave initial condition
   ampl = 1e-2
  kwave = 3.0*2π/g.Lx
  lwave = 2.0*2π/g.Ly
      ω = -p.beta*kwave/(kwave^2.0 + lwave^2.0)
     ζ0 = ampl*cos.(kwave*g.X).*cos.(lwave*g.Y)
    ζ0h = rfft(ζ0)

  BarotropicQG.set_zeta!(prob, ζ0)

  stepforward!(prob, nsteps)
  dealias!(s.sol, g)
  BarotropicQG.updatevars!(prob)

  ζ_theory = ampl*cos.(kwave*(g.X .- ω/kwave*s.t)) .* cos.(lwave*g.Y)

  isapprox(ζ_theory, v.zeta, rtol=g.nx*g.ny*nsteps*1e-12)
end

"""
    test_stochasticforcingbudgets(; kwargs...)

Tests if the energy budgets are closed for BarotropicQG with stochastic forcing.
"""
function test_stochasticforcingbudgets(; n=256, dt=0.01, L=2π, nu=1e-7, nnu=2, mu=1e-1, message=false)
  n, L  = 256, 2π
  nu, nnu = 1e-7, 2
  mu = 1e-1
  dt, tf = 0.005, 0.1/mu
  nt = round(Int, tf/dt)
  ns = 1

  # Forcing
  kf, dkf = 12.0, 2.0
  σ = 0.1
  gr  = TwoDGrid(n, L)

  force2k = zero(gr.Kr)
  @. force2k = exp.(-(sqrt.(gr.KKrsq)-kf).^2/(2*dkf^2))
  @. force2k[gr.KKrsq .< 2.0^2 ] = 0
  @. force2k[gr.KKrsq .> 20.0^2 ] = 0
  @. force2k[gr.Kr.<2π/L] = 0
  σ0 = FourierFlows.parsevalsum(force2k.*gr.invKKrsq/2.0, gr)/(gr.Lx*gr.Ly)
  force2k .= σ/σ0 * force2k

  Random.seed!(1234)

  function calcFq!(F, sol, t, s, v, p, g)
    eta = exp.(2π*im*rand(Float64, size(sol)))/sqrt(s.dt)
    eta[1, 1] = 0
    @. F = eta .* sqrt(force2k)
    nothing
  end

  prob = BarotropicQG.ForcedProblem(nx=n, Lx=L, nu=nu, nnu=nnu, mu=mu, dt=dt,
   stepper="RK4", calcFq=calcFq!)

  s, v, p, g, eq, ts = prob.state, prob.vars, prob.params, prob.grid, prob.eqn, prob.ts;

  BarotropicQG.set_zeta!(prob, 0*g.X)
  E = Diagnostic(FourierFlows.BarotropicQG.energy,      prob, nsteps=nt)
  D = Diagnostic(FourierFlows.BarotropicQG.dissipation, prob, nsteps=nt)
  R = Diagnostic(FourierFlows.BarotropicQG.drag,        prob, nsteps=nt)
  W = Diagnostic(FourierFlows.BarotropicQG.work,        prob, nsteps=nt)
  diags = [E, D, W, R]

  # Step forward

  stepforward!(prob, diags, round(Int, nt))

  BarotropicQG.updatevars!(prob)

  cfl = prob.ts.dt*maximum([maximum(v.v)/g.dx, maximum(v.u)/g.dy])

  E, D, W, R = diags

  t = round(mu*prob.state.t, digits=2)

  i₀ = 1
  dEdt = (E[(i₀+1):E.count] - E[i₀:E.count-1])/prob.ts.dt
  ii = (i₀):E.count-1
  ii2 = (i₀+1):E.count

  # dEdt = W - D - R?
  # If the Ito interpretation was used for the work
  # then we need to add the drift term
  # total = W[ii2]+σ - D[ii] - R[ii]      # Ito
  total = W[ii2] - D[ii] - R[ii]        # Stratonovich

  residual = dEdt - total

  if message
    println("step: %04d, t: %.1f, cfl: %.3f, time: %.2f s\n", prob.step, prob.t, cfl, tc)
  end
  # println(mean(abs.(residual)))
  isapprox(mean(abs.(residual)), 0, atol=1e-4)
end

"""
    testnonlineartermsQGPV(dt, stepper; kwargs...)

Tests the advection term in the twodturb module by timestepping a
test problem with timestep dt and timestepper identified by the string stepper.
The test problem is derived by picking a solution ζf (with associated
streamfunction ψf) for which the advection term J(ψf, ζf) is non-zero. Next, a
forcing Ff is derived according to Ff = ∂ζf/∂t + J(ψf, ζf) - nuΔζf. One solution
to the vorticity equation forced by this Ff is then ζf. (This solution may not
be realized, at least at long times, if it is unstable.)
"""
function testnonlineartermsQGPV(dt, stepper; n=128, L=2π, nu=1e-2, nnu=1, mu=0.0, message=false)
  n, L  = 128, 2π
  nu, nnu = 1e-2, 1
  mu = 0.0
  tf = 1.0
  nt = round(Int, tf/dt)

  gr  = TwoDGrid(n, L)
  x, y = gr.X, gr.Y

  psif = @. sin(2x)*cos(2y) + 2sin(x)*cos(3y)
  qf = @. -8sin(2x)*cos(2y) - 20sin(x)*cos(3y)

  Ff = @. -(
    nu*( 64sin(2x)*cos(2y) + 200sin(x)*cos(3y) )
    + 8*( cos(x)*cos(3y)*sin(2x)*sin(2y) - 3cos(2x)*cos(2y)*sin(x)*sin(3y) )
  )

  Ffh = rfft(Ff)

  # Forcing
  function calcFq!(Fqh, sol, t, s, v, p, g)
    Fqh .= Ffh
    nothing
  end

  prob = BarotropicQG.ForcedProblem(nx=n, Lx=L, nu=nu, nnu=nnu, mu=mu, dt=dt, stepper=stepper, calcFq=calcFq!)
  s, v, p, g, eq, ts = prob.state, prob.vars, prob.params, prob.grid, prob.eqn, prob.ts
  BarotropicQG.set_zeta!(prob, qf)

  # Step forward
  stepforward!(prob, round(Int, nt))
  BarotropicQG.updatevars!(prob)
  isapprox(v.q, qf, rtol=1e-13)
end

"""
    testnonlineartermsU(dt, stepper; kwargs...)

Tests the form stress term that forces the domain-averaged zonal flow U(t).
"""
function testnonlineartermsU(dt, stepper; n=128, L=2π, nu=0.0, nnu=1, mu=0.0, message=false)
  n, L  = 128, 2π
  nu, nnu = 1e-2, 1
  mu = 0.0
  tf = 1.0
  nt = 1



  gr  = TwoDGrid(n, L)
  x, y = gr.X, gr.Y

  zetai = -20*sin.(10*x).*cos.(10*y)
  topoPV(x, y) = cos.(10x).*cos.(10y)
  F(t) = 0 #no forcing

  answer = 0.25 # this is what <v*eta> should be

  prob = BarotropicQG.ForcedProblem(nx=n, Lx=L, nu=nu, nnu=nnu, mu=mu, dt=dt, stepper=stepper, eta=topoPV, calcFU = F, calcFq=nothing)
  s, v, p, g, eq, ts = prob.state, prob.vars, prob.params, prob.grid, prob.eqn, prob.ts
  BarotropicQG.set_zeta!(prob, zetai)
  BarotropicQG.updatevars!(prob)

  # Step forward
  stepforward!(prob, nt)
  isapprox(prob.ts.N[1, 1], answer, rtol=1e-13)
end

function testenergyenstrophy()
  nx, Lx  = 64, 2π
  ny, Ly  = 64, 3π
  g  = TwoDGrid(nx, Lx, ny, Ly)
  k0 = g.k[2] # fundamental wavenumber
  l0 = g.l[2] # fundamental wavenumber
  x, y = g.X, g.Y

    eta = @. cos(10*k0*x)*cos(10*l0*y)
   psi0 = @. sin(2*k0*x)*cos(2*l0*y) + 2sin(k0*x)*cos(3*l0*y)
  zeta0 = @. -((2*k0)^2+(2*l0)^2)*sin(2*k0*x)*cos(2*l0*y) - (k0^2+(3*l0)^2)*2sin(k0*x)*cos(3*l0*y)

  prob = BarotropicQG.InitialValueProblem(nx=nx, Lx=Lx, ny=ny, Ly=Ly, eta = eta, stepper="ForwardEuler")
  s, v, p, g, eq, ts = prob.state, prob.vars, prob.params, prob.grid, prob.eqn, prob.ts
  BarotropicQG.set_zeta!(prob, zeta0)
  BarotropicQG.updatevars!(prob)

  energyzeta0 = FourierFlows.BarotropicQG.energy(prob)
  enstrophyzeta0 = FourierFlows.BarotropicQG.enstrophy(prob)

  isapprox(energyzeta0, 29.0/9, rtol=1e-13) && isapprox(enstrophyzeta0, 2701.0/162, rtol=1e-13)
end

function testenergyenstrophy00()
  nx, Lx  = 64, 2π
  ny, Ly  = 96, 3π
  g  = TwoDGrid(nx, Lx, ny, Ly)
  k0 = g.k[2] # fundamental wavenumber
  l0 = g.l[2] # fundamental wavenumber
  x, y = g.X, g.Y

  calcFU(t) = 0.0
  eta(x, y) = cos.(10x).*cos.(10y)
  psi0 = @. sin(2*k0*x)*cos(2*l0*y) + 2sin(k0*x)*cos(3*l0*y)
 zeta0 = @. -((2*k0)^2+(2*l0)^2)*sin(2*k0*x)*cos(2*l0*y) - (k0^2+(3*l0)^2)*2sin(k0*x)*cos(3*l0*y)
  beta = 10.0
  U = 1.2

  prob = BarotropicQG.ForcedProblem(nx=nx, Lx=Lx, ny=ny, Ly=Ly, beta=beta, eta=eta, calcFU = calcFU, stepper="ForwardEuler")
  s, v, p, g, eq, ts = prob.state, prob.vars, prob.params, prob.grid, prob.eqn, prob.ts

  BarotropicQG.set_zeta!(prob, zeta0)
  BarotropicQG.set_U!(prob, U)
  BarotropicQG.updatevars!(prob)

  energyU = FourierFlows.BarotropicQG.energy00(prob)
  enstrophyU = FourierFlows.BarotropicQG.enstrophy00(prob)

  energyzeta0 = FourierFlows.BarotropicQG.energy(prob)
  enstrophyzeta0 = FourierFlows.BarotropicQG.enstrophy(prob)

  (isapprox(energyU, 0.5*U^2, rtol=1e-13) &&
    isapprox(enstrophyU, beta*U, rtol=1e-13) &&
    isapprox(energyzeta0, 29.0/9, rtol=1e-13) &&
    isapprox(enstrophyzeta0, 2701.0/162, rtol=1e-13)
  )
end

# -----------------------------------------------------------------------------
# Running the tests


dt, nsteps, stepper  = 1e-2, 20, "ETDRK4"
@test test_baroQG_RossbyWave("ETDRK4", dt, nsteps)

dt, nsteps  = 1e-2, 20;
@test test_baroQG_RossbyWave("FilteredETDRK4", dt, nsteps)

dt, nsteps  = 1e-2, 20
@test test_baroQG_RossbyWave("RK4", dt, nsteps)

dt, nsteps  = 1e-2, 20
@test test_baroQG_RossbyWave("FilteredRK4", dt, nsteps)

dt, nsteps  = 1e-3, 200
@test test_baroQG_RossbyWave("AB3", dt, nsteps)

dt, nsteps  = 1e-3, 200
@test test_baroQG_RossbyWave("FilteredAB3", dt, nsteps)

dt, nsteps  = 1e-4, 2000
@test test_baroQG_RossbyWave("ForwardEuler", dt, nsteps)

dt, nsteps  = 1e-4, 2000
@test test_baroQG_RossbyWave("FilteredForwardEuler", dt, nsteps)

@test test_stochasticforcingbudgets()

@test testnonlineartermsQGPV(0.0005, "ForwardEuler")
@test testnonlineartermsU(0.01, "ForwardEuler")

@test testenergyenstrophy()
@test testenergyenstrophy00()
