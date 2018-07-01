import FourierFlows.TwoDTurb

cfl(prob) = maximum([maximum(abs.(prob.vars.U)), maximum(abs.(prob.vars.V))]*prob.ts.dt/prob.grid.dx)

function lambdipoletest(n, dt; L=2π, Ue=1, Re=L/20, nu=0.0, nnu=1, ti=L/Ue*0.01, nm=3)
  nt = round(Int, ti/dt)
  prob = TwoDTurb.InitialValueProblem(nx=n, Lx=L, nu=nu, nnu=nnu, dt=dt, stepper="FilteredRK4")
  q0 = FourierFlows.lambdipole(Ue, Re, prob.grid)
  TwoDTurb.set_q!(prob, q0)

  xq = zeros(nm)   # centroid of abs(q)
  Ue_m = zeros(nm) # measured dipole speed
  x, y, q = prob.grid.X, prob.grid.Y, prob.vars.q # nicknames

  for i = 1:nm # step forward
    tic()
    stepforward!(prob, nt)
    TwoDTurb.updatevars!(prob)
    xq[i] = mean(abs.(q).*x) / mean(abs.(q))

    if i > 1
      Ue_m[i] = (xq[i]-xq[i-1]) / ((nt-1)*dt)
    end
  end
  isapprox(Ue, mean(Ue_m[2:end]), rtol=1e-2)
end

function stochasticforcingbudgetstest(; n=256, dt=0.01, L=2π, nu=1e-7, nnu=2, mu=1e-1, nmu=0, message=false)
  n, L  = 256, 2π
  nu, nnu = 1e-7, 2
  mu, nmu = 1e-1, 0
  dt, tf = 0.005, 0.1/mu
  nt = round(Int, tf/dt)
  ns = 1

  # Forcing
  kf, dkf = 12.0, 2.0
  σ = 0.1
  gr  = TwoDGrid(n, L)

  force2k = exp.(-(sqrt.(gr.KKrsq)-kf).^2/(2*dkf^2))
  force2k[gr.KKrsq .< 2.0^2 ] = 0
  force2k[gr.KKrsq .> 20.0^2 ] = 0
  force2k[gr.Kr.<2π/L] = 0
  σ0 = FourierFlows.parsevalsum(force2k.*gr.invKKrsq/2.0, gr)/(gr.Lx*gr.Ly)
  force2k .= σ/σ0 * force2k

  srand(1234)

  function calcF!(F, sol, t, s, v, p, g)
    eta = exp.(2π*im*rand(size(sol)))/sqrt(s.dt)
    eta[1, 1] = 0.0
    @. F = eta .* sqrt(force2k)
    nothing
  end

  prob = TwoDTurb.ForcedProblem(nx=n, Lx=L, nu=nu, nnu=nnu, mu=mu, nmu=nmu, dt=dt,
   stepper="RK4", calcF=calcF!)

  s, v, p, g, eq, ts = prob.state, prob.vars, prob.params, prob.grid, prob.eqn, prob.ts;

  TwoDTurb.set_q!(prob, 0*g.X)
  E = Diagnostic(FourierFlows.TwoDTurb.energy,      prob, nsteps=nt)
  D = Diagnostic(FourierFlows.TwoDTurb.dissipation, prob, nsteps=nt)
  R = Diagnostic(FourierFlows.TwoDTurb.drag,        prob, nsteps=nt)
  W = Diagnostic(FourierFlows.TwoDTurb.work,        prob, nsteps=nt)
  diags = [E, D, W, R]

  # Step forward
  stepforward!(prob, diags, round(Int, nt))

  TwoDTurb.updatevars!(prob)

  cfl = prob.ts.dt*maximum(
    [maximum(v.V)/g.dx, maximum(v.U)/g.dy])

  E, D, W, R = diags

  t = round(mu*prob.state.t, 2)

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
    @printf("step: %04d, t: %.1f, cfl: %.3f, time: %.2f s\n", prob.step, prob.t, cfl, tc)
  end
  # println(mean(abs.(residual)))
  isapprox(mean(abs.(residual)), 0, atol=1e-4)
end

"""
    testnonlinearterms(dt, stepper; kwargs...)

Tests the advection term in the twodturb module by timestepping a
test problem with timestep dt and timestepper identified by the string stepper.
The test problem is derived by picking a solution ζf (with associated
streamfunction ψf) for which the advection term J(ψf, ζf) is non-zero. Next, a
forcing Ff is derived according to Ff = ∂ζf/∂t + J(ψf, ζf) - nuΔζf. One solution
to the vorticity equation forced by this Ff is then ζf. (This solution may not
be realized, at least at long times, if it is unstable.)
"""
function testnonlinearterms(dt, stepper; n=128, L=2π, nu=1e-2, nnu=1, mu=0.0, nmu=0, message=false)
  n, L  = 128, 2π
  nu, nnu = 1e-2, 1
  mu, nmu = 0.0, 0
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
  function calcF!(Fh, sol, t, s, v, p, g)
    Fh .= Ffh
    nothing
  end

  prob = TwoDTurb.ForcedProblem(nx=n, Lx=L, nu=nu, nnu=nnu, mu=mu, nmu=nmu, dt=dt, stepper=stepper, calcF=calcF!)
  s, v, p, g, eq, ts = prob.state, prob.vars, prob.params, prob.grid, prob.eqn, prob.ts
  TwoDTurb.set_q!(prob, qf)

  # Step forward
  stepforward!(prob, round(Int, nt))
  TwoDTurb.updatevars!(prob)
  isapprox(v.q, qf, rtol=1e-13)
end

function testenergyenstrophy()
  nx, Lx  = 128, 2π
  ny, Ly  = 128, 3π
  g  = TwoDGrid(nx, Lx, ny, Ly)
  k0 = g.k[2] # fundamental wavenumber
  l0 = g.l[2] # fundamental wavenumber
  x, y = g.X, g.Y

  psi0 = @. sin(2*k0*x)*cos(2*l0*y) + 2sin(k0*x)*cos(3*l0*y)
    q0 = @. -((2*k0)^2+(2*l0)^2)*sin(2*k0*x)*cos(2*l0*y) - (k0^2+(3*l0)^2)*2sin(k0*x)*cos(3*l0*y)

  prob = TwoDTurb.InitialValueProblem(nx=nx, Lx=Lx, ny=ny, Ly=Ly, stepper="ForwardEuler")
  s, v, p, g, eq, ts = prob.state, prob.vars, prob.params, prob.grid, prob.eqn, prob.ts
  TwoDTurb.set_q!(prob, q0)
  TwoDTurb.updatevars!(prob)

  energyq0 = FourierFlows.TwoDTurb.energy(prob)
  enstrophyq0 = FourierFlows.TwoDTurb.enstrophy(prob)

  isapprox(energyq0, 29.0/9, rtol=1e-13) && isapprox(enstrophyq0, 2701.0/162, rtol=1e-13)
end



# Run the tests
@test testnonlinearterms(0.0005, "ForwardEuler")
@test lambdipoletest(256, 1e-3)
@test stochasticforcingbudgetstest()
@test testenergyenstrophy()
