import FourierFlows.VerticallyCosineBoussinesq

cfl(prob) = maximum([maximum(abs.(prob.vars.U)), maximum(abs.(prob.vars.V))]*
              prob.ts.dt/prob.grid.dx)

function test_lambdipole(n, dt; L=2π, Ue=1, Re=L/20, nu0=0, nnu0=1,
  ti=L/Ue*0.01, nm=3, message=false, atol=1e-2)

  nt = round(Int, ti/dt)

  prob = VerticallyCosineBoussinesq.Problem(nx=n, Lx=L,
    nu0=nu0, nnu0=nnu0, dt=dt, stepper="FilteredRK4")
  x, y, Z = prob.grid.X, prob.grid.Y, prob.vars.Z # nicknames

  Z0 = FourierFlows.lambdipole(Ue, Re, prob.grid)
  VerticallyCosineBoussinesq.set_Z!(prob, Z0)

  xZ = zeros(nm)   # centroid of abs(Z)
  Ue_m = zeros(nm) # measured dipole speed

  # Step forward
  for i = 1:nm
    stepforward!(prob, nt)
    VerticallyCosineBoussinesq.updatevars!(prob)
    xZ[i] = mean(abs.(Z).*x) / mean(abs.(Z))

    if i > 1
      Ue_m[i] = (xZ[i]-xZ[i-1]) / ((nt-1)*dt)
    else
      Ue_m[i] = 0.0
    end

    if message
      println("     step: %04d, t: %3.3f, time: %.3f, cfl: %.2f\n",
        prob.step, prob.t, cfl(prob))
    end
  end

  isapprox(Ue, mean(Ue_m[2:end]), atol=atol)
end

function testnonlinearterms(dt, stepper; n=128, L=2π, nu=1e-2, nnu=1, mu=0.0, nmu=0, message=false)
  n, L  = 128, 2π
  nu, nnu = 1e-2, 1
  mu, nmu = 0.0, 0
  tf = 1.0
  nt = round(Int, tf/dt)

  gr  = TwoDGrid(n, L)
  x, y = gr.X, gr.Y
  psif = @. sin(2x)*cos(2y) + 2sin(x)*cos(3y)
  Zf = @. -8sin(2x)*cos(2y) - 20sin(x)*cos(3y)

  Ff = @. -(
    nu*( 64sin(2x)*cos(2y) + 200sin(x)*cos(3y) )
    + 8*( cos(x)*cos(3y)*sin(2x)*sin(2y) - 3cos(2x)*cos(2y)*sin(x)*sin(3y) )
  )

  Ffh = rfft(Ff)

  # Forcing
  function calcF!(Fh, sol, t, s, v, p, g)
    Fh[:, :, 1] .= Ffh
    nothing
  end

  prob = VerticallyCosineBoussinesq.Problem(
    nx=n, Lx=L, nu0=nu, nnu0=nnu, mu0=mu, nmu0=nmu, dt=dt, stepper=stepper, calcF=calcF!)
  s, v, p, g, eq, ts = prob.state, prob.vars, prob.params, prob.grid, prob.eqn, prob.ts
  VerticallyCosineBoussinesq.set_Z!(prob, Zf)

  # Step forward
  stepforward!(prob, round(Int, nt))
  VerticallyCosineBoussinesq.updatevars!(prob)
  isapprox(v.Z, Zf, rtol=1e-13)
end


e1(u, v, p, m, N) = @. ( u^2 + v^2 + m^2*p^2/N^2 )/2
e1(prob) = e1(prob.vars.u, prob.vars.v, prob.vars.p, prob.params.m, prob.params.N)
wavecentroid(prob) = (FourierFlows.xmoment(e1(prob), prob.grid), FourierFlows.ymoment(e1(prob), prob.grid))

function test_groupvelocity(kw; n=128, L=2π, f=1.0, N=1.0, m=4.0, uw=1e-2, rtol=1e-3)
   σ = f*sqrt(1 + (N*kw/m)^2)
  tσ = 2π/σ
  dt = tσ/100
  nt = round(Int, 3tσ/(2dt)) # 3/2 a wave period
  cga = N^2*kw/(σ*m^2) # analytical group velocity
  δ = L/10 # envelopt scale

  prob = VerticallyCosineBoussinesq.Problem(f=f, N=N, m=m, nx=n, Lx=L, dt=dt, stepper="FilteredRK4")
  envelope(x, y) = exp(-x^2/(2*δ^2))
  VerticallyCosineBoussinesq.set_planewave!(prob, uw, kw; envelope=envelope)

  t₋₁ = prob.t
  xw₋₁, yw₋₁ = wavecentroid(prob)

  stepforward!(prob, nt)
  VerticallyCosineBoussinesq.updatevars!(prob)
  xw, yw = wavecentroid(prob)
  cgn = (xw-xw₋₁) / (prob.t-t₋₁)
  isapprox(cga, cgn, rtol=rtol)
end

@test testnonlinearterms(0.0005, "ForwardEuler")
@test test_lambdipole(256, 1e-3)
@test test_groupvelocity(16)
