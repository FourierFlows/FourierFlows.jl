import FourierFlows, FourierFlows.VerticallyFourierBoussinesq

cfl(prob) = maximum([maximum(abs.(prob.vars.U)), maximum(abs.(prob.vars.V))]*
              prob.ts.dt/prob.grid.dx)

e1(u, v, p, m, N) = @. abs2(u) + abs2(v) + m^2*abs2(p)/N^2
e1(prob) = e1(prob.vars.u, prob.vars.v, prob.vars.p, prob.params.m, prob.params.N)
wavecentroid(prob) = (FourierFlows.xmoment(e1(prob), prob.grid), FourierFlows.ymoment(e1(prob), prob.grid))

function lambdipoletest(n, dt; L=2π, Ue=1, Re=L/20, nu0=0, nnu0=1,
  ti=L/Ue*0.01, nm=3, message=false)

  nt = round(Int, ti/dt)

  prob = VerticallyFourierBoussinesq.Problem(nx=n, Lx=L, nu0=nu0, nnu0=nnu0, dt=dt, stepper="FilteredRK4")
  x, y, Z = prob.grid.X, prob.grid.Y, prob.vars.Z # nicknames

  Z0 = FourierFlows.lambdipole(Ue, Re, prob.grid)
  VerticallyFourierBoussinesq.set_Z!(prob, Z0)

  xZ = zeros(nm)   # centroid of abs(Z)
  Ue_m = zeros(nm) # measured dipole speed

  # Step forward
  for i = 1:nm
    stepforward!(prob, nt)
    VerticallyFourierBoussinesq.updatevars!(prob)
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

  isapprox(Ue, mean(Ue_m[2:end]), atol=0.02)
end

import FourierFlows.VerticallyFourierBoussinesq: mode1u

function test_groupvelocity(nkw; n=128, L=2π, f=1.0, N=1.0, m=4.0, uw=1e-2, rtol=1e-3, del=L/10)
  kw = nkw*2π/L
   σ = f*sqrt(1 + (N*kw/m)^2)
  tσ = 2π/σ
  dt = tσ/100
  nt = round(Int, 3tσ/(2dt)) # 3/2 a wave period
  cga = N^2*kw/(σ*m^2) # analytical group velocity

  prob = VerticallyFourierBoussinesq.Problem(f=f, N=N, m=m, nx=n, Lx=L, dt=dt, stepper="FilteredRK4")
  envelope(x, y) = exp(-x^2/(2*del^2))
  VerticallyFourierBoussinesq.set_planewave!(prob, uw, nkw; envelope=envelope)

  t₋₁ = prob.t
  xw₋₁, yw₋₁ = wavecentroid(prob)

  stepforward!(prob, nt)
  VerticallyFourierBoussinesq.updatevars!(prob)
  xw, yw = wavecentroid(prob)
  cgn = (xw-xw₋₁) / (prob.t-t₋₁)

  # close("all")
  # fig, ax = subplots()
  # imshow(mode1u(prob))

  isapprox(cga, cgn, rtol=rtol)
end

@test lambdipoletest(256, 1e-3)
@test test_groupvelocity(16)
