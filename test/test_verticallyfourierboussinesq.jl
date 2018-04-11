import FourierFlows.VerticallyFourierBoussinesq

cfl(prob) = maximum([maximum(abs.(prob.vars.U)), maximum(abs.(prob.vars.V))]*
              prob.ts.dt/prob.grid.dx)

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
    tic()
    stepforward!(prob, nt)
    VerticallyFourierBoussinesq.updatevars!(prob)  
    xZ[i] = mean(abs.(Z).*x) / mean(abs.(Z))

    if i > 1
      Ue_m[i] = (xZ[i]-xZ[i-1]) / ((nt-1)*dt)
    else
      Ue_m[i] = 0.0
    end

    if message
      @printf("     step: %04d, t: %3.3f, time: %.3f, cfl: %.2f\n",
        prob.step, prob.t, toq(), cfl(prob))
    end
  end

  isapprox(Ue, mean(Ue_m[2:end]), atol=0.02)
end

@test lambdipoletest(256, 1e-3)
