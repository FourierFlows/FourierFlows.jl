import FourierFlows.TwoDTurb

cfl(prob) = maximum([maximum(abs.(prob.vars.U)), maximum(abs.(prob.vars.V))]*
              prob.ts.dt/prob.grid.dx)

function lambdipoletest(n, dt; L=2π, Ue=1, Re=L/20, ν=0, nν=1, ti=L/Ue*0.01,
  nm=3, message=false)

  nt = round(Int, ti/dt)

  prob = TwoDTurb.InitialValueProblem(nx=n, Lx=L, ν=ν, nν=nν, dt=dt, stepper="FilteredRK4")
  x, y, q = prob.grid.X, prob.grid.Y, prob.vars.q # nicknames

  q0 = FourierFlows.lambdipole(Ue, Re, prob.grid)
  TwoDTurb.set_q!(prob, q0)

  xq = zeros(nm)   # centroid of abs(q)
  Ue_m = zeros(nm) # measured dipole speed

  # Step forward
  for i = 1:nm
    tic()
    stepforward!(prob, nt)
    TwoDTurb.updatevars!(prob)
    xq[i] = mean(abs.(q).*x) / mean(abs.(q))

    if i > 1
      Ue_m[i] = (xq[i]-xq[i-1]) / ((nt-1)*dt)
    end

    if message
      @printf("     step: %04d, t: %3.3f, time: %.3f, cfl: %.2f\n",
        prob.step, prob.t, toq(), cfl(prob))
    end
  end
  # println(Ue_m)
  # println(abs(Ue - mean(Ue_m[2:end]))/abs(Ue))
  isapprox(Ue, mean(Ue_m[2:end]), rtol=1e-2)
end

@test lambdipoletest(256, 1e-3)
