import FourierFlows.KuramotoSivashinsky

function growingwavetest()
  nx, Lx = 256, 4π
  dt, nt = 1e-6, 100
  prob = KuramotoSivashinsky.Problem(nx=nx, Lx=Lx, dt=dt, stepper="FilteredETDRK4")
  x = prob.grid.x

  # Initial condition
  a = 1e-4
  u0 = @. a*cos(x/2)
  ua(x, t) = a*exp(3t/16)*cos(x/2) + a^2*2/3*(exp(3t/8)-1)*sin(x)
  KuramotoSivashinsky.set_u!(prob, u0)

  stepforward!(prob, nt)
  KuramotoSivashinsky.updatevars!(prob)
  u, t = prob.vars.u, prob.t

  arbitraryerror = 1e-14 # what is an appropriate error?
  L1(u) = prob.grid.dx*sum(abs.(u))/Lx

  L1(u-ua.(x, t))/L1(ua.(x, 0)) < arbitraryerror # tests relative error
end

@test growingwavetest()
