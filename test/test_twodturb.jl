import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy, dissipation, injection, drag

# -----------------------------------------------------------------------------
# TWODTURB's TEST FUNCTIONS

cfl(prob) = maximum([maximum(abs.(prob.vars.U)), maximum(abs.(prob.vars.V))]*
              prob.ts.dt/prob.grid.dx)

# LAMB DIPOLE TEST
function lambdipoletest(n, dt; L=2π, Ue=1, Re=L/20, ν=0, nν=1, ti=L/Ue*0.01,
  nm=3, message=false)

  nt = round(Int, ti/dt)

  prob = TwoDTurb.InitialValueProblem(nx=n, Lx=L, ν=ν, nν=nν, dt=dt,
                                        stepper="FilteredRK4")
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


# ISOTROPIC RING FORCING BUDGETS
function stochasticforcingbudgetstest( n ; dt = 0.01, L=2π, ν=1e-7, nν=2,
                                            μ = 1e-1, nμ = 0, message=false)

  n, L  = 256, 2π
  ν, nν = 1e-7, 2
  μ, nμ = 1e-1, 0
  dt, tf = 0.005, 0.2/μ
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
    eta[1, 1] = 0
    @. F = eta .* sqrt(force2k)
    nothing
  end

  prob = TwoDTurb.ForcedProblem(nx=n, Lx=L, ν=ν, nν=nν, μ=μ, nμ=nμ, dt=dt,
   stepper="RK4", calcF=calcF!)

  s, v, p, g, eq, ts = prob.state, prob.vars, prob.params, prob.grid, prob.eqn, prob.ts;

  TwoDTurb.set_q!(prob, 0*g.X)
  E = Diagnostic(energy,      prob, nsteps=nt)
  D = Diagnostic(dissipation, prob, nsteps=nt)
  R = Diagnostic(drag,        prob, nsteps=nt)
  I = Diagnostic(injection,   prob, nsteps=nt)
  diags = [E, D, I, R]

  # Step forward

  stepforward!(prob, diags, round(Int, nt))

  TwoDTurb.updatevars!(prob)

  cfl = prob.ts.dt*maximum(
    [maximum(v.V)/g.dx, maximum(v.U)/g.dy])

  E, D, I, R = diags

  t = round(μ*prob.state.t, 2)

  i₀ = 1
  dEdt = (E[(i₀+1):E.count] - E[i₀:E.count-1])/prob.ts.dt
  ii = (i₀):E.count-1
  ii2 = (i₀+1):E.count

  # dEdt = I - D - R?
  # If the Ito interpretation was used for the work
  # then we need to add the drift term
  # total = I[ii2]+σ - D[ii] - R[ii]      # Ito
  total = I[ii2] - D[ii] - R[ii]        # Stratonovich

  residual = dEdt - total



  if message
    @printf("step: %04d, t: %.1f, cfl: %.3f, time: %.2f s\n",
            prob.step, prob.t, cfl, tc)
  end

  # println(mean(abs.(residual)))
  isapprox(mean(abs.(residual)), 0, atol=1e-4)
end

# -----------------------------------------------------------------------------
# Running the tests


@test lambdipoletest(256, 1e-3)

@test stochasticforcingbudgetstest(256)
