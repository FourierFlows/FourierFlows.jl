using PyPlot, FourierFlows
import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy, dissipation, work, drag

# This tests the nonlinear terms in the twodturb module; it works as follows
# assume a solution ψ_guess for the unforced case. Then we insert the solution
# in the equation and compute the remainder analytically, e.g.:
# remainder = -(∂ζ_guess/∂t + J(ψ_guess, ζ_guess) - νΔζ_guess)
# Finally we use remainder as forcing. Then ζ_guess should be a solution of
# the forced problem.

n, L  = 128, 2π
ν, nν = 1e-2, 1
μ, nμ = 0.0, 0
dt, tf = 0.005, 2
nt = round(Int, tf/dt)
ns = 4

gr  = TwoDGrid(n, L)

x, y = gr.X, gr.Y

psi_guess = sin.(2x).*cos.(2y) + 2*sin.(x).*cos.(3y)
q_guess = -8*sin.(2x).*cos.(2y) - 20*sin.(x).*cos.(3y)
remainder = -ν*(64*sin.(2x).*cos.(2y) + 200*sin.(x).*cos.(3y))-8* (cos.(x).*cos.(3y).*sin.(2x).*sin.(2y)- 3cos.(2x).*cos.(2y).*sin.(x).*sin.(3y))

remainderh = rfft(remainder);
# Forcing
function calcF!(Fh, sol, t, s, v, p, g)
   Fh .= remainderh
  nothing
end

prob = TwoDTurb.ForcedProblem(nx=n, Lx=L, ν=ν, nν=nν, μ=μ, nμ=nμ, dt=dt,
  stepper="RK4", calcF=calcF!)

s, v, p, g, eq, ts = prob.state, prob.vars, prob.params, prob.grid, prob.eqn, prob.ts;

TwoDTurb.set_q!(prob, q_guess)
E = Diagnostic(energy,      prob, nsteps=nt)
D = Diagnostic(dissipation, prob, nsteps=nt)
R = Diagnostic(drag,        prob, nsteps=nt)
W = Diagnostic(work,        prob, nsteps=nt)
diags = [E, D, W, R]


function makeplot(prob, diags)

  TwoDTurb.updatevars!(prob)

  E, D, W, R = diags

  x, y = prob.grid.X, prob.grid.Y

  t = round(prob.state.t, 2)
  sca(axs[1]); cla()
  pcolormesh(x, y, prob.vars.q - q_guess )
  xlabel(L"$x$")
  ylabel(L"$y$")
  title("\$\\nabla^2\\psi(x,y,\\mu t= $t )\$")
  axis("square")
  clim(-1,1)

  sca(axs[3]); cla()

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

  # If the Ito interpretation was used for the work
  # then we need to add the drift term: I[ii2] + σ
  plot(E.time[ii], W[ii2], label=L"work ($W$)")   # Ito
  plot(E.time[ii], -R[ii], label=L"drag ($D=2\mu E$)")
  plot(E.time[ii], 0*E.time[ii], "k:", linewidth=0.5)
  ylabel("Energy sources and sinks")
  xlabel(L"$t$")
  legend(fontsize=10)

  sca(axs[2]); cla()
  plot(E.time[ii], total[ii], label=L"computed $W-D$")
  plot(E.time[ii], dEdt, "--k", label=L"numerical $dE/dt$")
  ylabel(L"$dE/dt$")
  xlabel(L"$t$")
  legend(fontsize=10)

  sca(axs[4]); cla()
  plot(E.time[ii], residual, "c-", label=L"residual $dE/dt$ = computed $-$ numerical")
  xlabel(L"$\mu t$")
  legend(fontsize=10)

  residual
end

fig, axs = subplots(ncols=2, nrows=2, figsize=(12, 8))

TwoDTurb.updatevars!(prob)
# saveoutput(out)

cfl = prob.ts.dt*maximum(
  [maximum(v.V)/g.dx, maximum(v.U)/g.dy])

res = makeplot(prob, diags)
sca(axs[1])
clim(-1,1)
pause(0.01)


# Step forward
for i = 1:ns
  tic()

  stepforward!(prob, diags, round(Int, nt/ns))
  tc = toq()

  TwoDTurb.updatevars!(prob)
  # saveoutput(out)

  cfl = prob.ts.dt*maximum(
    [maximum(v.V)/g.dx, maximum(v.U)/g.dy])

  res = makeplot(prob, diags)
  pause(0.01)

  @printf("step: %04d, t: %.1f, cfl: %.3f, time: %.2f s\n",
    prob.step, prob.t, cfl, tc)

  # savename = @sprintf("./plots/stochastictest_kf%d_%06d.png", kf, prob.step)
  # savefig(savename, dpi=240)
end

# savename = @sprintf("./plots/stochastictest_kf%d_%06d.png", kf, prob.step)
# savefig(savename, dpi=240)

# savediagnostic(E, "energy", out.filename)
# savediagnostic(D, "dissipation", out.filename)
# savediagnostic(W, "work", out.filename)
# savediagnostic(R, "drag", out.filename)
