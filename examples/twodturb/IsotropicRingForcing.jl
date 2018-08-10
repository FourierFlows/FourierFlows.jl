using PyPlot, FourierFlows, Random
import Random:rand
import Printf.@printf
import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy, dissipation, work, drag

  n, L  = 256, 2π
nu, nnu = 1e-7, 2
mu, nmu = 1e-1, 0
dt, tf = 0.005, 0.2/mu
nt = round(Int, tf/dt)
ns = 4

# Forcing
kf, dkf = 12.0, 2.0
σ = 0.1

gr  = TwoDGrid(n, L)

force2k = exp.(-(sqrt.(gr.KKrsq) .- kf).^2/(2*dkf^2))
@. force2k[gr.KKrsq < 2.0^2 ] = 0
@. force2k[gr.KKrsq > 20.0^2 ] = 0
@. force2k[gr.Kr < 2π/L] = 0
σ0 = FourierFlows.parsevalsum(force2k.*gr.invKKrsq/2.0, gr)/(gr.Lx*gr.Ly)
force2k .= σ/σ0 * force2k

Random.seed!(1234)

function calcF!(F, sol, t, s, v, p, g)
  eta = exp.(2π*im*rand(eltype(L), size(sol)))/sqrt(s.dt)
  eta[1, 1] = 0
  @. F = eta .* sqrt(force2k)
  nothing
end

prob = TwoDTurb.ForcedProblem(nx=n, Lx=L, nu=nu, nnu=nnu, mu=mu, nmu=nmu, dt=dt, stepper="RK4", calcF=calcF!)
s, v, p, g, eq, ts = prob.state, prob.vars, prob.params, prob.grid, prob.eqn, prob.ts;

TwoDTurb.set_q!(prob, 0*g.X)
E = Diagnostic(energy,      prob, nsteps=nt)
D = Diagnostic(dissipation, prob, nsteps=nt)
R = Diagnostic(drag,        prob, nsteps=nt)
W = Diagnostic(work,        prob, nsteps=nt)
diags = [E, D, W, R]


function makeplot(prob, diags)
  TwoDTurb.updatevars!(prob)
  E, D, W, R = diags

  t = round(mu*prob.state.t; digits=2)
  sca(axs[1]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, prob.vars.q)
  xlabel(L"$x$")
  ylabel(L"$y$")
  title("\$\\nabla^2\\psi(x,y,\\mu t= $t )\$")
  axis("square")

  sca(axs[3]); cla()

  i₀ = 1
  ii = (i₀):E.count-1
  ii2 = (i₀+1):E.count
  dEdt = (E[ii2] .- E[ii])/prob.ts.dt
  # dEdt = W - D - R?

  # If the Ito interpretation was used for the work
  # then we need to add the drift term
  # total = W[ii2]+σ - D[ii] - R[ii]      # Ito
  total = W[ii2] .- D[ii] .- R[ii]        # Stratonovich
  residual = dEdt .- total

  # If the Ito interpretation was used for the work
  # then we need to add the drift term: I[ii2] + σ
  plot(mu*E.time[ii], W[ii2], label=L"work ($W$)")   # Ito
  # plot(mu*E.time[ii], W[ii2] , label=L"work ($W$)")      # Stratonovich
  plot(mu*E.time[ii], σ .+ 0*E.time[ii], "--", label=L"ensemble mean  work ($\langle W\rangle $)")
  # plot(mu*E.time[ii], -D[ii], label="dissipation (\$D\$)")
  plot(mu*E.time[ii], -R[ii], label=L"drag ($D=2\mu E$)")
  plot(mu*E.time[ii], 0*E.time[ii], "k:", linewidth=0.5)
  ylabel("Energy sources and sinks")
  xlabel(L"$\mu t$")
  legend(fontsize=10)

  sca(axs[2]); cla()
  plot(mu*E.time[ii], total[ii], label=L"computed $W-D$")
  plot(mu*E.time[ii], dEdt, "--k", label=L"numerical $dE/dt$")
  ylabel(L"$dE/dt$")
  xlabel(L"$\mu t$")
  legend(fontsize=10)

  sca(axs[4]); cla()
  plot(mu*E.time[ii], residual, "c-", label=L"residual $dE/dt$ = computed $-$ numerical")
  xlabel(L"$\mu t$")
  legend(fontsize=10)

  residual
end

fig, axs = subplots(ncols=2, nrows=2, figsize=(12, 8))

# Step forward
for i = 1:ns

  time = @elapsed stepforward!(prob, diags, round(Int, nt/ns))

  TwoDTurb.updatevars!(prob)
  # saveoutput(out)

  cfl = prob.ts.dt*maximum([maximum(v.U)/g.dx, maximum(v.V)/g.dy])
  res = makeplot(prob, diags)
  pause(0.01)

  @printf("step: %04d, t: %.1f, cfl: %.3f, time: %.2f s\n", prob.step, prob.t, cfl, time)

  # savename = @sprintf("./plots/stochastictest_kf%d_%06d.png", kf, prob.step)
  # savefig(savename, dpi=240)
end

# savename = @sprintf("./plots/stochastictest_kf%d_%06d.png", kf, prob.step)
# savefig(savename, dpi=240)

# savediagnostic(E, "energy", out.filename)
# savediagnostic(D, "dissipation", out.filename)
# savediagnostic(W, "work", out.filename)
# savediagnostic(R, "drag", out.filename)
