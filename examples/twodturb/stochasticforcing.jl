using PyPlot, FourierFlows
import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy, dissipation, injection, drag

 n, L  =  256, 2π
 ν, nν = 1e-3,  1
 μ, nμ = 1e-1, -1
dt, tf = 2e-3, 1000

nt = round(Int, tf/dt)
ns = 100

# Forcing
fi = 1.0
ki = 16
amplitude = fi*ki/sqrt(dt) * n^2/4

function calcF!(F, sol, t, s, v, p, g)
  if t == s.t # not a substep
    F .= 0.0

    θk = 2π*rand() 
    phase = 2π*im*rand()

    i₁ = round(Int, abs(ki*cos(θk))) + 1
    j₁ = round(Int, abs(ki*sin(θk))) + 1  # j₁ >= 1
    j₂ = g.nl + 2 - j₁                    # e.g. j₁ = 1 => j₂ = nl+1

    if j₁ != 1  # apply forcing to l = (+/-)l★ mode
      F[i₁, j₁] = amplitude*exp(phase)
      F[i₁, j₂] = amplitude*exp(phase)
    else        # apply forcing to l=0 mode
      F[i₁, 1] = 2amplitude*exp(phase)
    end
  end

  nothing
end

prob = TwoDTurb.ForcedProblem(nx=n, Lx=L, ν=ν, nν=nν, μ=μ, nμ=nμ, dt=dt, 
  calcF=calcF!, stepper="RK4")

E = Diagnostic(energy,      prob, nsteps=nt) 
D = Diagnostic(dissipation, prob, nsteps=nt)
R = Diagnostic(drag,        prob, nsteps=nt)
I = Diagnostic(injection,   prob, nsteps=nt)
diags = [E, D, I, R]

filename = @sprintf("stochastictest_ki%d.jld2", ki)
getsol(prob) = deepcopy(prob.state.sol)
out = Output(prob, filename, (:sol, getsol))

function makeplot(prob, diags)
  E, D, I, R = diags

  TwoDTurb.updatevars!(prob)  

  close("all")
  E, D, I, R = diags
  fig, axs = subplots(ncols=3, nrows=1, figsize=(13, 4))

  sca(axs[1]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, prob.vars.q)
  xlabel(L"x")
  ylabel(L"y")

  sca(axs[2]); cla()

  i₀ = 1
  dEdt = (E[(i₀+1):E.count] - E[i₀:E.count-1])/prob.ts.dt
  ii = (i₀+1):E.count

  # dEdt = I - D - R?
  total = I[ii] - D[ii] - R[ii]
  residual = dEdt - total

  plot(E.time[ii], I[ii], "o", markersize=0.5, label="injection (\$I\$)")
  plot(E.time[ii], -D[ii], label="dissipation (\$D\$)")
  plot(E.time[ii], -R[ii], label="drag (\$R\$)")
  plot(E.time[ii], residual, "c-", label="residual")

  ylabel("Energy sources and sinks")
  xlabel(L"t")
  legend(fontsize=10)

  sca(axs[3]); cla()
  plot(E.time[ii], E[ii])
  xlabel(L"t")
  ylabel(L"E")

  tight_layout()

  residual
end

# Step forward
for i = 1:ns
  tic()
  stepforward!(prob, diags, round(Int, nt/ns))
  tc = toq()

  TwoDTurb.updatevars!(prob)  
  saveoutput(out)

  cfl = prob.ts.dt*maximum(
    [maximum(prob.vars.V)/prob.grid.dx, maximum(prob.vars.U)/prob.grid.dy])

  res = makeplot(prob, diags)
  pause(0.1)

  @printf("step: %04d, t: %.1f, cfl: %.3f, time: %.2f s, mean(res) = %.3e\n", 
    prob.step, prob.t, cfl, tc, mean(res))
    
  savename = @sprintf("./plots/stochastictest_ki%d_%06d.png", ki, prob.step)
  savefig(savename, dpi=240)
end

savediagnostic(E, "energy", out.filename)
savediagnostic(D, "dissipation", out.filename)
savediagnostic(I, "injection", out.filename)
savediagnostic(R, "drag", out.filename)
