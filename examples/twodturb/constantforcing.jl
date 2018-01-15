using PyPlot, FourierFlows
import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, dissipation, injection, drag

 n, L  =  128, 2π
 ν, nν = 2e-4,  1
 μ, nμ = 1e-1, -1
dt, tf = 2e-1, 1000

nt = round(Int, tf/dt)

# Forcing
fi = 0.01
ki = 8
θk = π/4

i₁ = round(Int, abs(ki*cos(θk))) + 1
j₁ = round(Int, abs(ki*sin(θk))) + 1  # j₁ >= 1
j₂ = n + 2 - j₁                       # e.g. j₁ = 1 => j₂ = nl+1
amplitude = fi*ki * n^2/4

function calcF!(F, sol, t, s, v, p, g)
  F[i₁, j₁] = amplitude
  F[i₁, j₂] = amplitude
  nothing
end

prob = TwoDTurb.ForcedProblem(nx=n, Lx=L, ν=ν, nν=nν, μ=μ, nμ=nμ, dt=dt, 
  calcF=calcF!, stepper="RK4")

function runtest(prob, nt)
  TwoDTurb.set_q!(prob, rand(prob.grid.nx, prob.grid.ny))

  E = Diagnostic(energy,      prob, nsteps=nt)
  D = Diagnostic(dissipation, prob, nsteps=nt)
  R = Diagnostic(drag,        prob, nsteps=nt)
  I = Diagnostic(injection,   prob, nsteps=nt)
  diags = [E, D, I, R]

  tic()
  stepforward!(prob, diags, round(Int, nt))
  @printf("step: %04d, t: %.1f, time: %.2f s\n", prob.step, prob.t, toq())

  diags
end

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

  plot(E.time[ii], I[ii], label="injection (\$I\$)")
  plot(E.time[ii], -D[ii], label="dissipation (\$D\$)")
  plot(E.time[ii], -R[ii], label="drag (\$R\$)")
  plot(E.time[ii], total, label=L"I-D-R")
  plot(E.time[ii], dEdt, "k:", label=L"E_t")
  plot(E.time[ii], residual, "c-", label="residual")

  ylabel("Energy sources and sinks")
  xlabel(L"t")
  legend(fontsize=10)

  sca(axs[3]); cla()
  plot(E.time[ii], E[ii])
  xlabel(L"t")
  ylabel(L"E")

  tight_layout()
  show()

  residual
end

diags = runtest(prob, nt)
res = makeplot(prob, diags)

@printf "rms residual: %.3e" sqrt(mean(res.^2))
