using PyPlot, FourierFlows
import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy, dissipation, injection

  n = 128
  L = 2π
  μ = 0.0    # Bottom drag
  ν = 1e-3   # Laplacian viscosity
 nν = 1
 dt = 1e-2   # Time step
 tf = 10000  # final time
ndp = 5000   # Timesteps between plots

  fi = 0.5
  ki = 16
norm = n^2/4

function calcF!(F, sol, t, s, v, p, g)
  if t == s.t # not a substep
    F .= 0.0
    θ, ξ = 2π*rand(2) 
    i₁ = round(Int, abs(ki*cos(θ))) + 1
    j₁ = round(Int, abs(ki*sin(θ))) + 1 # j₁ >= 1
    j₂ = n + 2 - j₁ # e.g. j₁ = 1 => j₂ = nl+1

    if j₁ != 1
      F[i₁, j₁] = fi*ki/sqrt(s.dt)*norm*exp(im*ξ)
      F[i₁, j₂] = fi*ki/sqrt(s.dt)*norm*exp(im*ξ)
    else
      F[i₁, j₁] = fi*ki/sqrt(s.dt)*2norm*exp(im*ξ)
    end
  end

  nothing
end

prob = TwoDTurb.ForcedProblem(nx=n, Lx=L, ν=ν, nν=nν, μ=μ, dt=dt, 
  calcF=calcF!, 
  stepper="FilteredRK4")

E = Diagnostic(energy, prob, nsteps=round(Int, tf/dt))
Z = Diagnostic(enstrophy, prob, nsteps=round(Int, tf/dt))
D = Diagnostic(dissipation, prob, nsteps=round(Int, tf/dt))
I = Diagnostic(injection, prob, nsteps=round(Int, tf/dt))
diags = [E, Z, D, I]

# Step forward
fig, axs = subplots(ncols=3, figsize=(12, 4))
for i = 1:round(Int, tf/dt/ndp)

  tic()
  stepforward!(prob, diags, ndp)
  TwoDTurb.updatevars!(prob)  

  cfl = prob.ts.dt*maximum(
    [maximum(prob.vars.V)/prob.grid.dx, maximum(prob.vars.U)/prob.grid.dy])
  @printf("step: %04d, t: %.1f, cfl: %.2f, time: %.3f s\n", prob.step, prob.t, 
    cfl, toq())

  sca(axs[1]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, prob.vars.q)

  sca(axs[2]); cla()

  dEdt = (E[3:E.count] - E[1:E.count-2])/(2*dt)
  ii = 2:(E.count-1)

  plot(E.time[ii], dEdt)
  plot(E.time[ii], -D[ii])
  plot(E.time[ii], -μ*E[ii])
  plot(E.time[ii], I[ii], "k.", markersize=0.1)
  xlabel(L"t")

  sca(axs[3]); cla()
  plot(E.time[ii], E[ii])
  xlabel(L"t")
  ylabel(L"E")

  axs[1][:tick_params](bottom=false, labelbottom=false, 
    left=false, labelleft=false)

  pause(0.01)
end

show()
