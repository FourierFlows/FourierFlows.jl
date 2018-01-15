using PyPlot, FourierFlows
import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy

  n = 128
  L = 2π
  μ = 0.0    # Bottom drag
  ν = 1e-6   # Laplacian viscosity
 nν = 2
 dt = 2e-3   # Time step
 tf = 1000   # final time
ndp = 10000  # Timesteps between plots

  fi = 1.0
  ki = 16
norm = n^2/4

function calcF!(F, sol, t, v, p, g)
  θ = 2π*rand() 
  ξ = 2π*rand() 
  i₁ = round(Int, abs(ki*cos(θ))) + 1
  j₁ = round(Int, abs(ki*sin(θ))) + 1 # j₁ >= 1
  j₂ = n + 2 - j₁ # e.g. j₁ = 1 => j₂ = nl+1

  if j₁ != 1
    F[i₁, j₁] = fi*ki * norm*exp(im*ξ)
    F[i₁, j₂] = fi*ki * norm*exp(im*ξ)
  else
    F[i₁, j₁] = fi*ki * 2norm*exp(im*ξ)
  end

  nothing
end

prob = TwoDTurb.ForcedProblem(nx=n, Lx=L, ν=ν, nν=nν, dt=dt, calcF=calcF!)
#TwoDTurb.set_q!(prob, rand(n, n))

E = Diagnostic(energy, prob, nsteps=round(Int, tf/dt))
Z = Diagnostic(enstrophy, prob, nsteps=round(Int, tf/dt))
diags = [E, Z]

# Step forward
fig, axs = subplots(ncols=2, figsize=(8, 4))
tic()
for i = 1:round(Int, tf/dt/ndp)
  stepforward!(prob, diags, nsteps=ndp)
  TwoDTurb.updatevars!(prob)  

  cfl = maximum(prob.vars.U)*prob.grid.dx/prob.ts.dt
  @printf("step: %04d, t: %6.1f, cfl: %.2f, ", prob.step, prob.t, cfl)
  toc(); tic()

  sca(axs[1]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, prob.vars.q)

  sca(axs[2]); cla()
  plot(E.time[1:E.count], E.data[1:E.count])

  pause(0.01)

end

show()
