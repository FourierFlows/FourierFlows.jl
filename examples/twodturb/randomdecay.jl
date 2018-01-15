using PyPlot, FourierFlows
import FourierFlows.TwoDTurb

 n = 128
 L = 2π
 ν = 1e-4  # Laplacian viscosity
nν = 1
dt = 1e0   # Time step
nt = 100   # Number of time steps

prob = TwoDTurb.InitialValueProblem(nx=n, Lx=L, ν=ν, nν=nν, dt=dt,
  stepper="RK4")
TwoDTurb.set_q!(prob, rand(n, n))

# Step forward
fig = figure(); tic()
for i = 1:10
  stepforward!(prob, nt)
  TwoDTurb.updatevars!(prob)  

  cfl = maximum(prob.vars.U)*prob.grid.dx/prob.ts.dt
  @printf("step: %04d, t: %6.1f, cfl: %.2f, ", prob.step, prob.t, cfl)
  toc(); tic()

  clf(); imshow(prob.vars.q); pause(0.01)
end

show()
