using FourierFlows, FourierFlows.KuramotoSivashinsky, PyPlot

nx = 512
Lx = 32π
dt = 0.1
nt = 2000
prob = InitialValueProblem(nx=nx, Lx=Lx, dt=dt, stepper="ETDRK4")

x = prob.grid.x
u0 = @. cos(x/16) * (1 + sin(x/16))
set_u!(prob, u0)

T = [ dt*j for i=1:nx, j=1:nt ]
X = [ x[i] for i=1:nx, j=1:nt ]
u = zeros(nx, nt)

for i = 1:nt
  stepforward!(prob)
  updatevars!(prob)
  u[:, i] .= prob.vars.u
end

close("all")
fig, ax = subplots(figsize=(6, 8))
pcolormesh(X/2π, T, u)
ax[:tick_params](bottom=false, left=false, labelbottom=false, labelleft=false)
show()
