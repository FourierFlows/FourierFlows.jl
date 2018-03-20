using FourierFlows, FourierFlows.KuramotoSivashinsky, PyPlot

nx = 256
Lx = 32π
dt = 0.01
nt = 20000
prob = InitialValueProblem(nx=nx, Lx=Lx, dt=dt, stepper="ETDRK4")

x = prob.grid.x
u0 = @. cos(x/16) * (1 + sin(x/16))
set_u!(prob, u0)

sub = 10
T = [ dt*j for i=1:nx, j=sub:sub:nt ]
X = [ x[i] for i=1:nx, j=sub:sub:nt ]

u = zeros(nx, round(Int, nt/sub))

count = 1
for i = 1:nt
  println(i)
  stepforward!(prob)
  updatevars!(prob)
  if i % sub == 0
    u[:, count] .= prob.vars.u
    count += 1
  end
end

close("all")
fig, ax = subplots()
pcolormesh(X, T, u)
show()
