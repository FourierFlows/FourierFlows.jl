using FourierFlows, FourierFlows.KuramotoSivashinsky, PyPlot

nx, Lx = 256, 4π
dt, nt = 1e-6, 100
prob = InitialValueProblem(nx=nx, Lx=Lx, dt=dt, stepper="ETDRK4")
x = prob.grid.x

a = 1e-4 # amplitude of initial condition
u0 = @. a*cos(x/2) # initial condition
set_u!(prob, u0) 

ua(x, t) = a*exp(3t/16)*cos(x/2) + a^2*2/3*(exp(3t/8)-1)*sin(x) # analytical solution

t = dt*(1:nt)
T = [ t[j] for i=1:nx, j=1:nt ]
X = [ x[i] for i=1:nx, j=1:nt ]
U = zeros(nx, nt)
Ua = ua.(X, T)

for i = 1:nt
  stepforward!(prob)
  updatevars!(prob)
  U[:, i] .= prob.vars.u
end

L1(u) = prob.grid.dx*squeeze(sum(abs.(u), 1), 1)/Lx
L1(u::Array{T,1}) where {T} = prob.grid.dx*sum(abs.(u))/Lx

close("all")
fig, ax = subplots()
plot(t, L1(U-Ua)/L1(ua.(x, 0)))
ylabel(L"\mathrm{L1} = \frac{1}{L}\int |u-u_a| \, \mathrm{d} x")
xlabel(L"t")
show()
