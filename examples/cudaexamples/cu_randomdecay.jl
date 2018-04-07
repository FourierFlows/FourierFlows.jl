using CuArrays, PyPlot
using FourierFlows
import FourierFlows.TwoDTurb

   n, L = 512, 2π   # Domain
nu, nnu = 1e-6, 1  # Viscosity
 dt, nt = 1, 100   # Time step

q0 = CuArray(rand(n, n))

prob = TwoDTurb.CuProblem(nx=n, Lx=L, nu=nu, nnu=nnu, dt=dt, stepper="FilteredRK4")
TwoDTurb.set_q!(prob, q0)

# Step forward
fig = figure(); tic()
for i = 1:10
  stepforward!(prob, nt)
  TwoDTurb.updatevars!(prob)  

  cfl = maximum(prob.vars.U)*prob.grid.dx/prob.ts.dt
  @printf("step: %04d, t: %6.1f, cfl: %.2f, ", prob.step, prob.t, cfl)
  toc(); tic()

  q = Array(prob.vars.q)
  clf(); imshow(q); pause(0.01)
end

# Transfer grid and vars to device.
g = TwoDGrid(prob.grid)
U = Array(prob.vars.U)
V = Array(prob.vars.V)
q = Array(prob.vars.q)

# Calculate radial energy spectrum
E = @. 0.5*(U^2 + V^2) # energy density
Eh = rfft(E)
kr, Ehr = FourierFlows.radialspectrum(Eh, g, refinement=1)

fig, axs = subplots(ncols=2, figsize=(10, 4))

sca(axs[1]); cla()
pcolormesh(g.X, g.Y, q)
xlabel(L"x")
ylabel(L"y")
title("Vorticity")

sca(axs[2]); cla()
plot(kr, abs.(Ehr))
xlabel(L"k_r")
ylabel(L"\int | \hat{E} | k_r \mathrm{d} k_{\theta}")
title("Radial energy spectrum")

xlim(0, n/4)
axs[2][:set_yscale]("log")

tight_layout()
show()
