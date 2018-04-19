using CuArrays
using PyPlot, FourierFlows

import FourierFlows.TwoDTurb

# -----
# Setup
# -----

  n = 256
  L = 2π
 nu = 1e-6
nnu = 1
 dt = 1
 nt = 100

prob = TwoDTurb.CuProblem(nx=n, Lx=L, nu=nu, nnu=nnu, dt=dt, stepper="FilteredRK4")
TwoDTurb.set_q!(prob, rand(n, n))


# ---
# Run
# ---
cfl(prob) = maximum(prob.vars.U)*prob.grid.dx/prob.ts.dt

fig = figure(); tic()
for i = 1:10
  stepforward!(prob, nt)
  TwoDTurb.updatevars!(prob)  

  @printf("step: %04d, t: %6.1f, cfl: %.2f, ", prob.step, prob.t, cfl(prob))
  toc(); tic()

  clf()
  imshow(prob.vars.q)
  pause(0.01)
end


# ----
# Plot
# ----

 E = 0.5*(prob.vars.U.^2 + prob.vars.V.^2) # energy density
Eh = rfft(E)
kr, Ehr = FourierFlows.radialspectrum(Eh, prob.grid, refinement=1)

fig, axs = subplots(ncols=2, figsize=(10, 4))

sca(axs[1]); cla()
pcolormesh(prob.grid.X, prob.grid.Y, prob.vars.q)

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
