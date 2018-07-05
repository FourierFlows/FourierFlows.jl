using FourierFlows, PyPlot, JLD2

import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy

# Parameters
  n = 128
  L = 2π
nnu = 2
 nu = 0e-8
 dt = 5e-3
nsteps = 8000
nsubs = 200

# Files
filepath = "."
plotpath = "./plots"
plotname = "testplots"
filename = joinpath(filepath, "testdata.jld2")

# File management
if isfile(filename); rm(filename); end
if !isdir(plotpath); mkdir(plotpath); end

# Initialize problem
prob = TwoDTurb.Problem(; nx=n, Lx=L, ny=n, Ly=L, nu=nu, nnu=nnu, dt=dt, stepper="FilteredRK4")
g = prob.grid

# Initial condition closely following pyqg barotropic example
# that reproduces the results of the paper by McWilliams (1984)
srand(1234)
k0, E0 = 6, 0.5
qi = FourierFlows.peakedisotropicspectrum(g, k0, E0, mask = prob.ts.filter)
TwoDTurb.set_q!(prob, qi)

# Create Diagnostic -- "energy" is a function imported at the top.
E = Diagnostic(energy, prob; nsteps=nsteps)
Z = Diagnostic(enstrophy, prob; nsteps=nsteps)
diags = [E, Z] # A list of Diagnostics types passed to "stepforward!" will
# be updated every timestep. They should be efficient to calculate and
# have a small memory footprint. (For example, the domain-integrated kinetic
# energy is just a single number for each timestep). See the file in
# src/diagnostics.jl and the stepforward! function in timesteppers.jl.

# Create Output
get_sol(prob) = prob.vars.sol # extracts the Fourier-transformed solution
get_u(prob) = irfft(im*g.Lr.*g.invKKrsq.*prob.vars.sol, g.nx)
out = Output(prob, filename, (:sol, get_sol), (:u, get_u))


function plot_output(prob, fig, axs; drawcolorbar=false)
  # Plot the vorticity field and the evolution of energy and enstrophy.
  TwoDTurb.updatevars!(prob)
  sca(axs[1])
  pcolormesh(prob.grid.X, prob.grid.Y, prob.vars.q)
  clim(-40, 40)
  axis("off")
  axis("square")
  if drawcolorbar==true
    colorbar()
  end

  sca(axs[2])
  cla()
  plot(E.time[1:E.prob.step], E.data[1:prob.step]/E.data[1])
  plot(Z.time[1:Z.prob.step], Z.data[1:prob.step]/Z.data[1])
  xlabel(L"t")
  ylabel(L"\Delta E, \, \Delta Z")

  pause(0.01)
end

# Step forward
startwalltime = time()

fig, axs = subplots(ncols=2, nrows=1, figsize=(12, 4))
plot_output(prob, fig, axs; drawcolorbar=true)

while prob.step < nsteps
  stepforward!(prob, diags, nsubs)

  # Message
  log = @sprintf("step: %04d, t: %d, ΔE: %.4f, ΔZ: %.4f, τ: %.2f min",
    prob.step, prob.t, E.value/E.data[1], Z.value/Z.data[1], (time()-startwalltime)/60)

  println(log)
  plot_output(prob, fig, axs; drawcolorbar=false)
end

plot_output(prob, fig, axs; drawcolorbar=true)

savename = @sprintf("%s_%09d.png", joinpath(plotpath, plotname), prob.step)
savefig(savename, dpi=240)
