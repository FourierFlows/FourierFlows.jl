include("../../src/fourierflows.jl")

using FourierFlows, PyPlot, JLD2

import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy

# Physical parameters
 n = 256
 L = 2π
 nν = 4
  ν = 0e-8
withfilter = true

# Time-stepping
dt = 5e-3
nsteps = 8000
nsubs  = 200

# Files
filepath = "."
plotpath = "./plots"
plotname = "testplots"
filename = joinpath(filepath, "testdata.jld2")

# File management
if isfile(filename); rm(filename); end
if !isdir(plotpath); mkdir(plotpath); end

# Initialize with random numbers
prob = TwoDTurb.InitialValueProblem(n, L, ν, nν, dt, withfilter)


# Initial condition closely following pyqg barotropic example
# that reproduces the results of the paper by McWilliams (1984)

k0, E0 = 6, 0.5
modk = sqrt.(prob.grid.KKrsq)
psik = zeros(prob.grid.nk, prob.grid.nl)
psik =  (modk.^2 .* (1 + (modk/k0).^4)).^(-0.5)
psik[1, 1] = 0.0
psih = (randn(prob.grid.nkr, prob.grid.nl)+im*randn(prob.grid.nkr, prob.grid.nl)).*psik
psih = psih.*prob.ts.filter
Ein = real(sum(prob.grid.KKrsq.*abs2.(psih)/(prob.grid.nx*prob.grid.ny)^2))
psih = psih*sqrt(E0/Ein)
qi = -irfft(prob.grid.KKrsq.*psih, prob.grid.nx)
E0 = FourierFlows.parsevalsum(prob.grid.KKrsq.*abs2.(psih), prob.grid)


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

function get_u(prob) # extracts the physics-space x-velocity.
  irfft(im*prob.grid.Lr.*prob.grid.invKKrsq.*prob.vars.sol, prob.grid.nx)
end

# One way to create an Output type.
out = Output(prob, filename, (:sol, get_sol), (:u, get_u))


# Step forward
startwalltime = time()

fig, axs = subplots(ncols=2, nrows=1, figsize=(12, 4))

while prob.step < nsteps

  # axes(axs[1])
  # pcolormesh(abs.(prob.vars.sol))
  # # colorbar()
  # show()
  # draw()

  # Step forward
  stepforward!(prob, diags; nsteps=nsubs)

  # Message
  log = @sprintf("step: %04d, t: %d, ΔE: %.4f, ΔZ: %.4f, τ: %.2f min",
    prob.step, prob.t, E.value/E.data[1], Z.value/Z.data[1],
    (time()-startwalltime)/60)

  println(log)

  # axes(axs[2])
  # pcolormesh(abs.(prob.vars.sol))
  # # colorbar()
  # show()
  # draw()


  # Plot and save output
  if prob.step/nsubs % 1 == 0.0

    # TwoDTurb.updatevars!(prob)
    # saveoutput(out) # saves output to out.filename
    # figure(1)
    # pcolormesh(prob.vars.q)
    # show()
    # draw()

  end

end

TwoDTurb.updatevars!(prob)
fig, axs = subplots(ncols=2, nrows=1, figsize=(12, 4))

figure(1)
axes(axs[1])
pcolormesh(prob.grid.x, prob.grid.y, prob.vars.q)
axis("equal")
colorbar()
clim(-40, 40)
axs[1][:axis]("off")

axes(axs[2])
plot(E.time[1:E.prob.step], E.data[1:prob.step]/E.data[1])
plot(Z.time[1:E.prob.step], Z.data[1:prob.step]/Z.data[1])
xlabel(L"t")
ylabel(L"\Delta E, \, \Delta Z")

savename = @sprintf("%s_%09d.png", joinpath(plotpath, plotname), prob.step)
savefig(savename, dpi=240)
