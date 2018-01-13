include("../../src/fourierflows.jl")

using FourierFlows, PyPlot, JLD2

import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy

# Physical parameters
 n = 128
 L = 2π
 nν = 2
  ν = 0e-8

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
prob = TwoDTurb.InitialValueProblem(nx=n, Lx=L, ν=ν, nν=nν, dt=dt, 
  withfilter=true)
g = prob.grid


# Initial condition closely following pyqg barotropic example
# that reproduces the results of the paper by McWilliams (1984)
srand(1234)
k0, E0 = 6, 0.5
modk = sqrt.(g.KKrsq)
psik = zeros(g.nk, g.nl)
psik =  (modk.^2 .* (1 + (modk/k0).^4)).^(-0.5)
psik[1, 1] = 0.0
psih = (randn(g.nkr, g.nl)+im*randn(g.nkr, g.nl)).*psik
psih = psih.*prob.ts.filter
Ein = real(sum(g.KKrsq.*abs2.(psih)/(g.nx*g.ny)^2))
psih = psih*sqrt(E0/Ein)
qi = -irfft(g.KKrsq.*psih, g.nx)
E0 = FourierFlows.parsevalsum(g.KKrsq.*abs2.(psih), g)

TwoDTurb.set_q!(prob, qi)

# Create Diagnostic -- "energy" is a function imported at the top.
E = Diagnostic(energy, prob; nsteps=nsteps)
Z = Diagnostic(enstrophy, prob; nsteps=nsteps)
diags = [E, Z]

# Create Output
get_sol(prob) = prob.vars.sol # extracts the Fourier-transformed solution
get_u(prob) = irfft(im*g.Lr.*g.invKKrsq.*prob.vars.sol, g.nx)
out = Output(prob, filename, (:sol, get_sol), (:u, get_u))

# Step forward
startwalltime = time()

while prob.step < nsteps
  stepforward!(prob, diags, nsteps)

  # Message
  log = @sprintf("step: %04d, t: %d, ΔE: %.4f, ΔZ: %.4f, τ: %.2f min",
    prob.step, prob.t, E.value/E.data[1], Z.value/Z.data[1],
    (time()-startwalltime)/60)

  println(log)
end

TwoDTurb.updatevars!(prob)
fig, axs = subplots(ncols=2, nrows=1, figsize=(12, 4))

axes(axs[1])
pcolormesh(g.X, g.Y, prob.vars.q)
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
