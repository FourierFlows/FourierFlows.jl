using FourierFlows, PyPlot, JLD2


import FourierFlows.BarotropicQG
import FourierFlows.BarotropicQG: energy, enstrophy

nx  = 256
ν  = 0.0e-10
νn = 4
f0 = -1.0

β = 1.4015
Lx = 2π
μ = 1.0e-2
F = 0.0012

FU(t) = F

η(x, y) = 2*cos.(10*x).*cos.(10*y)


g  = BarotropicQG.Grid(nx, Lx)
p  = BarotropicQG.Params(g, f0, β, FU, η, μ, ν, νn)
v  = BarotropicQG.Vars(g)
eq = BarotropicQG.Equation(p, g)

# eta = 2*cos.(10*g.X).*cos.(10*g.Y)
# v = cos(10*g.X).*cos(10*g.Y) + cos(10*g.X-10*g.Y)/2
#
# mean(eta.*v)
#
# etah = fft(eta)
# vh = fft(v)
# println(sum(conj(vh).*etah).re / (g.nx^2.0*g.ny^2.0))
#
# pcolormesh(g.X, g.Y, eta)



# Time-stepping
dt  = 1e-2
nsteps = 200000
nsubs  = 1000


# ts = ETDRK4TimeStepper(dt, eq.LC)
ts = FilteredETDRK4TimeStepper(dt, eq.LC, g)
prob = FourierFlows.Problem(g, v, p, eq, ts)
s = prob.state


# Files
filepath = "."
plotpath = "./plots"
plotname = "testplots"
filename = joinpath(filepath, "testdata.jld2")

# File management
if isfile(filename); rm(filename); end
if !isdir(plotpath); mkdir(plotpath); end

# Initialize with zeros
BarotropicQG.set_zeta!(prob, 0*g.X)


# Create Diagnostic -- "energy" and "enstrophy" are functions imported
# at the top.
E = Diagnostic(energy, prob; nsteps=nsteps)
Q = Diagnostic(enstrophy, prob; nsteps=nsteps)
diags = [E, Q]

# Create Output
get_sol(prob) = prob.vars.sol # extracts the Fourier-transformed solution
get_u(prob) = irfft(im*g.Lr.*g.invKKrsq.*prob.vars.sol, g.nx)
out = Output(prob, filename, (:sol, get_sol), (:u, get_u))

# Step forward
startwalltime = time()

while prob.step < nsteps
  stepforward!(prob, diags; nsteps=nsubs)

  # Message
  log = @sprintf("step: %04d, t: %d, E: %.4f, Q: %.4f, τ: %.2f min",
    prob.step, prob.t, E.value, Q.value,
    (time()-startwalltime)/60)

  println(log)


  BarotropicQG.updatevars!(prob)
  fig, axs = subplots(ncols=2, nrows=1, figsize=(12, 4))

  axes(axs[1])
  pcolormesh(g.X, g.Y, prob.vars.zeta)
  axis("equal")
  xlim(0, 2)
  ylim(0, 2)
  colorbar()
  # clim(-40, 40)
  axs[1][:axis]("off")

  axes(axs[2])
  plot(mu*E.time[1:E.prob.step], E.data[1:prob.step])
  plot(mu*Q.time[1:E.prob.step], Q.data[1:prob.step])
  xlabel(L"\mu t")
  ylabel(L"E, \, Z")

end

BarotropicQG.updatevars!(prob)
fig, axs = subplots(ncols=2, nrows=1, figsize=(12, 4))

axes(axs[1])
pcolormesh(g.X, g.Y, prob.vars.zeta)
axis("equal")
colorbar()
# clim(-40, 40)
axs[1][:axis]("off")

axes(axs[2])
plot(mu*E.time[1:E.prob.step], E.data[1:prob.step])
plot(mu*Q.time[1:E.prob.step], Q.data[1:prob.step])
xlabel(L"\mu t")
ylabel(L"E, \, Z")

savename = @sprintf("%s_%09d.png", joinpath(plotpath, plotname), prob.step)
savefig(savename, dpi=240)
