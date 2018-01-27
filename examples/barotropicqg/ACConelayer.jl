using FourierFlows, PyPlot, JLD2


import FourierFlows.BarotropicQG
import FourierFlows.BarotropicQG: energy, energy00, enstrophy, enstrophy00

nx  = 512
ν  = 8.0e-10
νn = 2
f0 = -1.0

β = 1.4015
Lx = 2π
μ = 1.0e-2
F = 0.0012

FU(t) = F

η(x, y) = 2*cos.(10x).*cos.(10y)


g  = BarotropicQG.Grid(nx, Lx)
p  = BarotropicQG.Params(g, f0, β, FU, η, μ, ν, νn)
v  = BarotropicQG.Vars(g)
eq = BarotropicQG.Equation(p, g)


# Time-stepping
dt  = 2e-2
nsteps = 20000
nsubs  = 500


# ts = FilteredETDRK4TimeStepper(dt, eq.LC, g)
ts = ETDRK4TimeStepper(dt, eq.LC)
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
E00 = Diagnostic(energy00, prob; nsteps=nsteps)
Q00 = Diagnostic(enstrophy00, prob; nsteps=nsteps)
diags = [E, E00, Q, Q00]

# Create Output
get_sol(prob) = prob.state.sol # extracts the Fourier-transformed solution
get_u(prob) = irfft(im*g.Lr.*g.invKKrsq.*prob.state.sol, g.nx)
out = Output(prob, filename, (:sol, get_sol), (:u, get_u))



function plot_output(prob, fig, axs; drawcolorbar=false)
  # Plot the vorticity field and the evolution of energy and enstrophy.

  s, v, p, g = prob.state, prob.vars, prob.params, prob.grid
  BarotropicQG.updatevars!(prob)

  sca(axs[1])
  pcolormesh(g.X, g.Y, v.q)
  axis("square")
  xlim(0, 2)
  ylim(0, 2)
  # axis("off")
  if drawcolorbar==true
    colorbar()
  end

  sca(axs[2])
  cla()
  plot(μ*E.time[1:E.prob.step], E.data[1:prob.step], label=L"$E_{\psi}$")
  plot(μ*E.time[1:E00.prob.step], E00.data[1:prob.step], label=L"$E_U$")
  xlabel(L"\mu t")
  ylabel(L"E")
  legend()

  sca(axs[3])
  cla()
  plot(μ*Q.time[1:Q.prob.step], Q.data[1:prob.step], label=L"$Q_{\psi}$")
  plot(μ*Q00.time[1:Q00.prob.step], Q00.data[1:prob.step], label=L"$Q_U$")
  xlabel(L"\mu t")
  ylabel(L"Q")
  legend()
  pause(0.001)
end


fig, axs = subplots(ncols=3, nrows=1, figsize=(15, 4))
plot_output(prob, fig, axs; drawcolorbar=true)


# Step forward
startwalltime = time()

while prob.step < nsteps
  stepforward!(prob, diags, nsubs)

  # Message
  log = @sprintf("step: %04d, t: %d, E: %.4f, Q: %.4f, τ: %.2f min",
    prob.step, prob.t, E.value, Q.value,
    (time()-startwalltime)/60)

  println(log)

  plot_output(prob, fig, axs; drawcolorbar=false)

end
println((time()-startwalltime))

plot_output(prob, fig, axs; drawcolorbar=false)

savename = @sprintf("%s_%09d.png", joinpath(plotpath, plotname), prob.step)
savefig(savename, dpi=240)
