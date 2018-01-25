using FourierFlows, PyPlot, JLD2

import FourierFlows.BarotropicQG
import FourierFlows.BarotropicQG: energy, enstrophy

# Physical parameters
nx  = 256
ν  = 0.0
νn = 2
f0 = 1.0

β = 12.0
Lx = 1.0
μ = 1.0e-1
F = 0.0

FU(t) = F
η(x, y) = 0*x

g  = BarotropicQG.Grid(nx, Lx)
p  = BarotropicQG.Params(g, f0, β, FU, η, μ, ν, νn)
v  = BarotropicQG.Vars(g)
eq = BarotropicQG.Equation(p, g)



# Time-stepping
dt = 0.01
nsteps = 8000
nsubs  = 200


ts = FilteredETDRK4TimeStepper(dt, eq.LC)
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

BarotropicQG.set_zeta!(prob, qi)


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
  xlabel(L"\mu t")
  ylabel(L"E")
  legend()

  sca(axs[3])
  cla()
  plot(μ*Q.time[1:Q.prob.step], Q.data[1:prob.step], label=L"$Q_{\psi}$")
  xlabel(L"\mu t")
  ylabel(L"Q")
  legend()
  pause(0.001)
end



fig, axs = subplots(ncols=2, nrows=1, figsize=(12, 4))
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
