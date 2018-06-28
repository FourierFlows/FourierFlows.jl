using FourierFlows, PyPlot, JLD2


import FourierFlows.BarotropicQG
import FourierFlows.BarotropicQG: energy, energy00, enstrophy, enstrophy00

# Numerical parameters and time-stepping parameters
nx  = 512      # 2D resolution = nx^2
stepper = "FilteredETDRK4"   # timestepper
dt  = 2e-2     # timestep
nsteps = 20000 # total number of time-steps
nsubs  = 500   # number of time-steps for plotting
               # (nsteps must be multiple of nsubs)

# Physical parameters
Lx  = 2π       # domain size
nu  = 8.0e-10  # viscosity
nnu = 2        # viscosity order
f0  = -1.0     # Coriolis parameter
beta = 1.4015  # the y-gradient of planetary PV
mu   = 1.0e-2  # linear drag
   F = 0.0012  # normalized wind stress forcing on domain-averaged
               # zonal flow U(t) flow

# Topographic PV
topoPV(x, y) = 2*cos.(10x).*cos.(10y)

# Forcing on the domain-averaged U equation
calcFU(t) = F


# Initialize problem
prob = BarotropicQG.ForcedProblem(nx=nx, Lx=Lx, f0=f0, beta=beta, eta=topoPV,
                  calcFU=calcFU, nu=nu, nnu=nnu, mu=mu, dt=dt, stepper=stepper)
s, v, p, g, eq, ts = prob.state, prob.vars, prob.params, prob.grid, prob.eqn, prob.ts;


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

  # Plot the PV field and the evolution of energy and enstrophy.

  s, v, p, g = prob.state, prob.vars, prob.params, prob.grid
  BarotropicQG.updatevars!(prob)

  sca(axs[1])
  pcolormesh(g.X, g.Y, v.q)
  axis("square")
  xlim(0, 2)
  xticks(0:0.5:2)
  ylim(0, 2)
  yticks(0:0.5:2)
  title(L"$\nabla^2\psi + \eta$ (part of the domain)")
  if drawcolorbar==true
    colorbar()
  end

  sca(axs[2])
  cla()
  plot(mu*E.time[1:E.prob.step], E.data[1:prob.step], label=L"$E_{\psi}$")
  plot(mu*E.time[1:E00.prob.step], E00.data[1:prob.step], label=L"$E_U$")

  xlabel(L"\mu t")
  ylabel(L"E")
  legend()

  sca(axs[3])
  cla()
  plot(mu*Q.time[1:Q.prob.step], Q.data[1:prob.step], label=L"$Q_{\psi}$")
  plot(mu*Q00.time[1:Q00.prob.step], Q00.data[1:prob.step], label=L"$Q_U$")
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
