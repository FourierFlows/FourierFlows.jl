using FourierFlows, PyPlot, JLD2

import FourierFlows.BarotropicQG
import FourierFlows.BarotropicQG: energy, enstrophy

# Physical parameters
nx  = 256
Lx  = 2π     # domain size
nu  = 0e-05  # viscosity
nnu = 1      # viscosity order

f0   = 1.0    # Coriolis parameter
beta = 15.0    # planetary PV gradient
mu   = 0e-1    # bottom drag

# Time-stepping
stepper = "FilteredETDRK4"
dt = 0.02
nsteps = 8000
nsubs  = 500


prob = BarotropicQG.InitialValueProblem(nx=nx, Lx=Lx, f0=f0, beta=beta, nu=nu,
nnu=nnu, mu=mu, dt=dt, stepper=stepper)
s, v, p, g, eq, ts = prob.state, prob.vars, prob.params, prob.grid, prob.eqn, prob.ts;


# Files
filepath = "."
plotpath = "./plots"
plotname = "testplots"
filename = joinpath(filepath, "decayingbetaturb.jld2")

# File management
if isfile(filename); rm(filename); end
if !isdir(plotpath); mkdir(plotpath); end




# Initial condition that has power only at wavenumbers with
# 8<L/(2π)*sqrt(kx^2+ky^2)<10 and initial energy E0
srand(1234)
E0 = 0.1
modk = ones(g.nkr, g.nl)
modk[real(g.KKrsq).<(8*2*pi/g.Lx)^2]=0
modk[real(g.KKrsq).>(10*2*pi/g.Lx)^2]=0
modk[1, :]=0
psih = (randn(g.nkr, g.nl)+im*randn(g.nkr, g.nl)).*modk
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
  # Plot the vorticity and streamfunction fields as well as the zonal mean
  # vorticity and the zonal mean zonal velocity.

  s, v, p, g = prob.state, prob.vars, prob.params, prob.grid
  BarotropicQG.updatevars!(prob)

  sca(axs[1])
  cla()
  pcolormesh(g.X, g.Y, v.q)
  axis("square")
  # axis("off")
  title(L"vorticity $\zeta$")
  if drawcolorbar==true
    colorbar()
  end

  sca(axs[2])
  cla()
  pcolormesh(g.X, g.Y, v.psi)
  axis("square")
  # axis("off")
  title(L"streamfunction $\psi$")
  if drawcolorbar==true
    colorbar()
  end

  sca(axs[3])
  cla()
  plot(mean(v.zeta, 1).', g.Y[1,:])
  plot(0*g.Y[1,:], g.Y[1,:], "k--")
  ylim(-Lx/2, Lx/2)
  xlim(-2, 2)
  title(L"mean $\zeta$")

  sca(axs[4])
  cla()
  plot(mean(v.u, 1).', g.Y[1,:])
  plot(0*g.Y[1,:], g.Y[1,:], "k--")
  ylim(-Lx/2, Lx/2)
  xlim(-0.5, 0.5)
  title(L"mean $u$")

  pause(0.001)
end



fig, axs = subplots(ncols=2, nrows=2, figsize=(8, 8))
plot_output(prob, fig, axs; drawcolorbar=false)

# Step forward
startwalltime = time()

while prob.step < nsteps
  stepforward!(prob, diags, nsubs)

  # Message
  log = @sprintf("step: %04d, t: %d, E: %.4f, Q: %.4f, τ: %.2f min",
    prob.step, prob.t, E.value, Z.value,
    (time()-startwalltime)/60)

  println(log)

  plot_output(prob, fig, axs; drawcolorbar=false)
end

# how long did it take?
println((time()-startwalltime))

plot_output(prob, fig, axs; drawcolorbar=false)

# save the figure as png
savename = @sprintf("%s_%09d.png", joinpath(plotpath, plotname), prob.step)
savefig(savename)
