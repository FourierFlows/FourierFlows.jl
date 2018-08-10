using FourierFlows, PyPlot, JLD2, FFTW, Random, Statistics
import Printf.@sprintf
import Statistics: mean
import FourierFlows.BarotropicQG
import FourierFlows.BarotropicQG: energy, enstrophy

# Numerical parameters and time-stepping parameters
nx  = 256      # 2D resolution = nx^2
stepper = "FilteredETDRK4"   # timestepper
dt  = 0.02     # timestep
nsteps = 8000  # total number of time-steps
nsubs  = 500   # number of time-steps for plotting
               # (nsteps must be multiple of nsubs)

# Physical parameters
Lx  = 2π       # domain size
nu  = 0e-05    # viscosity
nnu = 1        # viscosity order
beta = 15.0    # planetary PV gradient
mu   = 0e-1    # bottom drag

# Initialize problem
prob = BarotropicQG.InitialValueProblem(nx=nx, Lx=Lx, beta=beta, nu=nu,
                                        nnu=nnu, mu=mu, dt=dt, stepper=stepper)
s, v, p, g, eq, ts = prob.state, prob.vars, prob.params, prob.grid, prob.eqn, prob.ts;


# Files
filepath = "."
plotpath = "./plots"
plotname = "snapshot"
filename = joinpath(filepath, "decayingbetaturb.jld2")

# File management
if isfile(filename); rm(filename); end
if !isdir(plotpath); mkdir(plotpath); end




# Initial condition that has power only at wavenumbers with
# 8<L/(2π)*sqrt(kx^2+ky^2)<10 and initial energy E0
Random.seed!(1234)
E0 = 0.1
modk = ones(g.nkr, g.nl)
@. modk[real(g.KKrsq) < (8*2*pi/g.Lx)^2] = 0
@. modk[real(g.KKrsq) > (10*2*pi/g.Lx)^2] = 0
@. modk[1, :] = 0
psih = (randn(g.nkr, g.nl) + im*randn(g.nkr, g.nl)).*modk
psih = psih.*prob.ts.filter
Ein = real(sum(g.KKrsq.*abs2.(psih)/(g.nx*g.ny)^2))
psih = psih*sqrt(E0/Ein)
qi = -irfft(g.KKrsq.*psih, g.nx)
E0 = FourierFlows.parsevalsum(g.KKrsq.*abs2.(psih), g)

BarotropicQG.set_zeta!(prob, qi)

# Create Diagnostic -- "energy" and "enstrophy" are functions imported at the top.
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
  xticks(-2:2:2)
  yticks(-2:2:2)
  title(L"vorticity $\zeta = \partial_x v - \partial_y u$")
  if drawcolorbar==true
    colorbar()
  end

  sca(axs[2])
  cla()
  pcolormesh(g.X, g.Y, v.psi)
  axis("square")
  xticks(-2:2:2)
  yticks(-2:2:2)
  title(L"streamfunction $\psi$")
  if drawcolorbar==true
    colorbar()
  end

  sca(axs[3])
  cla()
  plot(mean(v.zeta, dims=1).', g.y.')
  plot(0*g.Y[1,:], g.Y[1,:], "k--")
  ylim(-Lx/2, Lx/2)
  xlim(-2, 2)
  title(L"zonal mean $\zeta$")

  sca(axs[4])
  cla()
  plot(mean(v.u, dims=1).', g.y.')
  plot(0*g.Y[1,:], g.Y[1,:], "k--")
  ylim(-Lx/2, Lx/2)
  xlim(-0.5, 0.5)
  title(L"zonal mean $u$")

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
