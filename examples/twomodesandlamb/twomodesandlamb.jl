include("/Users/glwagner/Numerics/FourierFlows/src/fourierflows.jl")

using FourierFlows,
      PyPlot

import FourierFlows.TwoModeBoussinesq

include("./twomodeutils.jl")


# Physical parameters
Lx  = 2*pi*100e3              # Domain extent
f0  = 1e-4                    # Inertial or Coriolis frequency
nu0 = nu1 = 1e8               # Hyperviscosities
nnu = 4                       # Hyperviscous order
N0  = 5e-3                    # Buoyancy frequency
m   = 2*pi/750                # Vertical wavenumber

# Initial condition
Uw  = 5e-1                    # Wave speed
Ue  = 5e-2                    # Eddy speed
R   = Lx/10                   # Eddy radius
ke  = 2*pi/R                  # Inverse eddy scale
te  = 1/(Ue*ke)               # Eddy turn-over time
tf  = 2*pi/f0                 # Inertial period

# Numerical params
nx  = 128                     # Resolution
dt  = 1e-2 * tf               # Time-step
nsteps = Int(ceil(40*te/dt))  # Total number of time-steps
nsubs  = 100 #ceil(Int, 4*tf/dt)   # Number of steps between plots
nplots = ceil(Int, nsteps/nsubs)   # Number of plots


# Initialize problem
g  = TwoDGrid(nx, Lx)
pr = TwoModeBoussinesq.Params(nu0, nnu, nu1, nnu, f0, N0, m, -Ue, 0.0)
vs = TwoModeBoussinesq.Vars(g)
eq = TwoModeBoussinesq.Equation(pr, g)
ts = ETDRK4TimeStepper(dt, eq.LCc, eq.LCr)


# Make initial condition
Z0 = FourierFlows.lambdipole(Ue, R, g; center=(0.0, 0.0))
u0 = Uw/2 * ones(Complex{Float64}, g.nx, g.ny)
v0 = -im*Uw/2 * ones(Complex{Float64}, g.nx, g.ny)

TwoModeBoussinesq.set_zeta!(vs, pr, g, Z0)
TwoModeBoussinesq.set_uvp!(vs, pr, g, u0, v0, 0.0*u0)


# Save some initial properties
Z00 = maximum(Z0)
E0 = TwoModeBoussinesq.calc_energy(vs, pr, g)


# Plot
fig, axs = subplots(nrows=1, ncols=2, sharex=true, sharey=true,
  figsize=(12, 5))

d = 0.5 # domain fraction for plot
speed(vs, pr, g) = wavespeed(vs)
pl = TwoComponentPlot(
  g, vs, pr, TwoModeBoussinesq.calc_apv, speed, R,
  [-Z00, Z00], [0.0, 2.0*Uw],
  [-0.5*d*g.Lx, 0.5*d*g.Lx], [-0.5*d*g.Ly, 0.5*d*g.Ly],
  L"q/f_0", L"\sqrt{u^2+v^2}", 
  "RdBu_r", "YlGnBu_r",
  L"x/R", L"y/R"
)

makeplot!(axs, pl)
#niwqgplot(axs, vs, pr, g, Z00*te, 2*Uw, R, te, E0)

for i = 1:nplots

  @time stepforward!(vs, nsubs, ts, eq, pr, g)

  TwoModeBoussinesq.updatevars!(vs, pr, g)

  q = TwoModeBoussinesq.calc_apv(vs, pr, g)
  E = TwoModeBoussinesq.calc_energy(vs, pr, g)
  s = TwoModeBoussinesq.calc_speed1(vs) 

  @printf("
    step: %04d, t: %.3f, CFL: %.3f,
    max Ro: %.2e, max speed: %.2e, E: %.6f\n\n", 
    ts.r.step, vs.t/te, 
    maximum([abs.(2*vs.u); abs2.(2*vs.v); vs.U; vs.V])*ts.r.dt/g.dx, 
    maximum(q)/pr.f, maximum(s), E/E0
  )

  makeplot!(axs, pl)

end
