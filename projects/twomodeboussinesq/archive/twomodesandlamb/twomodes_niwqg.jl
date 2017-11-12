include("/Users/glwagner/Numerics/FourierFlows/src/fourierflows.jl")

using FourierFlows,
      PyPlot

import FourierFlows.TwoModeBoussinesq,
       FourierFlows.NIW-QG

include("./twomodeutils.jl")


# Physical parameters
Lx  = 2*pi*100e3              # Domain extent
f0  = 1e-4                    # Inertial or Coriolis frequency
nu0 = nu1 = 1e8               # Hyperviscosities
nnu = 4                       # Hyperviscous order
N0  = 5e-3                    # Buoyancy frequency
m   = 2*pi/325                # Vertical wavenumber
eta = N0^2/(m^2*f0)

# Initial condition
Uw  = 5e-2                    # Wave speed
Ue  = 5e-2                    # Eddy speed
R   = Lx/10                   # Eddy radius
ke  = 2*pi/R                  # Inverse eddy scale
te  = 1/(Ue*ke)               # Eddy turn-over time
tf  = 2*pi/f0                 # Inertial period

# Numerical params
nx  = 256                     # Resolution
dt  = 1e-2 * tf               # Time-step
nsteps = Int(ceil(20*te/dt))  # Total number of time-steps
nsubs  = 100 #ceil(Int, 4*tf/dt)   # Number of steps between plots
nplots = ceil(Int, nsteps/nsubs)   # Number of plots


# Initialize the two problems
g  = TwoDGrid(nx, Lx)

pr_bo = TwoModeBoussinesq.Params(nu0, nnu, nu1, nnu, f0, N0, m, -Ue, 0.0)
vs_bo = TwoModeBoussinesq.Vars(g)
eq_bo = TwoModeBoussinesq.Equation(pr, g)
ts_bo = ETDRK4TimeStepper(dt, eq.LCc, eq.LCr)

pr_qg = NIW-QG.Params(nu0, nnu, nu, nnu, eta, f0, -Ue, 0.0)
vs_qg = NIW-QG.Vars(g)
eq_qg = NIW-QG.Equation(pr, g)
ts_qg = NIW-QG.Timestepper(dt, eq)




# Make initial condition
Z0 = FourierFlows.lambdipole(Ue, R, g; center=(0.0, 0.0))
u0 = Uw/2 * ones(Complex{Float64}, g.nx, g.ny)
v0 = -im*Uw/2 * ones(Complex{Float64}, g.nx, g.ny)
phi0 = (u0 + v0)*sqrt(2)

TwoModeBoussinesq.set_zeta!(vs, pr, g, Z0)
NIW-QG.set_q!(vs, pr, g, Z0)

TwoModeBoussinesq.set_uvp!(vs, pr, g, u0, v0, 0.0*u0)
NIW-QG.set_phi!(vs, pr, g, phi0)


# Save some initial properties
Z00 = maximum(Z0)
E0 = TwoModeBoussinesq.calc_energy(vs, pr, g)


# Plot
fig, axs = subplots(nrows=1, ncols=2, sharex=true, sharey=true,
  figsize=(12, 5))

dfrac = 2 # domain fraction for plot
dom = [-g.Lx/(2*d), g.Lx/(2*d)]

speed_qg(vs_qg, pr_qg, g) = wavespeed(vs_qg)
speed_bo(vs_bo, pr_bo, g) = wavespeed(vs_bo)
getqgapv(vs_qg, pr_qg, g) = vs_qg.q

pl = ComparisonPlot(
  g, vs_bo, pr_bo, vs_bo, pr_bo, 
  TwoModeBoussinesq.calc_apv, speed_bo, getqgapv, speed_qg, 
  R,
  [-Z00*2, Z00*2], [0.0, 4.0*Uw],
  dom, dom,
  L"Boussinesq q/f_0", L"Boussinesq \sqrt{u^2+v^2}", 
  L"NIW-QG q/f_0", L"NIW-QG \sqrt{u^2+v^2}", 
  "RdBu_r", "YlGnBu_r",
  L"x/R", L"y/R"
)

makeplot!(axs, pl)

for i = 1:nplots

  @time stepforward!(vs_qg, nsubs, ts_qg, eq_qg, pr_qg, g)
  @time stepforward!(vs_bo, nsubs, ts_bo, eq_bo, pr_bo, g)

  TwoModeBoussinesq.updatevars!(vs_qg, pr_qg, g)
  NIWQG.updatevars!(vs_bo, pr_bo, g)

  q = TwoModeBoussinesq.calc_apv(vs, pr, g)                                        
  E_bo = TwoModeBoussinesq.calc_energy(vs, pr, g)                                     
  s_bo = TwoModeBoussinesq.calc_speed1(vs) 
  s_qg = abs.(vs_qg.phi)

  @printf("\n
    step: %04d, t: %.3f, QG CFL: %.3f, Bo CFL: %.3f
    max QG Ro: %.2e, max QG Ro: %.2e, max QG speed: %.2e, max Bo speed: %.2e, 
    Bo E: %.6f\n\n", 
    ts.r.step, vs.t/te, 
    maximum(vs_qg.U)*ts.r.dt/g.dx, maximum(vs_bo.U)*ts.r.dt/g.dx, 
    maximum(vs_qg.q)/pr.f, maximum(q)/pr.f, 
    maximum(s_bo), maximum(s_qg), 
    E_bo/E0
  )



  makeplot!(axs, pl)

end
