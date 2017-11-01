__precompile__()

include("../../src/fourierflows.jl")

using FourierFlows, FourierFlows.NIWQG

struct LambDipole
  R           # eddy radius
  Ue          # dipole strength
  Uw          # initial wave speed
  N0          # buoyancy frequency
  f0          # inertial frequency
  m           # wave vertical wavenumber
  kape        # QG biharmonic viscosity
  nuw         # NIW biharmonic viscosity
  Lx          # domain size
end


""" Get parameters for Rocha's lamb dipole simulation, per Rocha, Wagner, and 
Young (submitted JFM 2017). """
function getrochadipole(;
  Lx   = 2*pi*200e3,      # domain size
  R    = Lx/15,           # eddy radius
  Ue   = 5e-2,            # dipole strength
  Uw   = 5e-1,            # initial wave speed
  N0   = 5e-3,            # buoyancy frequency
  f0   = 1e-4,            # inertial frequency
  m    = 2*pi/325,        # wave vertical wavenumber
  kape = 5e7,             # QG biharmonic viscosity
  nuw  = 1e7              # NIW biharmonic viscosity
  )

  LambDipole(R, Ue, Uw, N0, f0, m, kape, nuw, Lx)
end



""" Construct the 'wave-Lamb-dipole' interaction problem. """
function makedipoleproblem(;nx=nx, dtfrac=0.0025,
  Lx   = 2*pi*200e3,      # domain size
  R    = Lx/15,           # eddy radius
  Ue   = 5e-2,            # dipole strength
  Uw   = 5e-1,            # initial wave speed
  N0   = 5e-3,            # buoyancy frequency
  f0   = 1e-4,            # inertial frequency
  m    = 2*pi/325,        # wave vertical wavenumber
  kape = 5e7,             # QG biharmonic viscosity
  nuw  = 1e7              # NIW biharmonic viscosity
  )

  ke = 2π/R
  te = 1.0/(Ue*ke)

  dt = dtfrac * te
  kw = N0*m/f0
  eta = N0^2.0/(f0*m^2.0)

  prob = NIWQG.NIWQGProblem(nx, Lx, dt, kape, 4, nuw, 4, eta, f0, 
    -Ue, 0.0)

  # Initial condition
  q0   = FourierFlows.lambdipole(Ue, R, prob.grid; center=(0.0, 0.0))
  phi0 = ((1.0+im)/sqrt(2)*Uw
    * ones(Complex{Float64}, prob.grid.nx, prob.grid.ny))

  NIWQG.set_q!(prob, q0)
  NIWQG.set_phi!(prob, phi0)

  prob
end
