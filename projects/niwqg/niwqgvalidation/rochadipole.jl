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
function get_rocha_dipole(;
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
function make_dipole_problem(;nx=nx, dtfrac=0.1,
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

  dt = dtfrac * 2*pi/f0
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

  return prob
end



#=
""" Construct the 'wave-Lamb-dipole' interaction problem. """
function make_dipole_problem(;nx=nx, wavescaling=1.0, dtfrac=0.1)

  di = get_rocha_dipole()

  dt = dtfrac * 2*pi/di.f0
  kw = di.N0*di.m/di.f0
  eta = 0.5*di.f0*(di.N0*di.m/di.f0)^2.0

  prob = NIWQG.NIWQGProblem(nx, di.Lx, dt, di.kape, 4, di.nuw, 4, eta, di.f0, 
    0.0, 0.0)

  # Initial condition
  q0   = FourierFlows.lambdipole(di.Ue, di.R, prob.grid; center=(0.0, 0.0))
  phi0 = ((1.0+im)/sqrt(2)*di.Uw*wavescaling
    * ones(Complex{Float64}, prob.grid.nx, prob.grid.ny))

  NIWQG.set_q!(prob, q0)
  NIWQG.set_phi!(prob, phi0)

  return prob
end
=#
