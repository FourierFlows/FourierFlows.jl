__precompile__()

include("../fourierflows.jl")

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# B E G I N    M O D U L E    B A R O T R O P I C Q G --------------------------

module BarotropicQG

using FourierFlowTypes, Domains, TimeSteppers

export Grid, Vars, Params, Equation
export ForwardEulerTimeStepper, ETDRK4TimeStepper
export RK4TimeStepper, AB3TimeStepper

export set_q!, updatevars!, calc_NL!, calc_NL, stepforward!
export southern_gaussian_mountain

Grid = TwoDGrid

# P A R A M S  T Y P E -------------------------------------------------------- 
type Params <: AbstractParams
  f0::Float64                       # Constant inertial frequency
  beta::Float64                     # Inertial frequency y-gradient
  FU::Function                      # Time-dependent forcing of domain average flow
  eta::Array{Float64, 2}            # Topographic PV
  etah::Array{Complex{Float64}, 2}  # Transform of topographic PV
  mu::Float64                       # Linear drag
  nu::Float64                       # Vorticity viscosity
  nun::Int                          # Vorticity hyperviscous order
end

# With element-wise topographic generating function
function Params(g::Grid, f0::Float64, beta::Float64, FU::Function, eta::Function,
  mu::Float64, nu::Float64, nun::Int)

  etagridded = eta.(g.X, g.Y)
  etah = rfft(etagridded)

  Params(f0, beta, FU, etagridded, etah, mu, nu, nun)
end

# Element-wise topographic generating function and constant forcing
function Params(g::Grid, f0::Float64, beta::Float64, FU::Float64, eta::Function,
  mu::Float64, nu::Float64, nun::Int)

  FUfunction(t::Float64) = FU

  Params(g, f0, beta, FUfunction, eta, mu, nu, nun)
end

# E Q U A T I O N  T Y P E ---------------------------------------------------- 
type Equation <: AbstractEquation
  LC::Array{Complex{Float64}, 2}  # Element-wise coeff of the eqn's linear part
  calcNL!::Function               # Function to calculate eqn's nonlinear part
end

function Equation(p::Params, g::Grid)
  # Function calcNL! is defined below.
  LC = -p.mu - p.nu * g.KKrsq.^(0.5*p.nun)
  LC[1, 1] = -p.mu
  Equation(LC, calcNL!)
end

# V A R S  T Y P E ------------------------------------------------------------ 
type Vars <: AbstractVars

  t::Float64
  sol::Array{Complex{Float64}, 2}

  # Auxiliary vars
  zet::Array{Float64, 2}
  U::Float64
  u::Array{Float64, 2}
  v::Array{Float64, 2}
  uUq::Array{Float64, 2}
  vq::Array{Float64, 2}
  psi::Array{Float64, 2}

  # Solution
  zeth::Array{Complex{Float64}, 2}
  uh::Array{Complex{Float64}, 2}
  vh::Array{Complex{Float64}, 2}
  uUqh::Array{Complex{Float64}, 2}
  vqh::Array{Complex{Float64}, 2}
  psih::Array{Complex{Float64}, 2}

end

function Vars(g::Grid)
  # Initialize with t=0
  t = 0.0
  sol  = zeros(Complex{Float64}, g.nkr, g.nl)

  # Vorticity auxiliary vars
  zet  = zeros(Float64, g.nx, g.ny)
  U    = 0.0
  u    = zeros(Float64, g.nx, g.ny)
  v    = zeros(Float64, g.nx, g.ny)
  uUq  = zeros(Float64, g.nx, g.ny)
  vq   = zeros(Float64, g.nx, g.ny)
  psi  = zeros(Float64, g.nx, g.ny)

  zeth = zeros(Complex{Float64}, g.nkr, g.nl)
  uh   = zeros(Complex{Float64}, g.nkr, g.nl)
  vh   = zeros(Complex{Float64}, g.nkr, g.nl)
  uUqh = zeros(Complex{Float64}, g.nkr, g.nl)
  vqh  = zeros(Complex{Float64}, g.nkr, g.nl)
  psih = zeros(Complex{Float64}, g.nkr, g.nl)

  # Random initial condition
  sol = exp.( 2.0*pi*im*rand(g.nkr, g.nl) )

  return Vars(t, sol, zet, U, u, v, uUq, vq, psi, zeth, uh, vh, uUqh, vqh, psih)
end


# D E F A U L T  P R O B L E M S ---------------------------------------------- 
function build_problem(nx::Int, Lx::Float64, nu::Float64, nun::Int)
  g  = Grid(nx, Lx)
  p  = Params(nu, nun, g)
  v  = Vars(g)
  eq = Equation(p, g)
  return eq, v, p, g
end

function sample_problem(nx::Int)

  # Default parameters
  Lx   = 2.0*pi
  f0   = 1.0
  beta = 0.0
  nu   = 1e-4
  nun  = 4
  mu   = 1.0
  FU   = 1.0

  g  = Grid(nx, Lx)
  p  = Params(f0, beta, FU, nu, nun)
  v  = Vars(g)
  eq = Equation(p, g)

  return eq, v, p, g
end

function southern_gaussian_mountain(nx::Int, betastar::Float64)

  # This function constructs a Barotropic QG simulation of flow over a
  # "Gaussian mountain" with parameters that approximately correspond
  # to the environment of the Southern Ocean.

  rho0, H, tau = 1035.0, 4e3, 0.2   # density, mean depth, wind stress

  Lx   = 2.0*pi*800e3               # Domain size (meters)
  mu   = 6.3e-6                     # s^(-1)
  FU   = tau/(rho0*H)               # Forcing

  nu   = 1e-4
  nun  = 4

  f0   = -1.26e-4                   # Central inertial frequency
  beta = 1.14e-11                   # Inertial frequency gradient

  # Construct topography generator
  Leta = Lx/20.0
  h    = f0*betastar / (beta*H*Leta)
  eta(x, y) = f0*h/H * exp( -(x^2+y^2)/(2*Leta) )
   
  g  = Grid(nx, Lx)
  p  = Params(g, f0, beta, FU, eta, mu, nu, nun)
  v  = Vars(g)
  eq = Equation(p, g)

  return g, p, v, eq
end





# -----------------------------------------------------------------------------
# Solver ----------------------------------------------------------------------
# -----------------------------------------------------------------------------
function calcNL!(NL::Array{Complex{Float64}, 2}, sol::Array{Complex{Float64}, 2},
  t::Float64, v::Vars, p::Params, g::Grid)

  # Note: U is stored in sol[1, 1]; the other elements of sol are zeth.

  A_mul_B!( v.zet, g.irfftplan, sol )

  v.uh .=    im .* g.Lr .* g.invKKrsq .* sol
  v.vh .= (-im) .* g.Kr .* g.invKKrsq .* sol

  A_mul_B!(v.u, g.irfftplan, v.uh )
  A_mul_B!(v.v, g.irfftplan, v.vh )

  v.uUq .= (sol[1, 1].re .+ v.u).*(v.zet .+ p.eta)
  v.vq  .= v.v.*(v.zet .+ p.eta)

  A_mul_B!(v.uUqh, g.rfftplan, v.uUq )
  A_mul_B!(v.vqh,  g.rfftplan, v.vq  )

  # Nonlinear term for zeta
  NL .= (-im) .* g.Kr.*v.uUqh .- im .* g.Lr.*v.vqh - p.beta.*v.vh

  # 'Nonlinear' term for U with topo correlation.
  # Note: < v*eta > = sum( vh*eta* )
  NL[1, 1] = p.FU(t) - sum(v.vh.*p.etah)

end

# -----------------------------------------------------------------------------
# Helper functions ------------------------------------------------------------
# -----------------------------------------------------------------------------
function updatevars!(v::Vars, g::Grid)

  v.U     = v.sol[1, 1]
  v.zeth .= v.sol
  v.zeth[1, 1] = 0.0 + 0.0*im

  A_mul_B!(v.zet, g.irfftplan, v.zeth )

  v.uh .=    im .* g.Lr .* g.invKKrsq .* sol
  v.vh .= (-im) .* g.Kr .* g.invKKrsq .* sol

  A_mul_B!( v.u, g.irfftplan, v.uh )
  A_mul_B!( v.v, g.irfftplan, v.vh )

  v.psih .= .- v.zeth .* g.invKKrsq

  A_mul_B!(v.psi, g.irfftplan, v.psih)

end

function set_zeta!(v::Vars, g::Grid, zeta::Array{Float64, 2})
  # Set vorticity
  A_mul_B!( v.sol, g.rfftplan, zeta )
  updatevars!(v, g)
end


# Include solver-related functions from solvers.jl.
#include("../solvers.jl")

end


# E N D    M O D U L E    T W O D T U R B --------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
