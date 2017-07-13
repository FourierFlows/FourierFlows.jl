__precompile__()

include("../fourierflows.jl")

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# B E G I N    M O D U L E    B A R O T R O P I C Q G --------------------------

module BarotropicQG

using FourierFlowTypes, Domains, TimeSteppers


export Grid, Vars, Params, ConstMeanParams, Equation
export ForwardEulerTimeStepper, ETDRK4TimeStepper
export RK4TimeStepper, AB3TimeStepper

export set_q!, updatevars!, calc_NL!, stepforward!

# This module always uses the TwoDGrid
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

# With constant forcing
function Params(g::Grid, f0::Float64, beta::Float64, FU::Real, eta::Array{Float64, 2},
  etah::Array{Complex{Float64}, 2}, mu::Float64, nu::Float64, nun::Int)
  FUfunction(t::Float64) = FU
  Params(f0, beta, FUfunction, eta, etah, mu, nu, nun)
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




# "ConstMeanParams for flows with fixed values of U --------------------------- 
type ConstMeanParams <: AbstractParams
  f0::Float64                       # Constant inertial frequency
  beta::Float64                     # Inertial frequency y-gradient
  U::Float64                        # Time-dependent forcing of domain average flow
  eta::Array{Float64, 2}            # Topographic PV
  etah::Array{Complex{Float64}, 2}  # Transform of topographic PV
  mu::Float64                       # Linear drag
  nu::Float64                       # Vorticity viscosity
  nun::Int                          # Vorticity hyperviscous order
end






# E Q U A T I O N  T Y P E ---------------------------------------------------- 
type Equation <: AbstractEquation
  LC::Array{Complex{Float64}, 2}  # Element-wise coeff of the eqn's linear part
  calcNL!::Function               # Function to calculate eqn's nonlinear part
end

function Equation(p::Params, g::Grid)
  LC = -p.mu - p.nu.*g.KKrsq.^(0.5*p.nun)
  Equation(LC, calcNL!)
end

# Constructor for the "constant mean flow" problem type
function Equation(p::ConstMeanParams, g::Grid)
  # Function calcNL! is defined below.
  LC = -p.mu - p.nu.*g.KKrsq.^(0.5*p.nun)
  Equation(LC, calc_const_mean_NL!)
end







# V A R S  T Y P E ------------------------------------------------------------ 
type Vars <: AbstractVars

  t::Float64
  sol::Array{Complex{Float64}, 2}

  # Auxiliary vars
  q::Array{Float64, 2}
  U::Float64
  u::Array{Float64, 2}
  v::Array{Float64, 2}
  uUq::Array{Float64, 2}
  vq::Array{Float64, 2}
  psi::Array{Float64, 2}
  sp::Array{Float64, 2}

  # Solution
  qh::Array{Complex{Float64}, 2}
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
  q  = zeros(Float64, g.nx, g.ny)
  U    = 0.0
  u    = zeros(Float64, g.nx, g.ny)
  v    = zeros(Float64, g.nx, g.ny)
  uUq  = zeros(Float64, g.nx, g.ny)
  vq   = zeros(Float64, g.nx, g.ny)
  psi  = zeros(Float64, g.nx, g.ny)
  sp   = zeros(Float64, g.nx, g.ny)

  qh = zeros(Complex{Float64}, g.nkr, g.nl)
  uh   = zeros(Complex{Float64}, g.nkr, g.nl)
  vh   = zeros(Complex{Float64}, g.nkr, g.nl)
  uUqh = zeros(Complex{Float64}, g.nkr, g.nl)
  vqh  = zeros(Complex{Float64}, g.nkr, g.nl)
  psih = zeros(Complex{Float64}, g.nkr, g.nl)

  return Vars(t, sol, q, U, u, v, uUq, vq, psi, sp, qh, uh, vh, uUqh, vqh, psih)
end




# S O L R V E R S -------------------------------------------------------------
function calcNL!(NL::Array{Complex{Float64}, 2}, sol::Array{Complex{Float64}, 2},
  t::Float64, v::Vars, p::Params, g::Grid)

  # Note: U is stored in sol[1, 1]; the other elements of sol are qh.
  v.U = sol[1, 1].re
  sol[1, 1] = 0.0
  
  A_mul_B!( v.q, g.irfftplan, sol )

  v.uh .=    im .* g.Lr .* g.invKKrsq .* (sol .- p.etah)
  v.vh .= (-im) .* g.Kr .* g.invKKrsq .* (sol .- p.etah)

  A_mul_B!(v.u, g.irfftplan, v.uh )
  A_mul_B!(v.v, g.irfftplan, v.vh )

  v.uUq .= (v.U .+ v.u).*v.q
  v.vq  .= v.v.*v.q

  A_mul_B!(v.uUqh, g.rfftplan, v.uUq )
  A_mul_B!(v.vqh,  g.rfftplan, v.vq  )

  # Nonlinear term for q
  NL .= (-im) .* g.Kr.*v.uUqh .- im .* g.Lr.*v.vqh - p.beta.*v.vh

  # 'Nonlinear' term for U with topo correlation.
  # Note: < v*eta > = sum( vh*eta* )
  NL[1, 1] = p.FU(t) - sum(v.vh.*p.etah)

end


# ----------------------------------------------------------------------------- 
function calc_const_mean_NL!(NL::Array{Complex{Float64}, 2}, 
  sol::Array{Complex{Float64}, 2}, t::Float64, v::Vars, 
  p::ConstMeanParams, g::Grid)

  # Note: U is stored in sol[1, 1]; the other elements of sol are qh.
  A_mul_B!( v.q, g.irfftplan, sol )

  v.uh .=    im .* g.Lr .* g.invKKrsq .* sol
  v.vh .= (-im) .* g.Kr .* g.invKKrsq .* sol

  A_mul_B!(v.u, g.irfftplan, v.uh )
  A_mul_B!(v.v, g.irfftplan, v.vh )

  v.uUq .= (p.U .+ v.u).*v.q
  v.vq  .= v.v.*v.q

  A_mul_B!(v.uUqh, g.rfftplan, v.uUq )
  A_mul_B!(v.vqh,  g.rfftplan, v.vq  )

  # Nonlinear term for q
  NL .= (-im) .* g.Kr.*v.uUqh .- im .* g.Lr.*v.vqh .- p.beta.*v.vh

end




# -----------------------------------------------------------------------------
# Helper functions ------------------------------------------------------------
# -----------------------------------------------------------------------------
function updatevars!(v::Vars, p::Params, g::Grid)
  v.U   = v.sol[1, 1].re
  v.qh .= v.sol
  v.qh[1, 1] = 0.0

  A_mul_B!(v.q, g.irfftplan, v.qh )

  v.uh .=    im .* g.Lr .* g.invKKrsq .* v.sol
  v.vh .= (-im) .* g.Kr .* g.invKKrsq .* v.sol

  A_mul_B!( v.u, g.irfftplan, v.uh )
  A_mul_B!( v.v, g.irfftplan, v.vh )

  v.psih .= .- v.qh .* g.invKKrsq

  A_mul_B!(v.psi, g.irfftplan, v.psih)

  v.sp .= sqrt.( (v.u.+v.U).^2.0 .+ v.v.^2.0 )
end



function updatevars!(v::Vars, p::ConstMeanParams, g::Grid)
  v.qh .= v.sol
  A_mul_B!(v.q, g.irfftplan, v.qh )

  v.uh .=    im .* g.Lr .* g.invKKrsq .* v.sol
  v.vh .= (-im) .* g.Kr .* g.invKKrsq .* v.sol

  A_mul_B!( v.u, g.irfftplan, v.uh )
  A_mul_B!( v.v, g.irfftplan, v.vh )

  v.psih .= .- v.qh .* g.invKKrsq

  A_mul_B!(v.psi, g.irfftplan, v.psih)

  v.sp .= sqrt.( (v.u.+v.U).^2.0 .+ v.v.^2.0 )
end



function set_q!(v::Vars, p::AbstractParams, g::Grid, q::Array{Float64, 2})
  # Set vorticity
  A_mul_B!( v.sol, g.rfftplan, q )
  updatevars!(v, p, g)
end





end
# E N D    M O D U L E    T W O D T U R B --------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------












# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# B E G I N    M O D U L E    B A R O T R O P I C Q G P R O B L E M S-----------
module BarotropicQGProblems

using BarotropicQG

export peaked_isotropic_spectrum
export southern_gaussian_mountain, nondimensional_southern_ocean
export twodturb

# Sample problems
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



function nondimensional_southern_ocean(nx::Int, etaRms::Float64,
  muStar::Float64, betaStar::Float64, FStar::Float64)

  # This function constructs a Barotropic QG simulation of flow over a
  # "Gaussian mountain" with parameters that approximately correspond
  # to the environment of the Southern Ocean.

  f0     = 1.0                  # Central inertial frequency

  mu     = muStar*etaRms    
  beta   = betaStar*etaRms
  FU     = FStar*mu*etaRms 

  Lx     = 32.0*pi              # Domain size (meters)
  nu     = 1e-6                 # Hyperviscosity
  nun    = 4                    # Order of the hyperviscosity

  # Construct topography
  Leta   = 1.0
  eta  = peaked_isotropic_spectrum(nx, 16.0; rms=etaRms)
  etah = rfft(eta)
   
  g  = Grid(nx, Lx)
  p  = ConstMeanParams(g, f0, beta, 1.0, eta, etah, mu, nu, nun)
  v  = Vars(g)
  eq = Equation(p, g)

  return g, p, v, eq
end


function twodturb(nx::Int, nu::Float64, nun::Int; beta=0.0, U=0.0, mu=0.0)

  Lx     = 2.0*pi               # Domain size (meters)
  f0     = 1.0                  # Central inertial frequency
  mu     = 0.0
  beta   = 0.0
  U      = 0.0
  eta    = zeros(nx, nx)
  etah   = rfft(eta)
   
  g  = Grid(nx, Lx)
  p  = ConstMeanParams(f0, beta, U, eta, etah, mu, nu, nun)
  v  = Vars(g)
  eq = Equation(p, g)

  set_q!(v, p, g, rand(nx, nx))

  return g, p, v, eq
end


end
# E N D    M O D U L E    B A R O T R O P I C Q G P R O B L E M S --------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
