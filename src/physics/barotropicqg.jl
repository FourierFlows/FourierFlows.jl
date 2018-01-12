# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# B A R O T R O P I C Q G >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

module BarotropicQG
using FourierFlows

# P A R A M S -----------------------------------------------------------------
struct Params <: AbstractParams
  f0::Float64                       # Constant planetary vorticity
  beta::Float64                     # Planetary vorticity y-gradient
  FU::Function                      # Time-dependent forcing of domain average flow
  etah::Array{Complex{Float64}, 2}  # Transform of topographic PV
  mu::Float64                       # Linear drag
  nu::Float64                       # Vorticity viscosity
  nun::Int                          # Vorticity hyperviscous order
end

function Params(g::TwoDGrid, f0::Float64, beta::Float64, FU::Real,
  etah::Array{Complex{Float64}, 2}, mu::Float64, nu::Float64, nun::Int)
  # Construct params for constant forcing
  FUfunction(t::Float64) = FU
  Params(f0, beta, FUfunction, etah, mu, nu, nun)
end

function Params(g::TwoDGrid, f0::Float64, beta::Float64, FU::Function,
  eta::Function, mu::Float64, nu::Float64, nun::Int)

  # Construct params given topographic generating function
  etagridded = eta.(g.X, g.Y)
  etah = rfft(etagridded)
  Params(f0, beta, FU, etah, mu, nu, nun)
end

# Element-wise topographic generating function and constant forcing
function Params(g::TwoDGrid, f0::Float64, beta::Float64, FU::Float64,
  eta::Function, mu::Float64, nu::Float64, nun::Int)

  # Construct params given topographic generating function and const forcing
  FUfunction(t::Float64) = FU

  Params(g, f0, beta, FUfunction, eta, mu, nu, nun)
end


# "ConstMeanParams for flows with fixed values of U
struct ConstMeanParams <: AbstractParams
  f0::Float64       # Constant planetary vorticity
  beta::Float64     # Planetary vorticity y-gradient
  U::Float64        # Time-dependent forcing of domain average flow
  etah::Array{Complex{Float64},2}  # Transform of topographic PV
  mu::Float64       # Linear drag
  nu::Float64       # Vorticity viscosity
  nun::Int          # Vorticity hyperviscous order
end

struct FreeDecayParams <: AbstractParams
  f0::Float64       # Constant planetary vorticity
  beta::Float64     # Planetary vorticity y-gradient
  etah::Array{Complex{Float64},2}  # Transform of topographic PV
  mu::Float64       # Linear drag
  nu::Float64       # Vorticity viscosity
  nun::Int          # Vorticity hyperviscous order
end

function FreeDecayParams(g::TwoDGrid, beta::Real, etah::Real, mu::Real,
  nu::Real, nun::Int)
  etaharray = zeros(Complex{Float64}, g.nkr, g.nl)
  FreeDecayParams(1.0, beta, etaharray, mu, nu, nun)
end





# E Q U A T I O N S -----------------------------------------------------------
function Equation(p::Params, g::TwoDGrid)
  LC = -p.mu - p.nu.*g.KKrsq.^(0.5*p.nun)
  FourierFlows.Equation{2}(LC, calcN!)
end

function Equation(p::ConstMeanParams, g::TwoDGrid)
  LC = -p.mu - p.nu.*g.KKrsq.^(0.5*p.nun)
  FourierFlows.Equation{2}(LC, calc_const_mean_N!)
end

function Equation(p::FreeDecayParams, g::TwoDGrid)
  LC = -p.mu - p.nu.*g.KKrsq.^(0.5*p.nun)
  FourierFlows.Equation{2}(LC, calc_free_decay_N!)
end










# V A R S ---------------------------------------------------------------------

struct Vars <: AbstractVars
  q::Array{Float64,2}
  U::Float64
  u::Array{Float64,2}
  v::Array{Float64,2}
  uUq::Array{Float64,2}
  vq::Array{Float64,2}
  psi::Array{Float64,2}
  zeta::Array{Float64,2}
  sp::Array{Float64,2}
  qh::Array{Complex{Float64},2}
  uh::Array{Complex{Float64},2}
  vh::Array{Complex{Float64},2}
  uUqh::Array{Complex{Float64},2}
  vqh::Array{Complex{Float64},2}
  psih::Array{Complex{Float64},2}
  zetah::Array{Complex{Float64},2}
end

function Vars(g::TwoDGrid)
  q     = zeros(Float64, g.nx, g.ny)
  U     = 0.0
  u     = zeros(Float64, g.nx, g.ny)
  v     = zeros(Float64, g.nx, g.ny)
  uUq   = zeros(Float64, g.nx, g.ny)
  vq    = zeros(Float64, g.nx, g.ny)
  psi   = zeros(Float64, g.nx, g.ny)
  zeta  = zeros(Float64, g.nx, g.ny)
  sp    = zeros(Float64, g.nx, g.ny)

  qh    = zeros(Complex{Float64}, g.nkr, g.nl)
  uh    = zeros(Complex{Float64}, g.nkr, g.nl)
  vh    = zeros(Complex{Float64}, g.nkr, g.nl)
  uUqh  = zeros(Complex{Float64}, g.nkr, g.nl)
  vqh   = zeros(Complex{Float64}, g.nkr, g.nl)
  psih  = zeros(Complex{Float64}, g.nkr, g.nl)
  zetah = zeros(Complex{Float64}, g.nkr, g.nl)

  return Vars(t, sol, q, U, u, v, uUq, vq, psi, zeta, sp, qh, uh, vh,
    uUqh, vqh, psih, zetah)
end






struct FreeDecayVars <: AbstractVars
  q::Array{Float64,2}
  u::Array{Float64,2}
  v::Array{Float64,2}
  uq::Array{Float64,2}
  vq::Array{Float64,2}
  psi::Array{Float64,2}
  zeta::Array{Float64,2}
  sp::Array{Float64,2}
  qh::Array{Complex{Float64},2}
  uh::Array{Complex{Float64},2}
  vh::Array{Complex{Float64},2}
  uqh::Array{Complex{Float64},2}
  vqh::Array{Complex{Float64},2}
  psih::Array{Complex{Float64},2}
  zetah::Array{Complex{Float64},2}
end



function FreeDecayVars(g::TwoDGrid)
  t     = 0.0
  sol   = zeros(Complex{Float64}, g.nkr, g.nl)

  q     = zeros(Float64, g.nx, g.ny)
  u     = zeros(Float64, g.nx, g.ny)
  v     = zeros(Float64, g.nx, g.ny)
  uq    = zeros(Float64, g.nx, g.ny)
  vq    = zeros(Float64, g.nx, g.ny)
  psi   = zeros(Float64, g.nx, g.ny)
  zeta  = zeros(Float64, g.nx, g.ny)
  sp    = zeros(Float64, g.nx, g.ny)

  qh    = zeros(Complex{Float64}, g.nkr, g.nl)
  uh    = zeros(Complex{Float64}, g.nkr, g.nl)
  vh    = zeros(Complex{Float64}, g.nkr, g.nl)
  uqh   = zeros(Complex{Float64}, g.nkr, g.nl)
  vqh   = zeros(Complex{Float64}, g.nkr, g.nl)
  psih  = zeros(Complex{Float64}, g.nkr, g.nl)
  zetah = zeros(Complex{Float64}, g.nkr, g.nl)

  return FreeDecayVars(t, sol, q, u, v, uq, vq, psi, zeta, sp, qh, uh, vh,
    uqh, vqh, psih, zetah)
end






# S O L V E R S ---------------------------------------------------------------

function calcN!(N::Array{Complex{Float64}, 2}, sol::Array{Complex{Float64}, 2},
  t::Float64, s::State, v::Vars, p::Params, g::TwoDGrid)
  # Calculate the explicit linear and nonlinear of two equations: one 2D 
  # equation governing the evolution of a barotropic QG flow, and a single
  # 0-dimensional equation for the time evolution of the zonal mean.

  # Note: U is stored in sol[1, 1]; the other elements of sol are qh.
  v.U = sol[1, 1].re
  sol[1, 1] = 0.0

  # This copy is necessary because FFTW's irfft destroys its input.
  v.qh .= sol
  A_mul_B!(v.q, g.irfftplan, v.qh)

  v.uh .=    im .* g.Lr .* g.invKKrsq .* (sol .- p.etah)
  v.vh .= (-im) .* g.Kr .* g.invKKrsq .* (sol .- p.etah)

  A_mul_B!(v.u, g.irfftplan, v.uh)
  A_mul_B!(v.v, g.irfftplan, v.vh)

  v.uUq .= (v.U .+ v.u).*v.q
  v.vq  .= v.v.*v.q

  A_mul_B!(v.uUqh, g.rfftplan, v.uUq)
  A_mul_B!(v.vqh,  g.rfftplan, v.vq)

  # Nonlinear term for q
  N .= (-im) .* g.Kr.*v.uUqh .- im .* g.Lr.*v.vqh - p.beta.*v.vh

  # 'Nonlinear' term for U with topo correlation.
  # Note: < v*eta > = sum( conj(vh)*eta* ) / (nx^2*ny^2)
  N[1, 1] = p.FU(t) - sum(conj(v.vh).*p.etah).re / (g.nx^2.0*g.ny^2.0)

end


function calc_const_mean_N!(N::Array{Complex{Float64}, 2},
  sol::Array{Complex{Float64}, 2}, t::Float64, v::Vars,
  p::ConstMeanParams, g::TwoDGrid)
  # Calculate the nonlinear part of a 2D equation
  # governing the evolution of a barotropic QG flow forced by
  # a constant zonal mean velocity.

  # This copy is necessary because FFTW's irfft destroys its input.
  v.qh .= sol
  A_mul_B!(v.q, g.irfftplan, v.qh)

  v.uh .=    im .* g.Lr .* g.invKKrsq .* (sol .- p.etah)
  v.vh .= (-im) .* g.Kr .* g.invKKrsq .* (sol .- p.etah)

  A_mul_B!(v.u, g.irfftplan, v.uh)
  A_mul_B!(v.v, g.irfftplan, v.vh)

  v.uUq .= (p.U .+ v.u).*v.q
  v.vq  .= v.v.*v.q

  A_mul_B!(v.uUqh, g.rfftplan, v.uUq)
  A_mul_B!(v.vqh,  g.rfftplan, v.vq)

  # Nonlinear term for q
  N .= (-im) .* g.Kr.*v.uUqh .- im .* g.Lr.*v.vqh .- p.beta.*v.vh

end



function calc_free_decay_N!(N::Array{Complex{Float64}, 2},
  sol::Array{Complex{Float64}, 2}, t::Float64, v::FreeDecayVars,
  p::FreeDecayParams, g::TwoDGrid)
  # Calculate the explicit linear and nonlinear part of a 2D equation
  # governing the unforced, free decay of a barotropic QG flow.

  # This copy is necessary because FFTW's irfft destroys its input.
  v.qh .= sol
  A_mul_B!(v.q, g.irfftplan, v.qh)

  v.uh .=    im .* g.Lr .* g.invKKrsq .* (sol .- p.etah)
  v.vh .= (-im) .* g.Kr .* g.invKKrsq .* (sol .- p.etah)

  A_mul_B!(v.u, g.irfftplan, v.uh)

  v.qh .= v.vh  # To preserve vh
  A_mul_B!(v.v, g.irfftplan, v.qh)

  v.uq .= v.u.*v.q
  v.vq .= v.v.*v.q

  A_mul_B!(v.uqh, g.rfftplan, v.uq)
  A_mul_B!(v.vqh, g.rfftplan, v.vq)

  # Nonlinear term for q
  N .= (-im) .* g.Kr.*v.uqh .- im .* g.Lr.*v.vqh .- p.beta.*v.vh

end





# H E L P E R  F U N C T I O N S ----------------------------------------------

function updatevars!(v::Vars, p::Params, g::TwoDGrid)
  # Update state variables, deriving model state from v.sol
  v.U   = v.sol[1, 1].re

  v.qh .= v.sol
  v.qh[1, 1] = 0.0
  v.zetah .= v.qh .- p.etah

  v.q    = irfft(v.qh, g.nx)
  v.zeta = irfft(v.zetah, g.nx)

  v.uh .=    im .* g.Lr .* g.invKKrsq .* v.zetah
  v.vh .= (-im) .* g.Kr .* g.invKKrsq .* v.zetah

  v.u = irfft(v.uh, g.nx)
  v.v = irfft(v.vh, g.nx)

  @. v.psih = -v.qh * g.invKKrsq
  v.psi = irfft(v.psih, g.nx)

  v.sp .= sqrt.( (v.u.+v.U).^2.0 .+ v.v.^2.0 )
end



function updatevars!(v::Vars, p::ConstMeanParams, g::TwoDGrid)
  # Update state variables, deriving model state from v.sol
  v.qh .= v.sol
  v.zetah .= v.qh .- p.etah

  v.q    = irfft(v.qh, g.nx)
  v.zeta = irfft(v.zetah, g.nx)

  v.uh .=    im .* g.Lr .* g.invKKrsq .* v.zetah
  v.vh .= (-im) .* g.Kr .* g.invKKrsq .* v.zetah

  v.u = irfft(v.uh, g.nx)
  v.v = irfft(v.vh, g.nx)

  v.psih .= .- v.zetah .* g.invKKrsq
  v.psi = irfft(v.psih, g.nx)

  v.sp .= sqrt.( (v.u.+p.U).^2.0 .+ v.v.^2.0 )
end




function updatevars!(v::FreeDecayVars, p::FreeDecayParams, g::TwoDGrid)
  # Update state variables, deriving model state from v.sol
  v.qh .= v.sol
  v.zetah .= v.qh .- p.etah

  v.q = irfft(v.qh, g.nx)
  v.zeta = irfft(v.zetah, g.nx)

  v.uh .=    im .* g.Lr .* g.invKKrsq .* v.zetah
  v.vh .= (-im) .* g.Kr .* g.invKKrsq .* v.zetah

  # Do not use A_mul_B here because it destroys its input.
  #A_mul_B!(v.u, g.irfftplan, v.uh)
  #A_mul_B!(v.v, g.irfftplan, v.vh)
  v.u = irfft(v.uh, g.nx)
  v.v = irfft(v.vh, g.nx)

  v.psih .= .- v.zetah .* g.invKKrsq

  # Do not use A_mul_B here because it destroys its input.
  #A_mul_B!(v.psi, g.irfftplan, v.psih)
  v.psi = irfft(v.psih, g.nx)

  v.sp .= sqrt.( v.u.^2.0 .+ v.v.^2.0 )
end





function set_zeta!(v::AbstractVars, p::AbstractParams, g::TwoDGrid,
  zeta::Array{Float64, 2})
  # Set relative vorticity and update model state.

  A_mul_B!(v.zetah, g.rfftplan, zeta)
  v.zetah[1, 1] = 0.0

  v.sol .= v.zetah .+ p.etah

  updatevars!(v, p, g)
end













# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# S E T U P S >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

module Setups

using FourierFlows.BarotropicQG

import FourierFlows: peaked_isotropic_spectrum

export idealizedACC,
       twodturb,
       gaussian_topo_turb


function idealizedACC(nx::Int, etaRms::Float64,
  muStar::Float64, betaStar::Float64, FStar::Float64)
  # Construct a barotropic QG problem with non-dimensional parameters that
  # roughly correspond to the Antarctic Circulpolar Current.

  mu     = muStar*etaRms
  beta   = betaStar*etaRms
  FU     = FStar*mu*etaRms
  f0     = 1.0                  # Central planetary vorticity
  Lx     = 32.0*pi              # Domain size (meters)
  nu     = 1e-6                 # Hyperviscosity
  nun    = 4                    # Order of the hyperviscosity

  # Construct topography
  Leta   = 1.0
  eta  = peaked_isotropic_spectrum(nx, 16.0; rms=etaRms)
  etah = rfft(eta)

  g  = Grid(nx, Lx)
  p  = ConstMeanParams(g, f0, beta, 1.0, etah, mu, nu, nun)
  v  = Vars(g)
  eq = Equation(p, g)

  return g, p, v, eq
end




function random_topo_turb(nx::Int; nu=1e-6, nun=4,
  beta=0.0, mu=0.0, etaRms=0.0, etaScale=8.0, etaMax=0.2)
  # Construct a barotropic QG problem that should reproduce results obtained
  # from a bare 2D turbulence simulation when beta=0.

  Lx     = 2.0*pi               # Domain size (meters)
  f0     = 1.0                  # Central planetary vorticity

  if etaRms > 0.0
    eta = f0*peaked_isotropic_spectrum(nx, etaScale; rms=etaRms)
  else
    eta = f0*peaked_isotropic_spectrum(nx, etaScale; maxval=etaMax)
  end

  etah   = rfft(eta)

  g  = Grid(nx, Lx)
  p  = FreeDecayParams(f0, beta, etah, mu, nu, nun)
  v  = FreeDecayVars(g)
  eq = Equation(p, g)

  # Random initial condition
  set_zeta!(v, p, g, rand(nx, nx))

  return g, p, v, eq
end




function gaussian_topo_turb(nx::Int; nu=1e-6, nun=4,
  beta=0.0, mu=0.0, hbump=0.2, dbump=0.1, topotype="mountain")
  # Construct a barotropic QG problem that should reproduce results obtained
  # from a bare 2D turbulence simulation when beta=0.

  Lx     = 2.0*pi               # Domain size (meters)
  f0     = 1.0                  # Central planetary vorticity
  g      = Grid(nx, Lx)

  Lbump  = dbump*Lx
  x0, y0 = mean(g.x), mean(g.y)

  if topotype == "mountain"
    eta = f0*hbump.*exp.(
      -( (g.X-x0).^2.0 .+ (g.Y-y0).^2.0 ) ./ (2.0*Lbump^2.0) )
  elseif topotype == "east-west ridge"
    eta = f0*hbump.*exp.(
      -( (g.Y-y0).^2.0 ) ./ (2.0*Lbump^2.0) )
  elseif topotype == "north-south ridge"
    eta = f0*hbump.*exp.(
      -( (g.X-x0).^2.0 ) ./ (2.0*Lbump^2.0) )
  end

  etah   = rfft(eta)

  p  = FreeDecayParams(f0, beta, etah, mu, nu, nun)
  v  = FreeDecayVars(g)
  eq = Equation(p, g)

  # Random initial condition
  set_zeta!(v, p, g, rand(nx, nx))

  return g, p, v, eq
end




function forced_gaussian_topo_turb(nx::Int; nu=1e-6, nun=4, U=0.0,
  beta=0.0, mu=0.0, hbump=0.2, dbump=0.1, topotype="mountain")
  # Construct a barotropic QG problem that should reproduce results obtained
  # from a bare 2D turbulence simulation when beta=0.

  Lx     = 2.0*pi               # Domain size (meters)
  f0     = 1.0                  # Central planetary vorticity
  g      = Grid(nx, Lx)

  Lbump  = dbump*Lx
  x0, y0 = mean(g.x), mean(g.y)

  if topotype == "mountain"
    eta    = f0*hbump.*exp.(
      -( (g.X-x0).^2.0 .+ (g.Y-y0).^2.0 ) ./ (2.0*Lbump^2.0) )
  elseif topotype == "east-west ridge"
    eta    = f0*hbump.*exp.(
      -( (g.X-x0).^2.0 ) ./ (2.0*Lbump^2.0) )
  elseif topotype == "north-south ridge"
    eta    = f0*hbump.*exp.(
      -( (g.Y-y0).^2.0 ) ./ (2.0*Lbump^2.0) )
  end

  etah   = rfft(eta)

  p  = ConstMeanParams(f0, beta, U, etah, mu, nu, nun)
  v  = Vars(g)
  eq = Equation(p, g)

  # Random initial condition
  set_zeta!(v, p, g, rand(nx, nx))

  return g, p, v, eq
end







end
# E N D   S E T U P S >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



end
# E N D   B A R O T R O P I C Q G >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
