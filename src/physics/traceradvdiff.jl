__precompile__()


module TracerAdvDiff

using FourierFlows

export Params,
       Vars,
       Equation

export set_c!, updatevars!




# Params
abstract type AbstractTracerParams <: AbstractParams end

type ConstDiffParams <: AbstractTracerParams
  eta::Float64                   # Constant isotropic horizontal diffusivity
  kap::Float64                   # Constant isotropic vertical diffusivity
  u::Function                    # Advecting horizontal velocity
  w::Function                    # Advecting vertical velocity (here z=y)
end

function ConstDiffParams(eta::Real, kap::Real, u::Real, w::Real)
  ufunc(x, z, t) = u
  wfunc(x, z, t) = w
  ConstDiffParams(eta, kap, ufunc, wfunc)
end

function ConstDiffParams(kap::Real, u::Function, w::Function)
  ConstDiffParams(kap, kap, u, w)
end

function ConstDiffParams(kap::Real, u::Real, w::Real)
  ufunc(x::Float64, z::Float64, t::Float64) = u
  wfunc(x::Float64, z::Float64, t::Float64) = w
  ConstDiffParams(kap, ufunc, wfunc)
end




# Equations
type Equation <: AbstractEquation
  LC::Array{Complex{Float64}, 2}  # Element-wise coeff of the eqn's linear part
  calcNL!::Function               # Function to calculate eqn's nonlinear part
end

""" Initialize an equation with constant diffusivity problem parameters p
and on a grid g. """
function Equation(p::ConstDiffParams, g::TwoDGrid)
  LC = -p.kap.*g.Kr.^2.0 - p.eta.*g.Lr.^2.0
  Equation(LC, calcNL!)
end




# Vars
type Vars <: AbstractVars
  t::Float64
  sol::Array{Complex{Float64}, 2}
  c::Array{Float64, 2}
  cu::Array{Float64, 2}
  cw::Array{Float64, 2}

  ch::Array{Complex{Float64}, 2}
  cuh::Array{Complex{Float64}, 2}
  cwh::Array{Complex{Float64}, 2}
end

""" Initialize the vars type on a grid g with zero'd arrays and t=0. """
function Vars(g::TwoDGrid)
  t     = 0.0
  sol   = zeros(Complex{Float64}, g.nkr, g.nl)

  c     = zeros(Float64, g.nx, g.ny)
  cu    = zeros(Float64, g.nx, g.ny)
  cw    = zeros(Float64, g.nx, g.ny)

  ch    = zeros(Complex{Float64}, g.nkr, g.nl)
  cuh   = zeros(Complex{Float64}, g.nkr, g.nl)
  cwh   = zeros(Complex{Float64}, g.nkr, g.nl)
  return Vars(t, sol, c, cu, cw, ch, cuh, cwh)
end




# Solvers
function calcNL!(NL::Array{Complex{Float64}, 2}, 
  sol::Array{Complex{Float64}, 2}, 
  t::Float64, v::Vars, p::ConstDiffParams, g::TwoDGrid)
  
  # Calculate the advective terms for a tracer equation with constant
  # diffusivity.

  # This copy is necessary because FFTW's irfft destroys its input.
  v.ch .= sol
  A_mul_B!(v.c, g.irfftplan, v.ch)

  v.cu .= p.u.(g.X, g.Y, v.t)
  v.cw .= p.w.(g.X, g.Y, v.t)

  v.cu .*= v.c
  v.cw .*= v.c

  A_mul_B!(v.cuh, g.rfftplan, v.cu)
  A_mul_B!(v.cwh, g.rfftplan, v.cw)

  NL .= (-im).*g.Kr.*v.cuh .- im.*g.Lr.*v.cwh

end




# Helper functions

""" Update state variables. """
function updatevars!(v::AbstractVars, p::AbstractTracerParams, g::TwoDGrid)
  v.ch  .= v.sol
  v.c    = irfft(v.ch, g.nx)
end

function updatevars!(prob::AbstractProblem)
  updatevars!(prob.vars, prob.params, prob.grid)
end


""" Set the concentration field of the model with an array. """
function set_c!(v::AbstractVars, p::AbstractTracerParams, g::TwoDGrid,
  c::Array{Float64, 2})
  A_mul_B!(v.ch, g.rfftplan, c)
  v.sol .= v.ch
  updatevars!(v, p, g)
end


""" Set the concentration field of the model with a function. """
function set_c!(v::AbstractVars, p::AbstractTracerParams, g::TwoDGrid,
  c::Function)
  cgrid = c.(g.X, g.Y)
  A_mul_B!(v.ch, g.rfftplan, cgrid)
  v.sol .= v.ch
  updatevars!(v, p, g)
end

function set_c!(prob::AbstractProblem, c::Function)
  set_c!(prob.vars, prob.params, prob.grid, c)
end

function set_c!(prob::AbstractProblem, c::Array{Float64, 2})
  set_c!(prob.vars, prob.params, prob.grid, c)
end






end
# end module TracerAdvDiff
