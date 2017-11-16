__precompile__()


module TracerAdvDiff

using FourierFlows

export Params,
       Vars,
       Equation

export set_c!, updatevars!




# Problems -------------------------------------------------------------------- 
function ConstDiffProblem(nx, Lx, ny=nx, Ly=Lx, η::Real, κ=η, u::Function, 
  v::Function)
  g  = TwoDGrid(nx, Lx, ny, Ly)
  vs = Vars(g)
  pr = ConstDiffParams(η, κ, u, v)
  eq = Equation(pr, g)
  ts = RK4TimeStepper(dt, eq.LC)

  FourierFlows.Problem(g, vs, pr, eq, ts)
end




function ConstDiffSteadyFlowProblem(nx::Int, Lx, ny=nx, Ly=Lx, η::Real, 
  κ=η, u, v)
  g  = TwoDGrid(nx, Lx, ny, Ly)
  vs = Vars(g)
  pr = ConstDiffSteadyFlowParams(η, κ, u, v)
  eq = Equation(pr, g)
  ts = RK4TimeStepper(dt, eq.LC)

  FourierFlows.Problem(g, vs, pr, eq, ts)
end



# Params ---------------------------------------------------------------------- 
abstract type AbstractTracerParams <: AbstractParams end

type ConstDiffParams <: AbstractTracerParams
  η::Float64                   # Constant isotropic horizontal diffusivity
  κ::Float64                   # Constant isotropic vertical diffusivity
  u::Function                    # Advecting x-velocity
  v::Function                    # Advecting y-velocity
end

function ConstDiffParams(η::Real, κ=η, u::Real, v::Real)
  ufunc(x, y, t) = u
  vfunc(x, y, t) = v
  ConstDiffParams(η, κ, ufunc, vfunc)
end




type ConstDiffSteadyFlowParams
  η::Float64                   # Constant isotropic horizontal diffusivity
  κ::Float64                   # Constant isotropic vertical diffusivity
  u::Array{Float64, 2}         # Advecting x-velocity
  v::Array{Float64, 2}         # Advecting y-velocity
end

function ConstDiffSteadyFlowParams(η, κ=η, u::Function, v::Function, 
  g::TwoDGrid)
  ugrid = u.(g.X, g.Y)
  vgrid = v.(g.X, g.Y)
  ConstDiffSteadyFlowParams(η, κ, ugrid, vgrid)
end




# Equations ------------------------------------------------------------------- 
type Equation <: AbstractEquation
  LC::Array{Complex{Float64}, 2}  # Element-wise coeff of the eqn's linear part
  calcNL!::Function               # Function to calculate eqn's nonlinear part
end

""" Initialize an equation with constant diffusivity problem parameters p
and on a grid g. """
function Equation(p::ConstDiffParams, g::TwoDGrid)
  LC = -p.κ.*g.Kr.^2.0 - p.η.*g.Lr.^2.0
  Equation(LC, calcNL!)
end

function Equation(p::ConstDiffSteadyFlowParams, g::TwoDGrid)
  LC = -p.κ.*g.Kr.^2.0 - p.η.*g.Lr.^2.0
  Equation(LC, calcNL_steadyflow!)
end




# Vars
type Vars <: AbstractVars
  t::Float64
  sol::Array{Complex{Float64}, 2}
  c::Array{Float64, 2}
  cu::Array{Float64, 2}
  cv::Array{Float64, 2}

  ch::Array{Complex{Float64}, 2}
  cuh::Array{Complex{Float64}, 2}
  cvh::Array{Complex{Float64}, 2}
end

""" Initialize the vars type on a grid g with zero'd arrays and t=0. """
function Vars(g::TwoDGrid)
  t     = 0.0
  sol   = zeros(Complex{Float64}, g.nkr, g.nl)

  c     = zeros(Float64, g.nx, g.ny)
  cu    = zeros(Float64, g.nx, g.ny)
  cv    = zeros(Float64, g.nx, g.ny)

  ch    = zeros(Complex{Float64}, g.nkr, g.nl)
  cuh   = zeros(Complex{Float64}, g.nkr, g.nl)
  cvh   = zeros(Complex{Float64}, g.nkr, g.nl)
  return Vars(t, sol, c, cu, cv, ch, cuh, cvh)
end




# Solvers --------------------------------------------------------------------- 
function calcNL!(NL::Array{Complex{Float64}, 2}, 
  sol::Array{Complex{Float64}, 2}, 
  t::Float64, v::Vars, p::ConstDiffParams, g::TwoDGrid)
  
  # Calculate the advective terms for a tracer equation with constant
  # diffusivity.

  # This copy is necessary because FFTW's irfft destroys its input.
  v.ch .= sol
  A_mul_B!(v.c, g.irfftplan, v.ch)

  v.cu .= p.u.(g.X, g.Y, v.t)
  v.cv .= p.v.(g.X, g.Y, v.t)

  v.cu .*= v.c
  v.cv .*= v.c

  A_mul_B!(v.cuh, g.rfftplan, v.cu)
  A_mul_B!(v.cvh, g.rfftplan, v.cv)

  NL .= (-im).*g.Kr.*v.cuh .- im.*g.Lr.*v.cvh

end




function calcNL_steadyflow!(NL::Array{Complex{Float64}, 2}, 
  sol::Array{Complex{Float64}, 2}, 
  t::Float64, v::Vars, p::ConstDiffParams, g::TwoDGrid)
  
  # Calculate the advective terms for a tracer equation with constant
  # diffusivity and time-constant flow.

  # This copy is necessary because FFTW's irfft destroys its input.
  v.ch .= sol
  A_mul_B!(v.c, g.irfftplan, v.ch)

  @. v.cu .= v.u*v.c
  @. v.cv .= v.v*v.c

  A_mul_B!(v.cuh, g.rfftplan, v.cu)
  A_mul_B!(v.cvh, g.rfftplan, v.cv)

  @. NL = -im*g.Kr*v.cuh - im*g.Lr*v.cvh

end








# Helper functions ------------------------------------------------------------ 
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
