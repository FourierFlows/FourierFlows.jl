__precompile__()

module FourierFlows

export AbstractGrid,
       AbstractParams,
       AbstractVars,
       AbstractEquation,
       AbstractTimeStepper,
       AbstractProblem,
       Problem, State, DualState

# Abstract supertypes
abstract type AbstractGrid end
abstract type AbstractParams end
abstract type AbstractVars end
abstract type AbstractEquation end
abstract type AbstractTimeStepper end
abstract type AbstractProblem end
abstract type AbstractState end


mutable struct State{T,dim} <: AbstractState
  t::Float64
  step::Int
  sol::Array{T,dim}
end

function State(T::DataType, sz::Tuple)
  sol = zeros(T, sz)
  State(0.0, 0, sol)
end

mutable struct DualState{T,dimc,dimr} <: AbstractState
  t::Float64
  step::Int
  solc::Array{T,dimc}
  solr::Array{T,dimr}
end

function DualState(T::DataType, sizec, sizer)
  solc = zeros(T, sizec)
  solr = zeros(T, sizer)
  DualState(0.0, 0, solc, solr)
end


"""
This type defines the linear implicit and explicit components of an equation.
The linear implicit part of an equation is defined by an array of coefficients
which multiply the solution. The explicit part of an equation is calculated
by a function that may define linear and nonlinear parts.
"""
struct Equation{dim} <: AbstractEquation
  LC::Array{Complex{Float64},dim} # Coeffs of the eqn's implicit linear part
  calcN!::Function # Function that calcs linear & nonlinear parts
end

struct DualEquation{dimc,dimr} <: AbstractEquation
  LCc::Array{Complex{Float64},dimc}
  LCr::Array{Complex{Float64},dimr}
  calcN!::Function
end


# Problem type and associated functions
mutable struct Problem <: AbstractProblem
  grid::AbstractGrid
  vars::AbstractVars
  params::AbstractParams
  eqn::AbstractEquation
  ts::AbstractTimeStepper
  state::AbstractState
  t::Real
  step::Int
end


# Include base functionality
include("domains.jl")
include("diagnostics.jl")
include("output.jl")
include("timesteppers.jl")
include("utils.jl")


function Problem(g::AbstractGrid, v::AbstractVars, p::AbstractParams,
  eq::AbstractEquation, ts::AbstractTimeStepper, st::AbstractState)
  Problem(g, v, p, eq, ts, st, 0.0, 0)
end

function Problem(g::TwoDGrid, v, p, eq, ts; nvars=1, realsol=true, 
                 T=Complex{Float64})
  if nvars == 1
    if realsol; sol = zeros(T, (g.nkr, g.nl))
    else;       sol = zeros(T, (g.nk, g.nl))
    end 
  else
    if realsol; sol = zeros(T, (g.nkr, g.nl, nvars))
    else;       sol = zeros(T, (g.nk, g.nl, nvars))
    end 
  end
  st = State{T,ndims(sol)}(0.0, 0, sol)

  Problem(g, v, p, eq, ts, st)
end


# Include physics modules
include("physics/twodturb.jl")
include("physics/barotropicqg.jl")
include("physics/twomodeboussinesq.jl")
include("physics/traceradvdiff.jl")
include("physics/tracerpatcheqn.jl")

end
