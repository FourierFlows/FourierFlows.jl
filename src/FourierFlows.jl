__precompile__()

module FourierFlows
import Base: getfield

export AbstractGrid,
       AbstractParams,
       AbstractVars,
       AbstractEquation,
       AbstractTimeStepper,
       AbstractProblem
       
export Equation, Problem, State, DualState

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
  t::Float64
  step::Int
end

#=
For v1.0 release:
  The `t` and `step` properties can be removed from Problem, instead
  overloading the `getfield` function to intercept attempts to access
  `t` and `step`, which will be redirected to prob.state.t and 
  prob.state.step (and perhaps sol as well). This'll make lots of things
  just a little bit nicer.

function getfield(prob::AbstractProblem, name)
  if name ∈ [:t, :step]
    return getfield(prob.state, name)
  else
    return getfield(prob, name)
  end
end
=#


# Time-steppers lists
steppers = [
  "ForwardEuler",
  "FilteredForwardEuler",
  "AB3",
  "RK4",
  "ETDRK4",
  "FilteredETDRK4",
]

filteredsteppers = [
  "FilteredForwardEuler",
  "FilteredETDRK4",
]

"""
Returns a time-stepper type defined by the prefix 'stepper', timestep dt
solution sol (used to construct variables with identical type and size as
the solution vector), and grid g.
"""
function autoconstructtimestepper(stepper, dt, sol, g::AbstractGrid=ZeroDGrid(1))
  fullsteppername = Symbol(stepper, :TimeStepper)
  if stepper ∈ filteredsteppers
    tsexpr = Expr(:call, fullsteppername, dt, sol, g)
  else
    tsexpr = Expr(:call, fullsteppername, dt, sol)
  end

  eval(tsexpr)
end

function autoconstructtimestepper(stepper, dt, solc, solr)
  fullsteppername = Symbol(stepper, :TimeStepper)
  tsexpr = Expr(:call, fullsteppername, dt, solc, solr)
  eval(tsexpr)
end

# Include base functionality
include("domains.jl")
include("diagnostics.jl")
include("output.jl")
include("timesteppers.jl")
include("utils.jl")


function Problem(g::AbstractGrid, v::AbstractVars, p::AbstractParams,
  eq::AbstractEquation, ts::AbstractTimeStepper, st::AbstractState)
  Problem(g, v, p, eq, ts, st, st.t, st.step)
end

function Problem(g::TwoDGrid, v, p, eq, ts)
  st = State{T,ndims(eq.LC)}(0.0, 0, zeros(eq.LC))
  Problem(g, v, p, eq, ts, st)
end


# Include physics modules
include("physics/twodturb.jl")
include("physics/barotropicqg.jl")
include("physics/twomodeboussinesq.jl")
include("physics/traceradvdiff.jl")
include("physics/tracerpatcheqn.jl")

end
