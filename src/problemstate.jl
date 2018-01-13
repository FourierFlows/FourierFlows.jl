export Equation, Problem, State, DualState

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

import Base: getfield

function getfield(prob::AbstractProblem, name)
  if name ∈ [:t, :step]
    return getfield(prob.state, name)
  else
    return getfield(prob, name)
  end
end
=#

function Problem(g, v, p, eq, ts, st)
  Problem(g, v, p, eq, ts, st, st.t, st.step)
end

function Problem(g, v, p, eq, ts)
  st = State{eltype(eq.LC),ndims(eq.LC)}(0.0, 0, ts.dt, zeros(eq.LC))
  Problem(g, v, p, eq, ts, st)
end


# Problem state 
mutable struct State{T,dim} <: AbstractState
  t::Float64
  step::Int
  dt::Float64
  sol::Array{T,dim}
end

function State(T::DataType, sz::Tuple, dt)
  sol = zeros(T, sz)
  State(0.0, 0, dt, sol)
end

mutable struct DualState{T,dimc,dimr} <: AbstractState
  t::Float64
  step::Int
  dt::Float64
  solc::Array{T,dimc}
  solr::Array{T,dimr}
end

function DualState(T::DataType, sizec, sizer, dt)
  solc = zeros(T, sizec)
  solr = zeros(T, sizer)
  DualState(0.0, 0, dt, solc, solr)
end


# Equation Composite Type
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
