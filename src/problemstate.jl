export Equation, Problem, State, DualState

mutable struct State{T,dim} <: AbstractState
  t::Float64
  step::Int
  dt::Float64
  sol::Array{T,dim}
end

State(T::DataType, sz::Tuple, dt) = State(0.0, 0, dt, zeros(T, sz))

mutable struct DualState{T,dimc,dimr} <: AbstractState
  t::Float64
  step::Int
  dt::Float64
  solc::Array{T,dimc}
  solr::Array{T,dimr}
end

DualState(T::DataType, sizec, sizer, dt) = DualState(0, 0, dt, zeros(T, sizec), zeros(T, sizer))

"""
This type defines the linear implicit and explicit components of an equation.
The linear implicit part of an equation is defined by an array of coefficients
which multiply the solution. The explicit part of an equation is calculated
by a function that may define linear and nonlinear parts.
"""
struct Equation{T,dim} <: AbstractEquation
  LC::Array{T,dim} # Coeffs of the eqn's implicit linear part
  calcN!::Function # Function that calcs linear & nonlinear parts
end

struct DualEquation{T,dimc,dimr} <: AbstractEquation
  LCc::Array{T,dimc}
  LCr::Array{T,dimr}
  calcN!::Function
end

"""
    Problem(g, v, p, eq, ts)

Initialize a FourierFlows problem on grid g, with variables v, parameters p,
equation eq, and timestepper ts.
"""
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

Problem(g, v, p, eq, ts, st) = Problem(g, v, p, eq, ts, st, st.t, st.step)

function Problem(g, v, p, eq, ts; cxsol=true) 
  Tsol = cxsol ? cxeltype(eq.LC) : eltype(eq.LC)
  Problem(g, v, p, eq, ts, State(Tsol, size(eq.LC), ts.dt))
end

#=
For Julia v1.0 release:
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
