mutable struct State{T,dim} <: AbstractState
  t::Float64
  step::Int
  dt::Float64
  sol::Array{T,dim}
end

State(T::DataType, sz::Tuple, dt) = State(0.0, 0, dt, zeros(T, sz))

mutable struct DualState{Tc,Tr,dimc,dimr} <: AbstractState
  t::Float64
  step::Int
  dt::Float64
  solc::Array{Tc,dimc}
  solr::Array{Tr,dimr}
end

DualState(Tc::DataType, Tr::DataType, sizec, sizer, dt) = DualState{Tc,Tr,length(sizec),length(sizer)}(
  0.0, 0, dt, zeros(Tc, sizec), zeros(Tr, sizer))

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

function Problem(g, v, p, eq::Equation, ts; cxsol=true)
  Tsol = cxsol ? cxeltype(eq.LC) : eltype(eq.LC)
  Problem(g, v, p, eq, ts, State(Tsol, size(eq.LC), ts.dt))
end

function Problem(g, v, p, eq::DualEquation, ts; cxsolc=true, cxsolr=true)
  Tc = cxsolc ? cxeltype(eq.LCc) : eltype(eq.LCc)
  Tr = cxsolr ? cxeltype(eq.LCr) : eltype(eq.LCr)
  Problem(g, v, p, eq, ts, DualState(Tc, Tr, size(eq.LCc), size(eq.LCr), ts.c.dt))
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
