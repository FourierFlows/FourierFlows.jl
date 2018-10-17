struct Equation{TL:<AbstractArray,T:<AbstractFloat}
  L::TL 
  calcN!::Function
  grid::AbstractGrid{T}
end

mutable struct Clock{T<:AbstractFloat}
  t::T
  dt::T
  step::Int
end

"""
    Problem(sol, clock, grid, eqn, vars, params, timestepper)

Initialize a FourierFlows problem on grid g, with variables v, parameters p,
equation eq, and timestepper ts.
"""
struct Problem{T:<AbstractFloat,TL:<AbstractArray,Ts:<AbstractArray}
  sol::Ts
  clock::Clock{T}
  grid::AbstractGrid{T}
  eqn::Equation{TL,T}
  vars::AbstractVars
  params::AbstractParams
  timestepper::AbstractTimeStepper{Ts}
end
