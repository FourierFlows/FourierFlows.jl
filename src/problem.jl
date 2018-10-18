struct Equation{T,TL<:AbstractArray,Tg<:AbstractFloat}
  L::TL 
  calcN!::Function
  grid::AbstractGrid{Tg}
  dims::Tuple
  T::T # element type of the solution
end

Equation(L, calcN!, grid; dims=size(L)) = Equation(L, calcN!, grid, dims)

mutable struct Clock{T<:AbstractFloat}
  dt::T
  t::T
  step::Int
end

Clock(dt) = Clock(dt, 0dt, 0)

"""
    Problem(sol, clock, grid, eqn, vars, params, timestepper)

Initialize a FourierFlows problem on grid g, with variables v, parameters p,
equation eq, and timestepper ts.
"""
struct Problem{T<:AbstractFloat,TL<:AbstractArray,Ta<:AbstractArray}
  sol::Ta
  clock::Clock{T}
  eqn::Equation{TL,T}
  grid::AbstractGrid{T}
  vars::AbstractVars
  params::AbstractParams
  timestepper::AbstractTimeStepper{Ta}
end

struct EmptyParams <: AbstractParams end
struct EmptyVars <: AbstractVars end

function Problem(eqn::Equation, stepper, dt, grid, vars=EmptyVars, params=EmptyParams)
  clock = Clock{fltype(eqn.T)}(dt)
  timestepper = TimeStepper(stepper, eqn, dt)
  Problem(clock, eqn, grid, vars, params, timestepper)
end

function Problem(clock, eqn, grid, vars, params, timestepper)
  sol = superzeros(eqn.T, eqn.dims)
  Problem(sol, clock, eqn, grid, vars, params, timestepper)
end
