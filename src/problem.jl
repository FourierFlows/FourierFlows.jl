struct Equation{T,TL,Tg<:AbstractFloat}
  L::TL
  calcN!::Function
  grid::AbstractGrid{Tg}
  dims::Tuple
  T::T # eltype or tuple of eltypes of sol and N
end

function Equation(L, calcN!, grid::AbstractGrid{Tg}; dims=supersize(L), T=nothing) where Tg
  T = T == nothing ? T = cxtype(Tg) : T
  Equation(L, calcN!, grid, dims, T)
end

mutable struct Clock{T<:AbstractFloat}
  dt::T
  t::T
  step::Int
end

"""
    Problem(sol, clock, grid, eqn, vars, params, timestepper)

Initialize a FourierFlows problem on grid g, with variables v, parameters p,
equation eq, and timestepper ts.
"""
struct Problem{T,Ta<:AbstractArray,Tg<:AbstractFloat,TL}
  sol::Ta
  clock::Clock{Tg}
  eqn::Equation{T,TL,Tg}
  grid::AbstractGrid{Tg}
  vars::AbstractVars
  params::AbstractParams
  timestepper::AbstractTimeStepper{Ta}
end

struct EmptyParams <: AbstractParams end
struct EmptyVars <: AbstractVars end

function Problem(eqn::Equation, stepper, dt, grid::AbstractGrid{T}, vars=EmptyVars, params=EmptyParams) where T
  clock = Clock{T}(dt, 0, 0)
  timestepper = TimeStepper(stepper, eqn, dt)
  sol = superzeros(eqn.T, eqn.dims)
  Problem(sol, clock, eqn, grid, vars, params, timestepper)
end
