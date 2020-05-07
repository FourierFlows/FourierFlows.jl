struct Equation{T, TL, G<:AbstractFloat}
       L :: TL
  calcN! :: Function
    grid :: AbstractGrid{G}
    dims :: Tuple
       T :: T # eltype or tuple of eltypes of sol and N
end

function Equation(L, calcN!, grid::AbstractGrid{G}; dims=supersize(L), T=nothing) where G
  T = T == nothing ? T = cxtype(G) : T
  Equation(L, calcN!, grid, dims, T)
end

mutable struct Clock{T<:AbstractFloat}
    dt :: T
     t :: T
  step :: Int
end

"""
    Problem(sol, clock, grid, eqn, vars, params, timestepper)

Initialize a FourierFlows problem on grid g, with variables v, parameters p,
equation eq, and timestepper ts.
"""
struct Problem{T, A<:AbstractArray, Tg<:AbstractFloat, TL}
          sol :: A
        clock :: Clock{Tg}
          eqn :: Equation{T, TL, Tg}
         grid :: AbstractGrid{Tg}
         vars :: AbstractVars
       params :: AbstractParams
  timestepper :: AbstractTimeStepper{A}
end

struct EmptyParams <: AbstractParams end
struct EmptyVars <: AbstractVars end

"""
    Problem(eqn::Equation, stepper, dt, grid::AbstractGrid,
            vars=EmptyVars, params=EmptyParams, dev::Device=CPU(); stepperkwargs...)

Construct a `Problem` for `eqn` using the time`stepper` with timestep `dt`, on `grid` and 
`dev`ice and with optional `vars`, and `params`. The `stepperkwargs` are passed to the time-stepper
constructor.
"""
function Problem(eqn::Equation, stepper, dt, grid::AbstractGrid{T, A}, 
                 vars=EmptyVars, params=EmptyParams, dev::Device=CPU(); stepperkwargs...) where {T, A}
                 
  clock = Clock{T}(dt, 0, 0)

  timestepper = TimeStepper(stepper, eqn, dt, dev; stepperkwargs...)

  sol = devzeros(dev, eqn.T, eqn.dims)

  return Problem(sol, clock, eqn, grid, vars, params, timestepper)
end

