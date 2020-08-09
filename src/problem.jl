"""
    Equation{T, TL, G<:AbstractFloat}
    
A struct that includes the equation to be solved `∂u/∂t = L*u + N(u)`. Array `L` 
includes the coefficients of the linear term `L*u` and `calcN!` is a function 
which computes the nonlinear term `N(u)`. The struct also includes the problem's
`grid` and the float type of the state vector (and consequently of `N(u)`).
"""
struct Equation{T, TL, G<:AbstractFloat}
       L :: TL
  calcN! :: Function
    grid :: AbstractGrid{G}
    dims :: Tuple
       T :: T # eltype or tuple of eltypes of sol and N
end

"""
    Equation(L, calcN!, grid; dims=supersize(L), T=nothing)
    
An equation constructor given the array `L` of the coefficients of the linear
term, the function `calcN!` that computes the nonlinear term and the problem's 
`grid`.
"""
function Equation(L, calcN!, grid::AbstractGrid{G}; dims=supersize(L), T=nothing) where G
  T = T == nothing ? T = cxtype(G) : T
  Equation(L, calcN!, grid, dims, T)
end

"""
    Clock{T<:AbstractFloat}
    
A struct containing the time-step `dt`, the time `t` and the `step`-number of the simulation.
"""
mutable struct Clock{T<:AbstractFloat}
    dt :: T
     t :: T
  step :: Int
end

"""
    Problem{T, A<:AbstractArray, Tg<:AbstractFloat, TL}
    
A struct including everything a FourierFlows problem requires: the state vector
`sol`, the `clock`, the equation `eqn`, the `grid`, all problem variables in 
`vars`, problem parameters in `params`, and the `timestepper`.
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

"""
    EmptyParams <: AbstractParams

A placeholder struct for parameters.
"""
struct EmptyParams <: AbstractParams end

"""
    EmptyVars <: AbstractVars

A placeholder struct for variables.
"""
struct EmptyVars <: AbstractVars end

"""
    Problem(eqn::Equation, stepper, dt, grid::AbstractGrid,
            vars=EmptyVars, params=EmptyParams, dev::Device=CPU(); stepperkwargs...)

Construct a `Problem` for equation `eqn` using the time`stepper` with timestep 
`dt`, on `grid` and `dev`ice. Optionally provide variables in `vars` and
parameters with `params`. The `stepperkwargs` are passed to the time-stepper
constructor.
"""
function Problem(eqn::Equation, stepper, dt, grid::AbstractGrid{T}, 
                 vars=EmptyVars, params=EmptyParams, dev::Device=CPU(); stepperkwargs...) where T
                 
  clock = Clock{T}(dt, 0, 0)

  timestepper = TimeStepper(stepper, eqn, dt, dev; stepperkwargs...)

  sol = devzeros(dev, eqn.T, eqn.dims)

  return Problem(sol, clock, eqn, grid, vars, params, timestepper)
end

show(io::IO, clock::FourierFlows.Clock) =
     print(io, "Clock\n",
               "  ├─── timestep dt: ", clock.dt, '\n',
               "  ├────────── step: ", clock.step, '\n',
               "  └─────────time t: ", clock.t)
