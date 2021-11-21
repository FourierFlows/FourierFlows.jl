"""
    struct Equation{T, TL, G<:AbstractFloat}
    
The equation to be solved `∂u/∂t = L*u + N(u)`. Array `L` includes the coefficients
of the linear term `L*u` and `calcN!` is a function which computes the nonlinear
term `N(u)`. The struct also includes the problem's `grid` and the float type of the
state vector (and consequently of `N(u)`).

$(TYPEDFIELDS)
"""
struct Equation{T, TL, G<:AbstractFloat}
    "array with coefficient for the linear part of the equation"
       L :: TL
    "function that computes the nonlinear part of the equation"
  calcN! :: Function
    "the grid"
    grid :: AbstractGrid{G}
    "the dimensions of `L`"
    dims :: Tuple
    "the float type for the state vector"
       T :: T # eltype or tuple of eltypes of sol and N
end

"""
    Equation(L, calcN!, grid; dims=supersize(L), T=nothing)
    
The equation constructor from the array `L` of the coefficients of the linear term, the function 
`calcN!` that computes the nonlinear term and the `grid` for the problem.
"""
function Equation(L, calcN!, grid::AbstractGrid{G}; dims=supersize(L), T=nothing) where G
  T = T == nothing ? T = cxtype(G) : T
  
  return Equation(L, calcN!, grid, dims, T)
end

"""
    mutable struct Clock{T<:AbstractFloat}
    
Represents the clock of a problem.

$(TYPEDFIELDS)
"""
mutable struct Clock{T<:AbstractFloat}
    "the time-step"
    dt :: T
    "the time"
     t :: T
    "the step number"
  step :: Int
end

"""
    struct Problem{T, A<:AbstractArray, Tg<:AbstractFloat, TL}
    
A problem that represents a partial differential equation.

$(TYPEDFIELDS)
"""
struct Problem{T, A<:AbstractArray, Tg<:AbstractFloat, TL}
    "the state vector"
          sol :: A
    "the problem's"
        clock :: Clock{Tg}
    "the equation"
          eqn :: Equation{T, TL, Tg}
    "the grid"
         grid :: AbstractGrid{Tg}
    "the variables"
         vars :: AbstractVars
    "the parameters"
       params :: AbstractParams
    "the timestepper"
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
    Problem(eqn::Equation, stepper, dt, grid::AbstractGrid{T}, 
            vars=EmptyVars, params=EmptyParams, dev::Device=CPU(); stepperkwargs...) where T

Construct a `Problem` for equation `eqn` using the time`stepper` with timestep 
`dt`, on `grid` and on `dev`ice. Optionally, use the keyword arguments to provide 
variables with `vars` and parameters with `params`. The `stepperkwargs` are passed
to the time-stepper constructor.
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
               "  └──────── time t: ", clock.t)

show(io::IO, eqn::FourierFlows.Equation) =
     print(io, "Equation\n",
               "  ├──────── linear coefficients: L", '\n',
               "  │                              ├───type: ", eltype(eqn.L), '\n',
               "  │                              └───size: ", size(eqn.L), '\n', 
               "  ├───────────── nonlinear term: calcN!()", '\n',
               "  └─── type of state vector sol: ", eqn.T)

show(io::IO, problem::FourierFlows.Problem) =
    print(io, "Problem\n",
              "  ├─────────── grid: grid (on " * FourierFlows.griddevice(problem.grid) * ")", '\n',
              "  ├───── parameters: params", '\n',
              "  ├────── variables: vars", '\n',
              "  ├─── state vector: sol", '\n',
              "  ├─────── equation: eqn", '\n',
              "  ├────────── clock: clock", '\n',
              "  └──── timestepper: ", string(nameof(typeof(problem.timestepper))))
