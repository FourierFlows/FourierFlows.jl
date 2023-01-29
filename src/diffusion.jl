"""
A test-bed module that solves the 1D diffusion equation.

# Exports
$(EXPORTS)
"""
module Diffusion

export
  Problem,
  updatevars!,
  set_c!

using
  DocStringExtensions

using FourierFlows

using LinearAlgebra: mul!, ldiv!

"""
    Problem(dev::Device = CPU();
                     nx = 128,
                     Lx = 2π,
                      κ = 0,
                     dt = 0.01,
                stepper = "RK4",           
       aliased_fraction = 0,
                      T = Float64)

Construct a constant diffusivity problem.

Keyword arguments
=================
  - `dev`: (required) `CPU()` or `GPU()`; computer architecture used to time-step `problem`.
  - `nx`: Number of grid points in ``x``-domain.
  - `Lx`: Extent of the ``x``-domain.
  - `κ`: Diffusivity coefficient.
  - `dt`: Time-step.
  - `stepper`: The extent of the ``y``-domain.
  - `aliased_fraction`: the fraction of high-wavenubers that are zero-ed out by `dealias!()`.
  - `T`: `Float64` or `Float32`; floating point type used for `problem` data.
"""
function Problem(dev::Device=CPU();
                         nx = 128,
                         Lx = 2π,
                          κ = 0,
                         dt = 0.01,
                    stepper = "RK4",           
           aliased_fraction = 0,
                          T = Float64)

      grid = OneDGrid(dev; nx, Lx, aliased_fraction, T)
    params = Params(grid, κ)
      vars = Vars(grid)
  equation = Equation(grid, params)

  return FourierFlows.Problem(equation, stepper, dt, grid, vars, params)
end

"""
    struct Params{T} <: AbstractParams

The parameters for diffusion problem:

$(TYPEDFIELDS)
"""
struct Params{T} <: AbstractParams
    "diffusivity coefficient"
    κ :: T
end

Params(grid, κ::Number) = Params(κ)
Params(grid, κ::AbstractArray) = Params(device_array(grid.device)(κ))

"""
    Equation(grid, params)

Return the equation for a constant diffusivity problem on `grid` with diffusivity
found inside `params`.
"""
function Equation(grid::AbstractGrid{T}, params::Params{D}) where {T, D<:Number}
  L = zeros(grid.device, T, grid.nkr)
  @. L = - params.κ * grid.kr^2
  
  return FourierFlows.Equation(L, calcN!, grid)
end

Equation(grid::AbstractGrid{T}, ::Params{D}) where {T, D<:AbstractArray} =
  FourierFlows.Equation(0, calcN!, grid; dims=(grid.nkr,), T=cxtype(T))

"""
    struct Vars{Aphys, Atrans} <: AbstractVars

The variables for diffusion problem:

$(FIELDS)
"""
struct Vars{Aphys, Atrans} <: AbstractVars
    "tracer concentration ``c``"
    c :: Aphys
    "tracer concentration derivative ``∂ₓc``"
   cx :: Aphys
    "Fourier transform of tracer concentration ``c``"
   ch :: Atrans
    "Fourier transform of tracer concentration derivative ``∂ₓc``"
  cxh :: Atrans
end

"""
    Vars(grid)

Return the variables `vars` for a constant diffusivity problem on `grid`.
"""
function Vars(grid::AbstractGrid{T}) where T
  Dev = typeof(grid.device)

  @devzeros Dev T grid.nx c cx
  @devzeros Dev Complex{T} grid.nkr ch cxh

  return Vars(c, cx, ch, cxh)
end

"""
    calcN!(N, sol, t, clock, vars, params, grid)

Calculate the nonlinear term for the 1D diffusion equation.
"""
function calcN!(N, sol, t, clock, vars, params::Params{D}, grid) where D<:Number
  @. N = 0
  
  return nothing
end

function calcN!(N, sol, t, clock, vars, params::Params{D}, grid) where D<:AbstractArray
  @. vars.cxh = im * grid.kr * sol
  ldiv!(vars.cx, grid.rfftplan, vars.cxh)
  @. vars.cx *= params.κ
  mul!(vars.cxh, grid.rfftplan, vars.cx)
  @. N = im * grid.kr * vars.cxh
  
  return nothing
end

"""
    updatevars!(vars, grid, sol)

Update the variables in `vars` on the `grid` with the solution in `sol`.
"""
function updatevars!(vars, grid, sol)
  @. vars.ch  = sol  
  @. vars.cxh = im * grid.kr * sol

  ldiv!(vars.c, grid.rfftplan, deepcopy(vars.ch))
  ldiv!(vars.cx, grid.rfftplan, deepcopy(vars.cxh))
  
  return nothing
end

updatevars!(prob) = updatevars!(prob.vars, prob.grid, prob.sol)

"""
    set_c!(prob, c)

Set the solution `sol` as the transform of `c` and update `vars`.
"""
function set_c!(prob, c)
  A = typeof(prob.vars.c)

  prob.vars.c .= A(c) # ensure that c is converted to the type prob.vars.c (e.g., convert to CuArray if user gives `c` as Array)
  mul!(prob.sol, prob.grid.rfftplan, prob.vars.c)

  updatevars!(prob)
  
  return nothing
end

end # module
