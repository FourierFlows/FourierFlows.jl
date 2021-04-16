module Diffusion

export
  Problem,
  updatevars!,
  set_c!

using Reexport

@reexport using FourierFlows

using LinearAlgebra: mul!, ldiv!

"""
    Problem(; parameters...)

Construct a constant diffusivity problem.
"""
function Problem(;
            nx = 128,
            Lx = 2π,
             κ = 0,
            dt = 0.01,
       stepper = "RK4",
             T = Float64,
           dev = CPU()
  )

      grid = OneDGrid(dev, nx, Lx; T=T)
    params = Params(dev, κ)
      vars = Vars(dev, grid)
  equation = Equation(dev, κ, grid)

  return FourierFlows.Problem(equation, stepper, dt, grid, vars, params, dev)
end

struct Params{T} <: AbstractParams
  κ :: T
end

Params(dev, κ::Number) = Params(κ)
Params(dev, κ::AbstractArray) = Params(ArrayType(dev)(κ))

"""
    Equation(dev, κ, grid)

Returns the equation for a constant diffusivity problem on `grid` with diffusivity κ.
"""
function Equation(dev::Device, κ::T, grid) where T<:Number
  L = zeros(dev, T, grid.nkr)
  @. L = - κ * grid.kr^2
  
  return FourierFlows.Equation(L, calcN!, grid)
end

Equation(dev::Device, κ::T, grid::AbstractGrid{Tg}) where {T<:AbstractArray, Tg} = 
  FourierFlows.Equation(0, calcN!, grid; dims=(grid.nkr,), T=cxtype(Tg))

struct Vars{Aphys, Atrans} <: AbstractVars
    c :: Aphys
   cx :: Aphys
   ch :: Atrans
  cxh :: Atrans
end

"""
    Vars(dev, grid)

Returns the variables `vars` for a constant diffusivity problem on `grid`.
"""
function Vars(::Dev, grid::AbstractGrid{T}) where {Dev, T}
  @devzeros Dev T grid.nx c cx
  @devzeros Dev Complex{T} grid.nkr ch cxh
  Vars(c, cx, ch, cxh)
end

"""
    calcN!(N, sol, t, clock, vars, params, grid)

Calculate the nonlinear term for the 1D heat equation.
"""
function calcN!(N, sol, t, clock, vars, params::Params{T}, grid) where T<:Number
  @. N = 0
  
  return nothing
end

function calcN!(N, sol, t, clock, vars, params::Params{T}, grid) where T<:AbstractArray
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
  ldiv!(vars.c, grid.rfftplan, deepcopy(sol))
  @. vars.ch = sol
  
  return nothing
end

updatevars!(prob) = updatevars!(prob.vars, prob.grid, prob.sol)

"""
    set_c!(prob, c)

Set the solution as the transform of `c` and update `vars`.
"""
function set_c!(prob, c)
  T = typeof(prob.vars.c)
  prob.vars.c .= T(c) # this makes sure that c is converted to the ArrayType used in prob.vars.c (e.g., convert to CuArray if user gives c as Arrray)
  mul!(prob.sol, prob.grid.rfftplan, prob.vars.c)
  updatevars!(prob)
  
  return nothing
end

end # module
