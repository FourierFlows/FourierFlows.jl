module Diffusion

export
  Problem,
  updatevars!,
  set_c!

using
  FFTW,
  Reexport

@reexport using FourierFlows

using FourierFlows: varsexpression
using LinearAlgebra: mul!, ldiv!

"""
    Problem(; parameters...)

Construct a constant diffusivity problem with steady or time-varying flow.
"""
function Problem(;
            nx = 128,
            Lx = 2π,
         kappa = 0,
            dt = 0.01,
       stepper = "RK4",
             T = Float64,
           dev = CPU()
  )

    grid = OneDGrid(dev, nx, Lx; T=T)
  params = Params(dev, kappa)
    vars = Vars(dev, grid)
     eqn = DiffusionEquation(dev, kappa, grid)

  FourierFlows.Problem(eqn, stepper, dt, grid, vars, params, dev)
end

struct Params{T} <: AbstractParams
  kappa::T
end

Params(dev, kappa::Number) = Params(kappa)
Params(dev, kappa::AbstractArray) = Params(ArrayType(dev)(kappa))

"""
    DiffusionEquation(dev, kappa, grid)

Returns the equation for constant diffusivity problem with diffusivity kappa and grid.
"""
function DiffusionEquation(dev::Device, kappa::T, grid) where T<:Number
  L = zeros(dev, T, grid.nkr)
  @. L = -kappa * grid.kr^2
  FourierFlows.Equation(L, calcN!, grid)
end

function DiffusionEquation(dev::Device, kappa::T, grid::AbstractGrid{Tg}) where {T<:AbstractArray, Tg}
  FourierFlows.Equation(0, calcN!, grid; dims=(grid.nkr,), T=cxtype(Tg))
end

struct Vars{Aphys, Atrans} <: AbstractVars
    c :: Aphys
   cx :: Aphys
   ch :: Atrans
  cxh :: Atrans
end

"""
    Vars(dev, grid)

Returns the vars for constant diffusivity problem on grid g.
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
function calcN!(N, sol, t, cl, v, p::Params{T}, g) where T<:Number
  @. N = 0
  nothing
end

function calcN!(N, sol, t, cl, v, p::Params{T}, g) where T<:AbstractArray
  @. v.cxh = im * g.kr * sol
  ldiv!(v.cx, g.rfftplan, v.cxh)
  @. v.cx *= p.kappa
  mul!(v.cxh, g.rfftplan, v.cx)
  @. N = im*g.kr*v.cxh
  nothing
end

"""
    updatevars!(v, g, sol)

Update the vars in v on the grid g with the solution in `sol`.
"""
function updatevars!(v, g, sol)
  ldiv!(v.c, g.rfftplan, deepcopy(sol))
  @. v.ch = sol
  nothing
end

updatevars!(prob) = updatevars!(prob.vars, prob.grid, prob.sol)

"""
    set_c!(prob, c)

Set the solution as the transform of `c`.
"""
function set_c!(prob, c)
  mul!(prob.sol, prob.grid.rfftplan, c)
  updatevars!(prob)
end

set_c!(prob, c::Function) = set_c!(prob, c.(prob.grid.x))

end # module
