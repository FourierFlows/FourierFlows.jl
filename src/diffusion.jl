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
             T = Float64
  )

    grid = OneDGrid(nx, Lx; T=T)
  params = Params(kappa)
    vars = Vars(grid)
     eqn = Equation(kappa, grid)

  FourierFlows.Problem(eqn, stepper, dt, grid, vars, params)
end

struct Params{T} <: AbstractParams
  kappa::T
end

"""
    Equation(p, g)

Returns the equation for constant diffusivity problem with params p and grid g.
"""
function Equation(kappa::T, g) where T<:Number
  FourierFlows.Equation(-kappa*g.kr.^2, calcN!, g)
end

function Equation(kappa::T, g::AbstractGrid{Tg}) where {T<:AbstractArray,Tg}
  FourierFlows.Equation(0, calcN!, g; dims=(g.nkr,), T=cxtype(Tg))
end

# Construct Vars types
const physicalvars = [:c, :cx]
const  fouriervars = [:ch, :cxh]

eval(varsexpression(:Vars, physicalvars, fouriervars))

"""
    Vars(g)

Returns the vars for constant diffusivity problem on grid g.
"""
function Vars(g::AbstractGrid{T}) where T
  @zeros T g.nx c cx
  @zeros Complex{T} g.nkr ch cxh
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
