module OneDHeatEquation

export
  Problem,
  updatevars!,
  set_c!

using 
  FourierFlows, 
  FFTW

using FourierFlows: varsexpression
using LinearAlgebra: mul!, ldiv!

"""
    Problem(; parameters...)

Construct a constant diffusivity problem with steady or time-varying flow.
"""
function OneDProblem(;
            nx = 128,
            Lx = 2π,
           kap = 0,
            dt = 0.01,
       stepper = "RK4",
             T = Float64
  )

   g = OneDGrid(nx, Lx; T=T)
   p = Params(kap)
   v = Vars(g)
  eq = Equation(p, g)

  FourierFlows.Problem(eq, stepper, dt, grid, vars, params)
end

struct Params{T} <: AbstractParams
  kap::T
end

"""
    Equation(p, g)

Returns the equation for constant diffusivity problem with params p and grid g.
"""
function Equation(p::Params{T}, g) where T<:Number
  FourierFlows.Equation(-g.kr.^2, calcN_implicit!, g)
end

Equation(p::Params{T}, g) where T<:AbstractArray = FourierFlows.Equation(0, calcN!, g)

# Construct Vars types
const physicalvars = [:c, :cx]
const  fouriervars = [:ch, :cxh]

eval(varsexpression(:Vars, physicalvars, fouriervars))

"""
    Vars(g)

Returns the vars for constant diffusivity problem on grid g.
"""
Vars(g::AbstractGrid{T}) where T = Vars(zeros(T, g.nx), zeros(Complex{T}, g.nkr))

"""
    calcN!(N, sol, t, clock, vars, params, grid)

Calculate the nonlinear term for the 1D heat equation.
"""
calcN!(N, sol, t, cl, v, p::Params{T}, g) where T<:Number = nothing
function calcN!(N, sol, t, cl, v, p::Params{T}, g) where T<:AbstractArray
  @. v.cxh = im * g.kr * sol
  ldiv!(v.cx, g.rfftplan, v.cxh)
  @. v.cx *= p.kap
  ldiv!(v.cxh, g.rfftplan, v.cx)
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
  mul!(prob.sol, g.rfftplan, c)
  updatevars(prob)
end

set_c!(prob, c::Function) = set_c!(prob, c.(prob.grid.x))

end # module
