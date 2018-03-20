module KuramotoSivashinsky
using FourierFlows

#=
Solves the KuramotoSivashinsky equation

∂t u + ν∂ₓ⁴u + λ∂ₓ²u + 1/2 * (∂ₓu)² = 0 ,

with ν=λ=1.

=#

export InitialValueProblem, updatevars!, set_u!

"""
    InitialValueProblem(; parameters...)

Construct an initial-value Kuramoto-Sivashinky problem.
"""
function InitialValueProblem(;
       nx = 256,
       Lx = 2π,
       dt = 0.01,
  stepper = "RK4"
  )

  g  = OneDGrid(nx, Lx)
  pr = Params()
  vs = Vars(g)
  eq = Equation(pr, g)
  ts = FourierFlows.autoconstructtimestepper(stepper, dt, eq.LC, g)

  FourierFlows.Problem(g, vs, pr, eq, ts)
end

struct Params <: AbstractParams; end

"""
    Equation(p, g)

Returns the equation for two-dimensional turbulence with params p and grid g.
"""
function Equation(p, g)
  LC = @. g.kr^2 - g.kr^4
  FourierFlows.Equation{1}(LC, calcN!)
end

# Construct Vars type for unforced two-dimensional turbulence
physvars = [:u, :ux, :uxsq]
transvars = [:uxh, :uxsqh]
expr = FourierFlows.structvarsexpr(:Vars, physvars, transvars; vardims=1)
eval(expr)

"""
    Vars(g)

Returns the vars for unforced two-dimensional turbulence with grid g.
"""
function Vars(g)
  @createarrays Float64 (g.nx,) u ux uxsq
  @createarrays Complex{Float64} (g.nkr,) uxh uxsqh
  Vars(u, ux, uxsq, uxh, uxsqh)
end

# Solvers
function calcN!(N, sol, t, s, v, p, g)
  @. v.uxh = im*g.kr*sol
  A_mul_B!(v.ux, g.irfftplan, v.uxh)
  @. v.uxsq = v.ux^2
  A_mul_B!(v.uxsqh, g.rfftplan, v.uxsq)
  @. N = -0.5*v.uxsqh
  dealias!(N, g)
  nothing
end

# Helper functions
"""
    updatevars!(v, s, g)

Update the vars in v on the grid g with the solution in s.sol.
"""
function updatevars!(v, s, g)
  v.uxsqh .= s.sol
  A_mul_B!(v.u, g.irfftplan, v.uxsqh)
  nothing
end

updatevars!(prob::AbstractProblem) = updatevars!(prob.vars, prob.state, prob.grid)


"""
    set_u!(prob, u)
    set_u!(s, v, g, u)

Set the solution prob.state.sol as the transform of u and update variables.
"""
function set_u!(s, v, g, u)
  A_mul_B!(s.sol, g.rfftplan, u)
  updatevars!(v, s, g)
  nothing
end

set_u!(prob::AbstractProblem, u) = set_u!(prob.state, prob.vars, prob.grid, u)

end # module
