module TwoDTurb
using FourierFlows, Requires

"""
    InitialValueProblem(; parameters...)

Construct an 2D turbulence problem.
"""
function Problem(; nx=256, Lx=2π, ny=nx, Ly=Lx, nu=0.0, nnu=1, mu=0.0, nmu=0, dt=0.01, stepper="RK4", calcF=nothing)
  if calcF != nothing # the problem is forced
     g = TwoDGrid(nx, Lx, ny, Ly)
    pr = TwoDTurb.ForcedParams(nu, nnu, mu, nmu, calcF)
    vs = TwoDTurb.ForcedVars(g)
    eq = TwoDTurb.Equation(pr, g)
    ts = FourierFlows.autoconstructtimestepper(stepper, dt, eq.LC, g)
  else # initial value problem
     g = TwoDGrid(nx, Lx, ny, Ly)
    pr = TwoDTurb.Params(nu, nnu, mu, nmu)
    vs = TwoDTurb.Vars(g)
    eq = TwoDTurb.Equation(pr, g)
    ts = FourierFlows.autoconstructtimestepper(stepper, dt, eq.LC, g)
  end
  FourierFlows.Problem(g, vs, pr, eq, ts)
end

InitialValueProblem(; kwargs...) = Problem(; kwargs...)
ForcedProblem(; kwargs...) = Problem(; kwargs...)

"""
    Params(nu, nnu, mu, nmu)

Returns the params for unforced two-dimensional turbulence.
"""
struct Params{T} <: AbstractParams
  nu::T      # Vorticity viscosity
  nnu::Int   # Vorticity hyperviscous order
  mu::T      # Bottom drag or hypoviscosity
  nmu::Int   # Order of hypodrag
end
Params(nu, nnu) = Params(nu, nnu, 0, 0)

"""
    ForcedParams(nu, nnu, mu, nmu, calcF!)

Returns the params for forced two-dimensional turbulence.
"""
struct ForcedParams{T} <: AbstractParams
  nu::T              # Vorticity viscosity
  nnu::Int           # Vorticity hyperviscous order
  mu::T              # Bottom drag or hypoviscosity
  nmu::Int           # Order of hypodrag
  calcF!::Function   # Function that calculates the forcing F
end

"""
    Equation(p, g)

Returns the equation for two-dimensional turbulence with params p and grid g.
"""
function Equation(p::Params, g)
  LC = -p.nu*g.KKrsq.^p.nnu .- p.mu*g.KKrsq.^p.nmu
  LC[1, 1] = 0
  FourierFlows.Equation{typeof(p.nu),2}(LC, calcN_advection!)
end

function Equation(p::ForcedParams, g)
  LC = -p.nu*g.KKrsq.^p.nnu - p.mu*g.KKrsq.^p.nmu
  LC[1, 1] = 0
  FourierFlows.Equation{typeof(p.nu),2}(LC, calcN_forced!)
end

# Construct Vars types
       physicalvars = [:q, :U, :V]
      transformvars = [:qh, :Uh, :Vh]
forcedtransformvars = [:qh, :Uh, :Vh, :Fh, :prevsol]

eval(FourierFlows.structvarsexpr(:Vars, physicalvars, transformvars))
eval(FourierFlows.structvarsexpr(:ForcedVars, physicalvars, forcedtransformvars))

"""
    Vars(g)

Returns the vars for unforced two-dimensional turbulence with grid g.
"""
function Vars(g)
  @createarrays typeof(g.Lx) (g.nx, g.ny) q U V
  @createarrays Complex{typeof(g.Lx)} (g.nkr, g.nl) sol qh Uh Vh
  Vars(q, U, V, qh, Uh, Vh)
end

"""
    ForcedVars(g)

Returns the vars for unforced two-dimensional turbulence with grid g.
"""
function ForcedVars(g)
  @createarrays typeof(g.Lx) (g.nx, g.ny) q U V
  @createarrays Complex{typeof(g.Lx)} (g.nkr, g.nl) qh Uh Vh Fh prevsol
  ForcedVars(q, U, V, qh, Uh, Vh, Fh, prevsol)
end


# ------------------
# CUDA functionality
# ------------------

@require CuArrays begin

function CuInitialValueProblem(; stepper="RK4", kwargs...)
  prob = InitialValueProblem(; kwargs...)
  dt = prob.ts.dt

   g = CuTwoDGrid(prob.grid)
  vs = CuVars(prob.vars)
  eq = CuEquation(prob.eqn)
  ts = FourierFlows.autoconstructtimestepper(stepper, dt, eq.LC, g)

  FourierFlows.Problem(g, vs, prob.params, eq, ts)
end

eval(FourierFlows.structvarsexpr(:CuVars, physicalvars, transformvars, arraytype=:CuArray))
eval(FourierFlows.structvarsexpr(:CuForcedVars, physicalvars, transformvars, arraytype=:CuArray))

CuVars(v::Vars) = CuVars(getfield.(v, fieldnames(v))...)
CuVars(g::AbstractGrid) = CuVars(Vars(g))

CuForcedVars(v::Vars) = CuForcedVars(getfield.(v, fieldnames(v))...)
CuForcedVars(g::AbstractGrid) = CuForcedVars(ForcedVars(g))

end # CUDA stuff


# -------
# Solvers
# -------

function calcN_advection!(N, sol, t, s, v, p, g)
  @. v.Uh =  im * g.l  * g.invKKrsq * sol
  @. v.Vh = -im * g.kr * g.invKKrsq * sol
  @. v.qh = sol

  A_mul_B!(v.U, g.irfftplan, v.Uh)
  A_mul_B!(v.V, g.irfftplan, v.Vh)
  A_mul_B!(v.q, g.irfftplan, v.qh)

  @. v.U *= v.q # U*q
  @. v.V *= v.q # V*q

  A_mul_B!(v.Uh, g.rfftplan, v.U) # \hat{U*q}
  A_mul_B!(v.Vh, g.rfftplan, v.V) # \hat{U*q}

  @. N = -im*g.kr*v.Uh - im*g.l*v.Vh
  nothing
end

function calcN_forced!(N, sol, t, s, v, p, g)
  calcN_advection!(N, sol, t, s, v, p, g)
  if t == s.t # not a substep
    v.prevsol .= s.sol # used to compute budgets when forcing is stochastic
    p.calcF!(v.Fh, sol, t, s, v, p, g)
  end
  @. N += v.Fh
  nothing
end


# ----------------
# Helper functions
# ----------------

"""
    updatevars!(v, s, g)

Update the vars in v on the grid g with the solution in s.sol.
"""
function updatevars!(v, s, g)
  v.qh .= s.sol
  @. v.Uh =  im * g.l  * g.invKKrsq * s.sol
  @. v.Vh = -im * g.kr * g.invKKrsq * s.sol
  A_mul_B!(v.q, g.irfftplan, deepcopy(v.qh))
  A_mul_B!(v.U, g.irfftplan, deepcopy(v.Uh))
  A_mul_B!(v.V, g.irfftplan, deepcopy(v.Vh))
  nothing
end

updatevars!(prob::AbstractProblem) = updatevars!(prob.vars, prob.state, prob.grid)

"""
    set_q!(prob, q)
    set_q!(s, v, g, q)

Set the solution s.sol as the transform of q and update variables v
on the grid g.
"""
function set_q!(s, v, g, q)
  A_mul_B!(s.sol, g.rfftplan, q)
  s.sol[1, 1] = 0 # zero domain average
  updatevars!(v, s, g)
  nothing
end
set_q!(prob::AbstractProblem, q) = set_q!(prob.state, prob.vars, prob.grid, q)

"""
    energy(prob)
    energy(s, v, g)

Returns the domain-averaged kinetic energy in the Fourier-transformed vorticity
solution s.sol.
"""
@inline function energy(s, v, g)
  @. v.Uh = g.invKKrsq * abs2(s.sol)
  1/(2*g.Lx*g.Ly)*FourierFlows.parsevalsum(v.Uh, g)
end

@inline energy(prob) = energy(prob.state, prob.vars, prob.grid)

"""
    enstrophy(prob)
    enstrophy(s, g)

Returns the domain-averaged enstrophy in the Fourier-transformed vorticity
solution s.sol.
"""
@inline function enstrophy(s, g)
  1/(2*g.Lx*g.Ly)*FourierFlows.parsevalsum2(s.sol, g)
end

@inline enstrophy(prob) = enstrophy(prob.state, prob.grid)

"""
    dissipation(prob)
    dissipation(s, v, p, g)

Returns the domain-averaged dissipation rate. nnu must be >= 1.
"""
@inline function dissipation(s, v, p, g)
  @. v.Uh = g.KKrsq^(p.nnu-1) * abs2(s.sol)
  @. v.Uh[1, 1] = 0
  p.nu/(g.Lx*g.Ly)*FourierFlows.parsevalsum(v.Uh, g)
end

@inline dissipation(prob::AbstractProblem) = dissipation(prob.state, prob.vars, prob.params, prob.grid)
  
"""
    work(prob)
    work(s, v, p, g)

Returns the domain-averaged rate of work of energy by the forcing Fh.
"""
@inline function work(s, v::ForcedVars, g)
  @. v.Uh = g.invKKrsq * (v.prevsol + s.sol)/2.0 * conj(v.Fh) # Stratonovich
  # @. v.Uh = g.invKKrsq * v.prevsol * conj(v.Fh)             # Ito
  1/(g.Lx*g.Ly)*FourierFlows.parsevalsum(v.Uh, g)
end

@inline work(prob::AbstractProblem) = work(prob.state, prob.vars, prob.grid)

"""
    drag(prob)
    drag(s, v, p, g)

Returns the extraction of domain-averaged energy by drag mu.
"""
@inline function drag(s, v, p, g)
  @. v.Uh = g.KKrsq^(p.nmu-1) * abs2(s.sol)
  @. v.Uh[1, 1] = 0
  p.mu/(g.Lx*g.Ly)*FourierFlows.parsevalsum(v.Uh, g)
end

@inline drag(prob::AbstractProblem) = drag(prob.state, prob.vars, prob.params, prob.grid)

end # module
