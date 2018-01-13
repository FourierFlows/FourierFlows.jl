module TwoDTurb

using FourierFlows
Grid = TwoDGrid

# Problem
"""
    InitialValueProblem(; twodturbparameters...)

Construct an initial-value 2D turbulence problem.
"""
function InitialValueProblem(;
     nx = 256,
     Lx = 2π,
     ny = nx,
     Ly = Lx,
      ν = 0.0,
     nν = 1,
     dt = 0.01,
stepper = "RK4"
  )

  g  = TwoDGrid(nx, Lx, ny, Ly)
  pr = TwoDTurb.Params(ν, nν)
  vs = TwoDTurb.Vars(g)
  eq = TwoDTurb.Equation(pr, g)
  ts = FourierFlows.autoconstructtimestepper(stepper, dt, eq.LC, g)
  
  FourierFlows.Problem(g, vs, pr, eq, ts)
end

"""
    ForcedProblem(; twodturbparameters...)

Construct a forced 2D turbulence problem.
"""
function ForcedProblem(;
        nx = 256,
        Lx = 2π,
        ny = nx,
        Ly = Lx,
         ν = 0.0,
        nν = 1,
         μ = 0.0,
        dt = 0.01,
   stepper = "RK4",
     calcF = nothing
  )

  if calcF == nothing; _calcF(F, sol, t, s, v, p, g) = nothing
  else;                _calcF = calcF
  end

  g  = TwoDGrid(nx, Lx, ny, Ly)
  pr = TwoDTurb.ForcedParams(ν, nν, μ, _calcF)
  vs = TwoDTurb.ForcedVars(g)
  eq = TwoDTurb.Equation(pr, g)
  ts = FourierFlows.autoconstructtimestepper(stepper, dt, eq.LC, g)

  FourierFlows.Problem(g, vs, pr, eq, ts)
end

# Params
struct Params <: AbstractParams
  ν::Float64        # Vorticity viscosity
  nν::Int           # Vorticity hyperviscous order
end

struct ForcedParams <: AbstractParams
  ν::Float64        # Vorticity viscosity
  nν::Int           # Vorticity hyperviscous order
  μ::Float64        # Bottom drag
  calcF!::Function  # Function that calculates the forcing F
end

# Equations
function Equation(p::Params, g::TwoDGrid)
  LC = -p.ν * g.KKrsq.^p.nν
  FourierFlows.Equation{2}(LC, calcN_advection!)
end

function Equation(p::ForcedParams, g::TwoDGrid)
  LC = -p.ν*g.KKrsq.^p.nν - p.μ
  FourierFlows.Equation{2}(LC, calcN_forced!)
end




# Vars
physvars = [:q, :U, :V, :Uq, :Vq, :psi]
transvars = [:qh, :Uh, :Vh, :Uqh, :Vqh, :psih]

expr = FourierFlows.getexpr_varstype(:Vars, physvars, transvars)
eval(expr)

function Vars(g::TwoDGrid)
  @createarrays Float64 (g.nx, g.ny) q U V Uq Vq psi
  @createarrays Complex{Float64} (g.nkr, g.nl) sol qh Uh Vh Uqh Vqh psih
  Vars(q, U, V, Uq, Vq, psi, qh, Uh, Vh, Uqh, Vqh, psih)
end


forcedtransvars = [:qh, :Uh, :Vh, :Uqh, :Vqh, :psih, :F]
expr = FourierFlows.getexpr_varstype(:ForcedVars, physvars, forcedtransvars)
eval(expr)

function ForcedVars(g::TwoDGrid)
  @createarrays Float64 (g.nx, g.ny) q U V Uq Vq psi
  @createarrays Complex{Float64} (g.nkr, g.nl) sol qh Uh Vh Uqh Vqh psih F
  ForcedVars(q, U, V, Uq, Vq, psi, qh, Uh, Vh, Uqh, Vqh, psih, F)
end




# Solvers
function calcN_advection!(
  N::Array{Complex{Float64},2}, sol::Array{Complex{Float64},2},
  t::Float64, s::State, v::AbstractVars, p::AbstractParams, g::TwoDGrid)

  v.qh .= sol
  A_mul_B!(v.q, g.irfftplan, v.qh) # destroys qh when using fftw

  @. v.Uh =  im * g.l  * g.invKKrsq * sol
  @. v.Vh = -im * g.kr * g.invKKrsq * sol

  A_mul_B!(v.U, g.irfftplan, v.Uh)
  A_mul_B!(v.V, g.irfftplan, v.Vh)

  @. v.Uq = v.U * v.q
  @. v.Vq = v.V * v.q

  A_mul_B!(v.Uqh, g.rfftplan, v.Uq)
  A_mul_B!(v.Vqh, g.rfftplan, v.Vq)

  @. N = -im*g.kr*v.Uqh - im*g.l*v.Vqh
  nothing
end

function calcN_forced!(N::Array{Complex{Float64}, 2}, 
                sol::Array{Complex{Float64}, 2}, t::Float64, 
                s::State, v::ForcedVars, p::ForcedParams, g::TwoDGrid)

  calcN_advection!(N, sol, t, s, v, p, g)
  p.calcF!(v.F, sol, t, s, v, p, g)

  @. N += v.F
  nothing
end


# Helper functions
"""
Update solution variables using qh in v.sol.
"""
function updatevars!(s, v, g)
  v.qh .= s.sol
  @. v.psih = -g.invKKrsq * v.qh
  @. v.Uh = -im*g.l  * v.psih
  @. v.Vh =  im*g.kr * v.psih

  qh1 = deepcopy(v.qh)
  Uh1 = deepcopy(v.Uh)
  Vh1 = deepcopy(v.Vh)

  A_mul_B!(v.q, g.irfftplan, qh1)
  A_mul_B!(v.U, g.irfftplan, Uh1)
  A_mul_B!(v.V, g.irfftplan, Vh1)
  nothing
end

function updatevars!(prob::AbstractProblem)
  updatevars!(prob.state, prob.vars, prob.grid)
end


"""
    set_q!(s, v, g, q)

Set the solution s.sol as the transform of q and update variables v 
on the grid g.
"""
function set_q!(s, v, g, q)
  A_mul_B!(s.sol, g.rfftplan, q)
  updatevars!(s, v, g)
end

"""
    set_q!(prob, q)

Set the solution prob.state.sol as the transform of q and update variables.
"""
function set_q!(prob::AbstractProblem, q)
  set_q!(prob.state, prob.vars, prob.grid, q)
end


"""
    energy(s, v, g)

Returns the domain-averaged kinetic energy in the Fourier-transformed vorticity
solution s.sol.
"""
@inline function energy(s, v, g)
  @. v.Uh =  im * g.l  * g.invKKrsq * s.sol
  @. v.Vh = -im * g.kr * g.invKKrsq * s.sol
  1 / (2*g.Lx*g.Ly) * (
    FourierFlows.parsevalsum2(v.Uh, g)+FourierFlows.parsevalsum2(v.Vh, g))
end

@inline function energy(prob)
  energy(prob.state, prob.vars, prob.grid)
end


"""
    enstrophy(s, g)

Returns the domain-averaged enstrophy in the Fourier-transformed vorticity
solution s.sol.
"""
@inline function enstrophy(s, g)
  1/(2*g.Lx*g.Ly)*FourierFlows.parsevalsum2(s.sol, g)
end

@inline function enstrophy(prob)
  enstrophy(prob.state, prob.grid)
end


"""
    dissipation(s, v, p, g)

Returns the domain-averaged dissipation rate.
"""
@inline function dissipation(s, v, p, g)
  @. v.Uqh = (g.kr^(2*(p.nν-1)) + g.l^(2*(p.nν-1))) * abs2(s.sol)
  p.ν*FourierFlows.parsevalsum(v.Uqh, g)
end

@inline function dissipation(prob::AbstractProblem)
  dissipation(prob.state, prob.vars, prob.params, prob.grid)
end


"""
    injection(s, v, p, g)

Returns the domain-averaged rate of injection of energy by the forcing F.
"""
@inline function injection(s, v::ForcedVars, g)
  @. v.psih = -g.invKKrsq * s.sol
  @. v.Uqh = -v.psih*conj(v.F)
  FourierFlows.parsevalsum(v.Uqh, g)
end

@inline function injection(prob::AbstractProblem)
  injection(prob.state, prob.vars, prob.grid)
end


end # module
