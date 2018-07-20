module TracerAdvDiff
using FourierFlows, Requires

abstract type AbstractTracerParams <: AbstractParams end
abstract type AbstractConstDiffParams <: AbstractParams end
abstract type AbstractSteadyFlowParams <: AbstractParams end

# --
# Problems
# --

"""
    Problem(; parameters...)

Construct a constant diffusivity problem with steady or time-varying flow.
"""
Problem(; kwargs...) = ConstDiffProblem(; kwargs...) # only problem defined for now

function ConstDiffProblem(;
  grid = nothing,
    nx = 128,
    Lx = 2π,
    ny = nx,
    Ly = Lx,
   kap = 0.1,
   eta = kap,
     u = nothing,
     v = nothing,
    dt = 0.01,
  stepper = "RK4",
  steadyflow = false
  )

  zerofunction(args...) = 0.0

  # Defaults
  if u != nothing;   uin = u
  else; uin = zerofunction
  end

  if v != nothing;   vin = v
  else; vin = zerofunction
  end

  if grid == nothing; g = TwoDGrid(nx, Lx, ny, Ly)
  else;               g = grid
  end

  vs = TracerAdvDiff.Vars(g)
  if steadyflow; pr = TracerAdvDiff.ConstDiffSteadyFlowParams(eta, kap, uin, vin, g)
  else;          pr = TracerAdvDiff.ConstDiffParams(eta, kap, uin, vin)
  end
  eq = TracerAdvDiff.Equation(pr, g)
  ts = FourierFlows.autoconstructtimestepper(stepper, dt, eq.LC, g)
  FourierFlows.Problem(g, vs, pr, eq, ts)
end


# --
# Params
# --

"""
    ConstDiffParams(eta, kap, kaph, nkaph, u, v)
    ConstDiffParams(eta, kap, u, v)

Returns the params for constant diffusivity problem with time-varying flow.
"""
struct ConstDiffParams{T} <: AbstractConstDiffParams
  eta::T                 # Constant isotropic horizontal diffusivity
  kap::T                 # Constant isotropic vertical diffusivity
  kaph::T                # Constant isotropic hyperdiffusivity
  nkaph::Int             # Constant isotropic hyperdiffusivity order
  u::Function            # Advecting x-velocity
  v::Function            # Advecting y-velocity
end
ConstDiffParams(eta, kap, u, v) = ConstDiffParams(eta, kap, 0eta, 0, u, v)

"""
    ConstDiffSteadyFlowParams(eta, kap, kaph, nkaph, u, v, g)
    ConstDiffSteadyFlowParams(eta, kap, u, v, g)

Returns the params for constant diffusivity problem with time-steady flow.
"""
struct ConstDiffSteadyFlowParams{T} <: AbstractSteadyFlowParams
  eta::T                 # Constant horizontal diffusivity
  kap::T                 # Constant vertical diffusivity
  kaph::T                # Constant isotropic hyperdiffusivity
  nkaph::Int             # Constant isotropic hyperdiffusivity order
  u::Array{T,2}          # Advecting x-velocity
  v::Array{T,2}          # Advecting y-velocity
end

function ConstDiffSteadyFlowParams(eta, kap, kaph, nkaph, u::Function, v::Function, g)
  ugrid = u.(g.X, g.Y)
  vgrid = v.(g.X, g.Y)
  ConstDiffSteadyFlowParams(eta, kap, kaph, nkaph, ugrid, vgrid)
end

ConstDiffSteadyFlowParams(eta, kap, kaph, nkaph, u, v,
                            g) = ConstDiffSteadyFlowParams{typeof(eta)}(eta, kap, kaph, nkaph, u, v)
ConstDiffSteadyFlowParams(eta, kap, u, v, g) = ConstDiffSteadyFlowParams(eta, kap, 0eta, 0, u, v, g)


# --
# Equations
# --

"""
    Equation(p, g)

Returns the equation for constant diffusivity problem with params p and grid g.
"""
function Equation(p::ConstDiffParams, g)
  LC = zeros(g.Kr)
  @. LC = -p.eta*g.kr^2 - p.kap*g.l^2 - p.kaph*g.KKrsq^p.nkaph
  FourierFlows.Equation{typeof(LC[1, 1]),2}(LC, calcN!)
end

function Equation(p::ConstDiffSteadyFlowParams, g)
  LC = zeros(g.Kr)
  @. LC = -p.eta*g.kr^2 - p.kap*g.l^2 - p.kaph*g.KKrsq^p.nkaph
  FourierFlows.Equation{typeof(LC[1, 1]),2}(LC, calcN_steadyflow!)
end


# --
# Vars
# --

# Construct Vars types
 physicalvars = [:c, :cx, :cy]
transformvars = [:ch, :cxh, :cyh]

eval(FourierFlows.structvarsexpr(:Vars, physicalvars, transformvars))

"""
    Vars(g)

Returns the vars for constant diffusivity problem on grid g.
"""
function Vars(g)
  @createarrays typeof(g.Lx) (g.nx, g.ny) c cx cy
  @createarrays Complex{typeof(g.Lx)} (g.nkr, g.nl) ch cxh cyh
  Vars(c, cx, cy, ch, cxh, cyh)
end



# --
# Solvers
# --

"""
    calcN!(N, sol, t, s, v, p, g)

Calculate the advective terms for a tracer equation with constant diffusivity and time-varying flow.
"""
function calcN!(N, sol, t, s, v, p::AbstractConstDiffParams, g)
  @. v.cxh = im*g.kr*sol
  @. v.cyh = im*g.l*sol

  A_mul_B!(v.cx, g.irfftplan, v.cxh) # destroys v.cxh when using fftw
  A_mul_B!(v.cy, g.irfftplan, v.cyh) # destroys v.cyh when using fftw

  @. v.cx = -p.u(g.X, g.Y, s.t)*v.cx - p.v(g.X, g.Y, s.t)*v.cy # copies over v.cx so v.cx = N in physical space
  A_mul_B!(N, g.rfftplan, v.cx)
  nothing
end


"""
    calcN_steadyflow!(N, sol, t, s, v, p, g)

Calculate the advective terms for a tracer equation with constant diffusivity and time-constant flow.
"""
function calcN_steadyflow!(N, sol, t, s, v, p::AbstractSteadyFlowParams, g)
  @. v.cxh = im*g.kr*sol
  @. v.cyh = im*g.l*sol

  A_mul_B!(v.cx, g.irfftplan, v.cxh) # destroys v.cxh when using fftw
  A_mul_B!(v.cy, g.irfftplan, v.cyh) # destroys v.cyh when using fftw

  @. v.cx = -p.u*v.cx - p.v*v.cy # copies over v.cx so v.cx = N in physical space
  A_mul_B!(N, g.rfftplan, v.cx)
  nothing
end


# --
# Helper functions
# --

"""
    updatevars!(v, s, g)

Update the vars in v on the grid g with the solution in s.sol.
"""
function updatevars!(s, v, g)
  v.ch .= s.sol
  ch1 = deepcopy(v.ch)
  A_mul_B!(v.c, g.irfftplan, ch1)
  nothing
end

updatevars!(prob) = updatevars!(prob.state, prob.vars, prob.grid)


"""
    set_c!(s, v, g, c)
    set_c!(s, v, g, c::Function)
    set_c!(prob, c)

Set the solution s.sol as the transform of c and update variables v
on the grid g.
"""
function set_c!(s, v, g, c)
  A_mul_B!(s.sol, g.rfftplan, c)
  updatevars!(s, v, g)
  nothing
end

function set_c!(s, v, g, c::Function)
  cgrid = c.(g.X, g.Y)
  A_mul_B!(s.sol, g.rfftplan, cgrid)
  updatevars!(s, v, g)
  nothing
end

set_c!(prob, c) = set_c!(prob.state, prob.vars, prob.grid, c)


end # module
