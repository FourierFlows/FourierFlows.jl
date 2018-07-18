module TracerAdvDiff

using FourierFlows
export Params, Vars, Equation, set_c!, updatevars!

abstract type AbstractTracerParams <: AbstractParams end

# Problems
function ConstDiffProblem(;
  grid = nothing,
    nx = 128,
    Lx = 2π,
    ny = nx,
    Ly = Lx,
   kap = 1.0,
   eta = kap,
     u = nothing,
     v = nothing,
    dt = 0.01,
  stepper = "RK4",
  steadyflow = false,
  )

  # Defaults
  if u != nothing;   uin = u 
  elseif steadyflow; uin(x, y) = 0.0
  else;              uin(x, y, t) = 0.0
  end

  if v != nothing;   vin = v 
  elseif steadyflow; vin(x, y) = 0.0
  else;              vin(x, y, t) = 0.0
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


# Params
struct ConstDiffParams{T} <: AbstractTracerParams
  eta::T                   # Constant isotropic horizontal diffusivity
  kap::T                   # Constant isotropic vertical diffusivity
  kaph::T                  # Constant isotropic hyperdiffusivity
  nkaph::T                 # Constant isotropic hyperdiffusivity order
  u::Function            # Advecting x-velocity
  v::Function            # Advecting y-velocity
end
ConstDiffParams(eta, kap, u, v) = ConstDiffParams(eta, kap, 0eta, 0eta, u, v)

struct ConstDiffSteadyFlowParams{T} <: AbstractTracerParams
  eta::T                   # Constant horizontal diffusivity
  kap::T                   # Constant vertical diffusivity
  kaph::T                  # Constant isotropic hyperdiffusivity
  nkaph::T                 # Constant isotropic hyperdiffusivity order
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
ConstDiffSteadyFlowParams(eta, kap, u, v, g) = ConstDiffSteadyFlowParams(eta, kap, 0eta, 0eta, u, v, g)


"""
Initialize an equation with constant diffusivity problem parameters p
and on a grid g.
"""
function Equation(p::ConstDiffParams, g)
  LC = zeros(g.Kr)
  @. LC = -p.eta*g.kr^2 - p.kap*g.l^2
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

struct Vars{T} <: AbstractVars
  c::Array{T,2}
  cx::Array{T,2}
  cy::Array{T,2}
  ch::Array{Complex{T},2}
  cxh::Array{Complex{T},2}
  cyh::Array{Complex{T},2}
end

""" Initialize the vars type on a grid g with zero'd arrays and t=0. """
function Vars(g::TwoDGrid)
  @createarrays typeof(g.Lx) (g.nx, g.ny) c cx cy
  @createarrays Complex{typeof(g.Lx)} (g.nkr, g.nl) ch cxh cyh
  Vars(c, cx, cy, ch, cxh, cyh)
end


# --
# Solvers
# --

"""
Calculate the advective terms for a tracer equation with constant
diffusivity.
"""
function calcN!(N, sol, t, s, v, p::ConstDiffParams, g)
  @. v.cxh = im*g.kr*sol
  @. v.cyh = im*g.l*sol

  A_mul_B!(v.cx, g.irfftplan, v.cxh) # destroys v.cxh when using fftw
  A_mul_B!(v.cy, g.irfftplan, v.cyh) # destroys v.cyh when using fftw

  @. v.cx = -p.u(g.X, g.Y, s.t)*v.cx - p.v(g.X, g.Y, s.t)*v.cy # copies over v.cx so v.cx = N in physical space
  A_mul_B!(N, g.rfftplan, v.cx)
  nothing
end


"""
Calculate the advective terms for a tracer equation with constant
diffusivity and time-constant flow.
"""
function calcN_steadyflow!(N, sol, t, s, v, p::ConstDiffSteadyFlowParams, g)
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

""" Update state variables. """
function updatevars!(s, v, g)
  v.ch .= s.sol
  ch1 = deepcopy(v.ch)
  A_mul_B!(v.c, g.irfftplan, ch1)
  nothing
end

updatevars!(prob) = updatevars!(prob.state, prob.vars, prob.grid)


""" Set the concentration field of the model with an array. """
function set_c!(s, v, g, c)
  A_mul_B!(s.sol, g.rfftplan, c)
  updatevars!(s, v, g)
  nothing
end


""" Set the concentration field of the model with a function. """
function set_c!(s, v, g, c::Function)
  cgrid = c.(g.X, g.Y)
  A_mul_B!(s.sol, g.rfftplan, cgrid)
  updatevars!(s, v, g)
  nothing
end

set_c!(prob, c) = set_c!(prob.state, prob.vars, prob.grid, c)


end # module
