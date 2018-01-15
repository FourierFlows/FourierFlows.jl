module TracerAdvDiff

using FourierFlows
export Params, Vars, Equation, set_c!, updatevars!

abstract type AbstractTracerParams <: AbstractParams end

# Problems
function ConstDiffSteadyFlowProblem(;
  grid = nothing,
    nx = 128,
    Lx = 2π,
    ny = nothing,
    Ly = nothing,
     κ = 1.0,
     η = nothing,
     u = nothing,
     v = nothing,
    dt = 0.01,
  timestepper = "RK4"
  )

  # Defaults
  if ny == nothing; ny = nx;            end
  if Ly == nothing; Ly = Lx;            end
  if  η == nothing; η = κ;              end

  if u == nothing; uin(x, y) = 0.0
  else;            uin = u
  end

  if v == nothing; vin(x, y) = 0.0
  else;            vin = v
  end

  if grid == nothing;
    grid = TwoDGrid(nx, Lx, ny, Ly)
  end

  vs = Vars(grid)
  pr = ConstDiffSteadyFlowParams(η, κ, uin, vin, grid)
  eq = Equation(pr, grid)

  if     timestepper == "RK4";          ts = ETDRK4TimeStepper(dt, eq.LC)
  elseif timestepper == "ETDRK4";       ts = RK4TimeStepper(dt, eq.LC)
  elseif timestepper == "ForwardEuler"; ts = ForwardEulerTimeStepper(dt, eq.LC)
  end

  FourierFlows.Problem(grid, vs, pr, eq, ts; realsol=true, nvars=1)
end


# Params
struct ConstDiffParams <: AbstractTracerParams
  η::Float64                   # Constant isotropic horizontal diffusivity
  κ::Float64                   # Constant isotropic vertical diffusivity
  κh::Float64                  # Constant isotropic hyperdiffusivity
  nκh::Float64                 # Constant isotropic hyperdiffusivity order
  u::Function                  # Advecting x-velocity
  v::Function                  # Advecting y-velocity
end

function ConstDiffParams(η, κ, u, v)
  ConstDiffParams(η, κ, 0, 0, u, v)
end



struct ConstDiffSteadyFlowParams <: AbstractTracerParams
  η::Float64                   # Constant horizontal diffusivity
  κ::Float64                   # Constant vertical diffusivity
  κh::Float64                  # Constant isotropic hyperdiffusivity
  nκh::Float64                 # Constant isotropic hyperdiffusivity order
  u::Array{Float64,2}          # Advecting x-velocity
  v::Array{Float64,2}          # Advecting y-velocity
end

function ConstDiffSteadyFlowParams(η, κ, κh, nκh,
                                   u::Function, v::Function, g::TwoDGrid)
  ugrid = u.(g.X, g.Y)
  vgrid = v.(g.X, g.Y)
  ConstDiffSteadyFlowParams(η, κ, κh, nκh, ugrid, vgrid)
end

function ConstDiffSteadyFlowParams(η, κ, u, v, g)
  ConstDiffSteadyFlowParams(η, κ, 0, 0, u, v)
end




""" Initialize an equation with constant diffusivity problem parameters p
and on a grid g. """
function Equation(p::ConstDiffParams, g::TwoDGrid)
  LC = zeros(g.Kr)
  @. LC = -p.η*g.kr^2 - p.κ*g.l^2
  FourierFlows.Equation{2}(LC, calcN!)
end

function Equation(p::ConstDiffSteadyFlowParams, g::TwoDGrid)
  LC = zeros(g.Kr)
  @. LC = -p.η*g.kr^2 - p.κ*g.l^2 - p.κh*g.KKrsq^p.nκh
  FourierFlows.Equation{2}(LC, calcN_steadyflow!)
end




# Vars
struct Vars <: AbstractVars
  c::Array{Float64,2}
  cu::Array{Float64,2}
  cv::Array{Float64,2}
  ch::Array{Complex{Float64},2}
  cuh::Array{Complex{Float64},2}
  cvh::Array{Complex{Float64},2}
end

""" Initialize the vars type on a grid g with zero'd arrays and t=0. """
function Vars(g::TwoDGrid)
  @createarrays Float64 (g.nx, g.ny) c cu cv
  @createarrays Complex{Float64} (g.nkr, g.nl) ch cuh cvh
  Vars(c, cu, cv, ch, cuh, cvh)
end




# Solvers ---------------------------------------------------------------------
"""
Calculate the advective terms for a tracer equation with constant
diffusivity.
"""
function calcN!(
  N::Array{Complex{Float64},2}, sol::Array{Complex{Float64},2},
  t::Float64, s::State, v::Vars, p::ConstDiffParams, g::TwoDGrid)

  v.ch .= sol
  A_mul_B!(v.c, g.irfftplan, v.ch) # destroys v.ch when using fftw

  v.cu .= p.u.(g.X, g.Y, v.t)
  v.cv .= p.v.(g.X, g.Y, v.t)

  v.cu .*= v.c
  v.cv .*= v.c

  A_mul_B!(v.cuh, g.rfftplan, v.cu)
  A_mul_B!(v.cvh, g.rfftplan, v.cv)

  @. N = -im*g.kr*v.cuh - im*g.l*v.cvh
  nothing
end


"""
Calculate the advective terms for a tracer equation with constant
diffusivity and time-constant flow.
"""
function calcN_steadyflow!(
  N::Array{Complex{Float64},2}, sol::Array{Complex{Float64},2},
  t::Float64, s::State, v::Vars, p::ConstDiffSteadyFlowParams, g::TwoDGrid)

  v.ch .= sol
  A_mul_B!(v.c, g.irfftplan, v.ch) # destroys v.ch when using fftw

  @. v.cu = p.u * v.c
  @. v.cv = p.v * v.c

  A_mul_B!(v.cuh, g.rfftplan, v.cu)
  A_mul_B!(v.cvh, g.rfftplan, v.cv)

  @. N = -im*g.kr*v.cuh - im*g.l*v.cvh
  nothing
end








# Helper functions ------------------------------------------------------------
""" Update state variables. """
function updatevars!(s::State, v::AbstractVars, g::TwoDGrid)
  v.ch .= s.sol
  ch1 = deepcopy(v.ch)
  A_mul_B!(v.c, g.irfftplan, ch1)
  nothing
end

function updatevars!(prob::AbstractProblem)
  updatevars!(prob.state, prob.vars, prob.grid)
end


""" Set the concentration field of the model with an array. """
function set_c!(s::State, v::AbstractVars, g::TwoDGrid, c::Array{Float64, 2})
  A_mul_B!(s.sol, g.rfftplan, c)
  updatevars!(s, v, g)
  nothing
end


""" Set the concentration field of the model with a function. """
function set_c!(s::State, v::AbstractVars, g::TwoDGrid, c::Function)
  cgrid = c.(g.X, g.Y)
  A_mul_B!(s.sol, g.rfftplan, cgrid)
  updatevars!(s, v, g)
  nothing
end

function set_c!(prob::AbstractProblem, c)
  set_c!(prob.state, prob.vars, prob.grid, c)
end


end # module
