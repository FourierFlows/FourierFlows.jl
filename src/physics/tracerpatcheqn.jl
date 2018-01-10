module TracerPatchEqn

using FourierFlows

abstract type AbstractTracerParams <: AbstractParams end


# Problem initializer ---------------------------------------------------------
function AnalyticFlowProblem(;
   κ = 1.0,
   η = nothing,
   u = nothing,
   v = nothing,
  ux = nothing,
  uy = nothing,
  vx = nothing,
  dt = 0.01,
  timestepper = "RK4"
  )

  if η == nothing; η = κ; end

  # Set flow and flow gradients to zero if functions are not provided.
   if u == nothing;  u★(x, y, t) = 0.0; else  u★ = u;  end
   if v == nothing;  v★(x, y, t) = 0.0; else  v★ = v;  end
  if ux == nothing; ux★(x, y, t) = 0.0; else ux★ = ux; end
  if uy == nothing; uy★(x, y, t) = 0.0; else uy★ = uy; end
  if vx == nothing; vx★(x, y, t) = 0.0; else vx★ = vx; end

   g = FourierFlows.ZeroDGrid(5)
   p = AnalyticFlowParams(κ, η, u★, v★, ux★, uy★, vx★)
   v = Vars()
  eq = Equation()

  if timestepper == "RK4"
    ts = RK4TimeStepper(dt, eq.LC)
  elseif timestepper == "ForwardEuler"
    ts = ForwardEulerTimeStepper(dt, eq.LC)
  end

  FourierFlows.Problem(g, v, p, eq, ts)
end


# Params ----------------------------------------------------------------------
type AnalyticFlowParams <: AbstractTracerParams
  κ::Function
  η::Function
  u::Function
  v::Function
  ux::Function
  uy::Function
  vx::Function
end

function AnalyticFlowParams(κ::Float64, η::Float64, u, v, ux, uy, vx)
  κf(x) = κ
  ηf(z) = η
  AnalyticFlowParams(κf, ηf, u, v, ux, uy, vx)
end

function AnalyticFlowParams(κ::Float64, u, v, ux, uy, vx)
  AnalyticFlowParams(κ, κ, u, v, ux, uy, vx)
end

function AnalyticFlowParams(κ::Float64, η, u, v, ux, uy, vx)
  κf(x) = κ
  AnalyticFlowParams(κf, η, u, v, ux, uy, vx)
end

function AnalyticFlowParams(κ, η::Float64, u, v, ux, uy, vx)
  ηf(z) = η
  AnalyticFlowParams(κ, ηf, u, v, ux, uy, vx)
end




# Equation --------------------------------------------------------------------
type Equation <: AbstractEquation
  LC::Array{Float64, 1}
  calcNL!::Function
end

function Equation()
  Equation(zeros(5), calcNL!)
end


# Vars ------------------------------------------------------------------------
type Vars <: AbstractVars
  t::Float64
  sol::Array{Float64, 1} # x, y, u, v, ux, uy, vx.
end

function Vars()
  Vars(0.0, zeros(5))
end


# Solver ----------------------------------------------------------------------
function calcNL!(NL::Array{Float64, 1}, sol::Array{Float64, 1}, t::Float64,
  v::Vars, p::AnalyticFlowParams, g::AbstractGrid)

  # Key:
  # ξ   = sol[1]
  # ζ   = sol[2]
  # m₁₁ = sol[3]
  # m₁₂ = sol[4]
  # m₂₂ = sol[5]

  # Recall
  # \dot α = ux(ξ, ζ, t)
  # \dot β = uy(ξ, ζ, t)
  # \dot γ = vx(ξ, ζ, t)
  # κ = κ(ξ)
  # η = η(ζ)

  # ∂t ξ = u(ξ, ζ)
  # ∂t ζ = v(ξ, ζ)
  NL[1] = p.u(sol[1], sol[2], t)
  NL[2] = p.v(sol[1], sol[2], t)

  # ∂t m₁₁ = 2η + 2m₁₁dα + 2m₁₂dβ
  NL[3] = (
    2.0*p.η(sol[2])
    + 2.0*sol[3]*p.ux(sol[1], sol[2], t)
    + 2.0*sol[4]*p.uy(sol[1], sol[2], t)
  )

  # ∂t m₁₂ = m₂₂dβ + m₁₁dγ
  NL[4] = sol[5]*p.uy(sol[1], sol[2], t) + sol[3]*p.vx(sol[1], sol[2], t)

  # ∂t m₂₂ = 2κ - 2m₂₂dα + 2m₁₂dγ
  NL[5] = (
    2.0*p.κ(sol[1])
    - 2.0*sol[5]*p.ux(sol[1], sol[2], t)
    + 2.0*sol[4]*p.vx(sol[1], sol[2], t)
  )

  nothing
end


# Helper functions ------------------------------------------------------------
function set_position!(prob::AbstractProblem, ξ₀, ζ₀)
  prob.vars.sol[1:2] = [ξ₀, ζ₀]
  nothing
end

function set_moments!(prob::AbstractProblem, m₁₁, m₁₂, m₂₂)
  prob.vars.sol[3:5] = [m₁₁, m₁₂, m₂₂]
  nothing
end

function set_sol!(prob::AbstractProblem, sol)
  prob.vars.sol = sol
  nothing
end


# end module
end
