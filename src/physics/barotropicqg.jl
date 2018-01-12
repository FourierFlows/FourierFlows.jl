module BarotropicQG
using FourierFlows
Grid = TwoDGrid

# Params
struct Params <: AbstractParams
  f0::Float64      # Constant planetary vorticity
  beta::Float64    # Planetary vorticity y-gradient
  FU::Function     # Time-dependent forcing of domain average flow
  eta::Array{Float64,2}  # Topographic PV
  mu::Float64      # Linear drag
  nu::Float64      # Vorticity viscosity
  nun::Int         # Vorticity hyperviscous order
end

"""
Constructor that accepts generating function for the topographic height, eta.
"""
function Params(g::TwoDGrid, f0::Float64, beta::Float64, FU::Function,
  eta::Function, mu::Float64, nu::Float64, nun::Int)
  etagrid = eta.(g.X, g.Y)
  Params(f0, beta, FU, eta, mu, nu, nun)
end

# Equations
function Equation(p::Params, g::TwoDGrid)
  LC = -p.mu - p.nu * g.KKrsq.^(0.5*p.nun)
  FourierFlows.Equation{2}(LC, calcN!)
end

# Vars
struct Vars <: AbstractVars
  q::Array{Float64,2}
  U::Float64
  u::Array{Float64,2}
  v::Array{Float64,2}
  uUq::Array{Float64,2}
  vq::Array{Float64,2}
  psi::Array{Float64,2}
  zeta::Array{Float64,2}
  qh::Array{Complex{Float64},2}
  uh::Array{Complex{Float64},2}
  vh::Array{Complex{Float64},2}
  uUqh::Array{Complex{Float64},2}
  vqh::Array{Complex{Float64},2}
  psih::Array{Complex{Float64},2}
  zetah::Array{Complex{Float64},2}
end

function Vars(g::TwoDGrid)
  U     = 0.0
  @createarrays Float64 (g.nx, g.ny) q u v uUq vq psi zeta 
  @createarrays Complex{Float64} (g.nkr, g.nl) qh uh vh uUqh vqh psih zetah
 return Vars(q, U, u, v, uUq, vq, psi, zeta, qh, uh, vh,
    uUqh, vqh, psih, zetah)
end

# Solvers
function calcN!(N::Array{Complex{Float64}, 2}, sol::Array{Complex{Float64}, 2},
  t::Float64, s::State, v::Vars, p::Params, g::TwoDGrid)

  # Note that U = sol[1, 1]. For all other elements ζ = sol
  v.U = sol[1, 1].re
  @. v.zetah = sol
  v.zetah[1, 1] = 0

  @. v.uh =  im * g.l  * g.invKKrsq * v.zetah
  @. v.vh = -im * g.kr * g.invKKrsq * v.zetah

  A_mul_B!(v.zeta, g.irfftplan, v.zetah)
  A_mul_B!(v.u, g.irfftplan, v.uh)
  A_mul_B!(v.v, g.irfftplan, v.vh)

  @. v.q = v.zeta + v.eta
  @. v.uUq = (v.U + v.u)*v.q
  @. v.vq  = v.v*v.q

  A_mul_B!(v.uUqh, g.rfftplan, v.uUq)
  A_mul_B!(v.vqh,  g.rfftplan, v.vq)

  # Nonlinear term for q
  @. N = -im*g.kr*v.uUqh - im*g.l*v.vqh - p.beta*v.vh

  # 'Nonlinear' term for U with topographic correlation.
  # Note: < v*eta > = sum(v*eta)*dx*dy/(Lx*Ly)
  N[1, 1] = p.FU(t) - sum(v.v*v.eta)*g.dx*g.dy/(g.Lx*g.Ly)

end

# Helper functions
function updatevars!(s::State, v::Vars, p::Params, g::TwoDGrid)
  v.U = s.sol[1, 1].re
  @. v.zetah = s.sol
  v.zetah[1, 1] = 0.0

  @. v.psih = -v.zetah * g.invKKrsq
  @. v.uh = -im * g.l  * v.psih
  @. v.vh =  im * g.kr * v.psih

  zetah1 = deepcopy(v.zetah)
  uh1 = deepcopy(v.uh)
  vh1 = deepcopy(v.vh)
  
  A_mul_B!(v.zeta, g.irfftplan, zetah1)
  A_mul_B!(v.psi, g.irfftplan, v.psih)
  A_mul_B!(v.u, g.irfftplan, uh1)
  A_mul_B!(v.v, g.irfftplan, vh1)

  @. v.q = v.zeta + v.eta
  nothing
end

function set_zeta!(s::State, v::AbstractVars, p::AbstractParams, g::TwoDGrid,
  zeta::Array{Float64, 2})

  A_mul_B!(v.zetah, g.rfftplan, zeta)
  v.zetah[1, 1] = 0.0
  @. s.sol = v.zetah

  updatevars!(s, v, p, g)
  nothing
end

function set_zeta!(prob::AbstractProblem, zeta)
  set_zeta!(prob.state, prob.vars, prob.params, prob.grid, zeta)
  nothing
end

end # module
