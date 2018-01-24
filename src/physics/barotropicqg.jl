module BarotropicQG
using FourierFlows
Grid = TwoDGrid

# Params
struct Params <: AbstractParams
  f0::Float64      # Constant planetary vorticity
  beta::Float64    # Planetary vorticity y-gradient
  FU::Function     # Time-dependent forcing of domain average flow
  eta::Array{Float64,2}            # Topographic PV
  etah::Array{Complex{Float64},2}  # FFT of Topographic PV
  μ::Float64       # Linear drag
  ν::Float64       # Viscosity coefficient
  νn::Int          # Hyperviscous order (νn=1 is plain old viscosity)
end

"""
Constructor that accepts generating function for the topographic height, eta.
"""
function Params(g::TwoDGrid, f0::Float64, beta::Float64, FU::Function,
  eta::Function, μ::Float64, ν::Float64, νn::Int)
  etagrid = eta(g.X, g.Y)
  etah = rfft(etagrid)
  Params(f0, beta, FU, etagrid, etah, μ, ν, νn)
end

# Equations
function Equation(p::Params, g::TwoDGrid)
  LC = -p.μ - p.ν * g.KKrsq.^p.νn + im*p.beta*g.Kr.*g.invKKrsq
  FourierFlows.Equation{2}(LC, calcN!)
end

# Vars
mutable struct Vars <: AbstractVars
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
  vh = deepcopy(v.vh)   # FFTW's irfft destroys its input; v.vh is needed for N
  A_mul_B!(v.v, g.irfftplan, vh)


  @. v.q = v.zeta + p.eta
  @. v.uUq = (v.U + v.u)*v.q
  @. v.vq  = v.v*v.q

  A_mul_B!(v.uUqh, g.rfftplan, v.uUq)
  A_mul_B!(v.vqh,  g.rfftplan, v.vq)

  # Nonlinear term for q
  @. N = -im*g.kr*v.uUqh - im*g.l*v.vqh

  # 'Nonlinear' term for U with topographic correlation.
  # Note: < v*eta > = sum( conj(vh)*eta ) / (nx^2*ny^2) if fft is used
  # while < v*eta > = 2*sum( conj(vh)*eta ) / (nx^2*ny^2) if rfft is used
  if size(sol)[1] == g.nkr
    N[1, 1] = p.FU(t) + 2*sum(conj(v.vh).*p.etah).re / (g.nx^2.0*g.ny^2.0)
  else
    N[1, 1] = p.FU(t) + sum(conj(v.vh).*p.etah).re / (g.nx^2.0*g.ny^2.0)
  end

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

  @. v.q = v.zeta + p.eta
  nothing
end

function updatevars!(prob::AbstractProblem)
  s, v, p, g = prob.state, prob.vars, prob.params, prob.grid
  updatevars!(s, v, p, g)
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


"""
Calculate the domain-averaged kinetic energy.
"""
function energy(prob::AbstractProblem)
  s, g = prob.state, prob.grid
  0.5*(FourierFlows.parsevalsum2(g.Kr.*g.invKKrsq.*s.sol, g)
        + FourierFlows.parsevalsum2(g.Lr.*g.invKKrsq.*s.sol, g))/(g.Lx*g.Ly)
end


"""
Returns the domain-averaged enstrophy.
"""
function enstrophy(prob::AbstractProblem)
  s, g = prob.state, prob.grid
  0.5*FourierFlows.parsevalsum2(s.sol, g)/(g.Lx*g.Ly)
end


end # module
