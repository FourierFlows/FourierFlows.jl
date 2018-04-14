module BarotropicQG
using FourierFlows
Grid = TwoDGrid

# Params
struct Params{T} <: AbstractParams
  f0::T                      # Constant planetary vorticity
  beta::T                    # Planetary vorticity y-gradient
  FU::Function               # Time-dependent forcing of domain average flow
  eta::Array{T,2}            # Topographic PV
  etah::Array{Complex{T},2}  # FFT of Topographic PV
  mu::T                      # Linear drag
  nu::T                      # Viscosity coefficient
  nun::Int                   # Hyperviscous order (nun=1 is plain old viscosity)
end

"""
Constructor that accepts generating function for the topographic height, eta.
"""
function Params(g::TwoDGrid, f0, beta, FU, eta, mu, nu, nun)
  etagrid = eta(g.X, g.Y)
  etah = rfft(etagrid)
  Params(f0, beta, FU, etagrid, etah, mu, nu, nun)
end

# Equations
function Equation(p, g)
  LC = @. -p.mu - p.nu*g.KKrsq^p.nun + im*p.beta*g.kr*g.invKKrsq
  FourierFlows.Equation(LC, calcN!)
end

# Vars
mutable struct Vars{T} <: AbstractVars
  q::Array{T,2}
  U::T
  u::Array{T,2}
  v::Array{T,2}
  uUq::Array{T,2}
  vq::Array{T,2}
  psi::Array{T,2}
  zeta::Array{T,2}
  qh::Array{Complex{T},2}
  uh::Array{Complex{T},2}
  vh::Array{Complex{T},2}
  uUqh::Array{Complex{T},2}
  vqh::Array{Complex{T},2}
  psih::Array{Complex{T},2}
  zetah::Array{Complex{T},2}
end

function Vars(g::TwoDGrid)
  T = typeof(g.Lx)
  @createarrays T (g.nx, g.ny) q u v uUq vq psi zeta
  @createarrays Complex{T} (g.nkr, g.nl) qh uh vh uUqh vqh psih zetah
 Vars(q, 0.0, u, v, uUq, vq, psi, zeta, qh, uh, vh, uUqh, vqh, psih, zetah)
end

# Solvers
function calcN!(N, sol, t, s, v, p, g)
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
function updatevars!(s, v, p, g)
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

updatevars!(prob) = updatevars!(prob.state, prob.vars, prob.params, prob.grid)


function set_zeta!(s, v, p, g, zeta)
  #zeta = similar(v.u)

  A_mul_B!(v.zetah, g.rfftplan, zeta)
  v.zetah[1, 1] = 0.0
  @. s.sol = v.zetah

  updatevars!(s, v, p, g)
  nothing
end

set_zeta!(prob::AbstractProblem, zeta) = set_zeta!(prob.state, prob.vars, prob.params, prob.grid, zeta)


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
function enstrophy(prob)
  s, g = prob.state, prob.grid
  0.5*FourierFlows.parsevalsum2(s.sol, g)/(g.Lx*g.Ly)
end

"""
Returns the domain-averaged enstrophy.
"""
U00(prob) = s.sol[1, 1]
energy00(prob) = 0.5*prob.state.sol[1, 1].^2
enstrophy00(prob) = prob.params.beta*prob.state.sol[1, 1]



end # module
