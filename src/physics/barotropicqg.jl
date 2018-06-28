module BarotropicQG
using FourierFlows
Grid = TwoDGrid

"""
    Problem(; parameters...)

Construct a BarotropicQG turbulence problem.
"""
function Problem(; nx=256, Lx=2π, ny=nx, Ly=Lx, f0 = 1.0, beta=0.0, eta=nothing,
    nu=0.0, nnu=1, mu=0.0, dt=0.01, stepper="RK4", calcFU=nothing, calcFq=nothing)

  # the grid
  g  = BarotropicQG.Grid(nx, Lx, ny, Ly)

  # topographic PV
  if eta==nothing
    eta = 0*g.X
    etah = rfft(eta)
  end

  # construct params depending whether the problem is forced or not
  if calcFU != nothing || calcFq != nothing # the problem is forced
    if calcFU==nothing
      calcFU_fun(t) = 0.0
    else
      calcFU_fun = calcFU
    end
    if calcFq==nothing
      function calcFq!(Fh, sol, t, s, v, p, g)
        nothing
      end
      calcFq_fun = calcFq!
    else
      calcFq_fun = calcFq
    end
    if typeof(eta)!=Array{Float64,2} #this is true if eta was passes in Problem as a function
      pr = BarotropicQG.ForcedParams(g, f0, beta, eta, mu, nu, nnu, calcFU_fun, calcFq_fun)
    else
      pr = BarotropicQG.ForcedParams(f0, beta, eta, etah, mu, nu, nnu, calcFU_fun, calcFq_fun)
    end
    vs = BarotropicQG.ForcedVars(g)
    eq = BarotropicQG.Equation(pr, g)
    ts = FourierFlows.autoconstructtimestepper(stepper, dt, eq.LC, g)
  else # initial value problem
    if typeof(eta)!=Array{Float64,2} #this is true if eta was passes in Problem as a function
      pr = BarotropicQG.Params(g, f0, beta, eta, mu, nu, nnu)
    else
      pr = BarotropicQG.Params(f0, beta, eta, etah, mu, nu, nnu)
    end
    vs = BarotropicQG.Vars(g)
    eq = BarotropicQG.Equation(pr, g)
    ts = FourierFlows.autoconstructtimestepper(stepper, dt, eq.LC, g)
  end
  FourierFlows.Problem(g, vs, pr, eq, ts)
end

InitialValueProblem(; kwargs...) = Problem(; kwargs...)
ForcedProblem(; kwargs...) = Problem(; kwargs...)

"""
    Params(g::TwoDGrid, f0, beta, FU, eta, mu, nu, nun)

Returns the params for an unforced two-dimensional barotropic QG problem.
"""
struct Params{T} <: AbstractParams
  f0::T                      # Constant planetary vorticity
  beta::T                    # Planetary vorticity y-gradient
  eta::Array{T,2}            # Topographic PV
  etah::Array{Complex{T},2}  # FFT of Topographic PV
  mu::T                      # Linear drag
  nu::T                      # Viscosity coefficient
  nun::Int                   # Hyperviscous order (nun=1 is plain old viscosity)
end

"""
    Params(g::TwoDGrid, f0, beta, eta::Function, mu, nu, nun)

Constructor for Params that accepts a generating function for the topographic PV.
"""
function Params(g::TwoDGrid, f0, beta, eta::Function, mu, nu, nun)
  etagrid = eta(g.X, g.Y)
  etah = rfft(etagrid)
  Params(f0, beta, etagrid, etah, mu, nu, nun)
end

"""
    ForcedParams(g::TwoDGrid, f0, beta, FU, eta, mu, nu, nun)

Returns the params for an forced two-dimensional barotropic QG problem.
"""
struct ForcedParams{T} <: AbstractParams
  f0::T                      # Constant planetary vorticity
  beta::T                    # Planetary vorticity y-gradient
  eta::Array{T,2}            # Topographic PV
  etah::Array{Complex{T},2}  # FFT of Topographic PV
  mu::T                      # Linear drag
  nu::T                      # Viscosity coefficient
  nun::Int                   # Hyperviscous order (nun=1 is plain old viscosity)
  calcFU::Function   # Function that calculates the forcing F(t) on
                      # domain-averaged zonal flow U(t)
  calcFq!::Function   # Function that calculates the forcing on QGPV q
end

"""
    ForcedParams(g::TwoDGrid, f0, beta, eta::Function, mu, nu, nun, calcFU, calcFq)

Constructor for Params that accepts a generating function for the topographic PV.
"""
function ForcedParams(g::TwoDGrid, f0, beta, eta::Function, mu, nu, nun, calcFU::Function, calcFq::Function)
  etagrid = eta(g.X, g.Y)
  etah = rfft(etagrid)
  ForcedParams(f0, beta, etagrid, etah, mu, nu, nun, calcFU, calcFq)
end


"""
    Equation(p, g)

Returns the equation for two-dimensional barotropic QG problem with params p and grid g.
"""
function Equation(p::Params, g)
  LC = @. -p.mu - p.nu*g.KKrsq^p.nun + im*p.beta*g.kr*g.invKKrsq
  LC[1, 1] = 0
  FourierFlows.Equation{typeof(LC[1, 1]),2}(LC, calcN_advection!)
end

function Equation(p::ForcedParams, g)
  LC = @. -p.mu - p.nu*g.KKrsq^p.nun + im*p.beta*g.kr*g.invKKrsq
  LC[1, 1] = 0
  FourierFlows.Equation{typeof(LC[1, 1]),2}(LC, calcN_forced!)
end

"""
    Vars(g)

Returns the vars for unforced two-dimensional barotropic QG problem with grid g.
"""
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

"""
    ForcedVars(g)

Returns the vars for forced two-dimensional barotropic QG problem with grid g.
"""
mutable struct ForcedVars{T} <: AbstractVars
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
  Fqh::Array{Complex{T},2}
  prevsol::Array{Complex{T},2}
end

function ForcedVars(g::TwoDGrid)
  T = typeof(g.Lx)
  @createarrays T (g.nx, g.ny) q u v uUq vq psi zeta
  @createarrays Complex{T} (g.nkr, g.nl) qh uh vh uUqh vqh psih zetah Fqh prevsol
  ForcedVars(q, 0.0, u, v, uUq, vq, psi, zeta, qh, uh, vh, uUqh, vqh, psih, zetah, Fqh, prevsol)
end


# -------
# Solvers
# -------

function calcN_advection!(N, sol, t, s, v, p, g)
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

  # Nonlinear advection term for q
  @. N = -im*g.kr*v.uUqh - im*g.l*v.vqh
end

function calcN_forced!(N, sol, t, s, v, p, g)
  calcN_advection!(N, sol, t, s, v, p, g)

  # 'Nonlinear' term for U with topographic correlation.
  # Note: < v*eta > = sum( conj(vh)*eta ) / (nx^2*ny^2) if fft is used
  # while < v*eta > = 2*sum( conj(vh)*eta ) / (nx^2*ny^2) if rfft is used
  if size(sol)[1] == g.nkr
    N[1, 1] = p.calcFU(t) + 2*sum(conj(v.vh).*p.etah).re / (g.nx^2.0*g.ny^2.0)
  else
    N[1, 1] = p.calcFU(t) + sum(conj(v.vh).*p.etah).re / (g.nx^2.0*g.ny^2.0)
  end
  if t == s.t # not a substep
    v.prevsol .= s.sol # used to compute budgets when forcing is stochastic
    p.calcFq!(v.Fqh, sol, t, s, v, p, g)
  end
  @. N += v.Fqh
  nothing
end


# ----------------
# Helper functions
# ----------------

"""
    updatevars!(v, s, g)

Update the vars in v on the grid g with the solution in s.sol.
"""
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

"""
    set_zeta!(prob, zeta)
    set_q!(s, v, g, zeta)

Set the solution s.sol as the transform of zeta and update variables v
on the grid g.
"""
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
U00(prob) = real(s.sol[1, 1])
energy00(prob) = real(0.5*prob.state.sol[1, 1].^2)
enstrophy00(prob) = real(prob.params.beta*prob.state.sol[1, 1])



end # module
