module BarotropicQG
using FourierFlows, FFTW
import LinearAlgebra: mul!, ldiv!

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
    if calcFU == nothing
      calcFU_fun(t) = nothing
    else
      calcFU_fun = calcFU
    end
    if calcFq==nothing
      function calcFq!(Fqh, sol, t, s, v, p, g)
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
      etah = rfft(eta)
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
    Params(g::TwoDGrid, f0, beta, FU, eta, mu, nu, nnu)

Returns the params for an unforced two-dimensional barotropic QG problem.
"""
struct Params{T} <: AbstractParams
  f0::T                      # Constant planetary vorticity
  beta::T                    # Planetary vorticity y-gradient
  eta::Array{T,2}            # Topographic PV
  etah::Array{Complex{T},2}  # FFT of Topographic PV
  mu::T                      # Linear drag
  nu::T                      # Viscosity coefficient
  nnu::Int                   # Hyperviscous order (nnu=1 is plain old viscosity)
end

"""
    Params(g::TwoDGrid, f0, beta, eta::Function, mu, nu, nnu)

Constructor for Params that accepts a generating function for the topographic PV.
"""
function Params(g::TwoDGrid, f0, beta, eta::Function, mu, nu, nnu)
  etagrid = eta(g.X, g.Y)
  etah = rfft(etagrid)
  Params(f0, beta, etagrid, etah, mu, nu, nnu)
end

"""
    ForcedParams(g::TwoDGrid, f0, beta, FU, eta, mu, nu, nnu)

Returns the params for an forced two-dimensional barotropic QG problem.
"""
struct ForcedParams{T} <: AbstractParams
  f0::T                      # Constant planetary vorticity
  beta::T                    # Planetary vorticity y-gradient
  eta::Array{T,2}            # Topographic PV
  etah::Array{Complex{T},2}  # FFT of Topographic PV
  mu::T                      # Linear drag
  nu::T                      # Viscosity coefficient
  nnu::Int                   # Hyperviscous order (nnu=1 is plain old viscosity)
  calcFU::Function    # Function that calculates the forcing F(t) on
                      # domain-averaged zonal flow U(t)
  calcFq!::Function   # Function that calculates the forcing on QGPV q
end

"""
    ForcedParams(g::TwoDGrid, f0, beta, eta::Function, mu, nu, nnu, calcFU, calcFq)

Constructor for Params that accepts a generating function for the topographic PV.
"""
function ForcedParams(g::TwoDGrid, f0, beta, eta::Function, mu, nu, nnu, calcFU::Function, calcFq::Function)
  etagrid = eta(g.X, g.Y)
  etah = rfft(etagrid)
  ForcedParams(f0, beta, etagrid, etah, mu, nu, nnu, calcFU, calcFq)
end


"""
    Equation(p, g)

Returns the equation for two-dimensional barotropic QG problem with params p and grid g.
"""
function Equation(p::Params, g)
  LC = @. -p.mu - p.nu*g.KKrsq^p.nnu + im*p.beta*g.kr*g.invKKrsq
  LC[1, 1] = 0
  FourierFlows.Equation{typeof(LC[1, 1]),2}(LC, calcN_advection!)
end

function Equation(p::ForcedParams, g)
  LC = @. -p.mu - p.nu*g.KKrsq^p.nnu + im*p.beta*g.kr*g.invKKrsq
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
  psi::Array{T,2}
  zeta::Array{T,2}
  qh::Array{Complex{T},2}
  uh::Array{Complex{T},2}
  vh::Array{Complex{T},2}
  psih::Array{Complex{T},2}
  zetah::Array{Complex{T},2}
end

function Vars(g::TwoDGrid)
  T = typeof(g.Lx)
  @createarrays T (g.nx, g.ny) q u v psi zeta
  @createarrays Complex{T} (g.nkr, g.nl) qh uh vh psih zetah
  Vars(q, 0.0, u, v, psi, zeta, qh, uh, vh, psih, zetah)
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
  psi::Array{T,2}
  zeta::Array{T,2}
  qh::Array{Complex{T},2}
  uh::Array{Complex{T},2}
  vh::Array{Complex{T},2}
  psih::Array{Complex{T},2}
  zetah::Array{Complex{T},2}
  Fqh::Array{Complex{T},2}
  prevsol::Array{Complex{T},2}
end

function ForcedVars(g::TwoDGrid)
  T = typeof(g.Lx)
  @createarrays T (g.nx, g.ny) q u v psi zeta
  @createarrays Complex{T} (g.nkr, g.nl) qh uh vh psih zetah Fqh prevsol
  ForcedVars(q, 0.0, u, v, psi, zeta, qh, uh, vh, psih, zetah, Fqh, prevsol)
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

  ldiv!(v.zeta, g.rfftplan, v.zetah)
  ldiv!(v.u, g.rfftplan, v.uh)
  v.psih = deepcopy(v.vh) # FFTW's irfft destroys its input; v.vh is needed for N[1, 1]
  ldiv!(v.v, g.rfftplan, v.psih)

  @. v.q = v.zeta + p.eta
  @. v.u = (v.U + v.u)*v.q # (U+u)*q
  @. v.v = v.v*v.q # v*q

  mul!(v.uh, g.rfftplan, v.u) # \hat{(u+U)*q}
  # Nonlinear advection term for q (part 1)
  @. N = -im*g.kr*v.uh # -∂[(U+u)q]/∂x
  mul!(v.uh, g.rfftplan, v.v) # \hat{v*q}
  @. N += - im*g.l*v.uh # -∂[vq]/∂y
end

function calcN_forced!(N, sol, t, s, v, p, g)
  calcN_advection!(N, sol, t, s, v, p, g)
  if p.calcFU(t) != nothing
    # 'Nonlinear' term for U with topographic correlation.
    # Note: < v*eta > = sum( conj(vh)*eta ) / (nx^2*ny^2) if fft is used
    # while < v*eta > = 2*sum( conj(vh)*eta ) / (nx^2*ny^2) if rfft is used
    N[1, 1] = p.calcFU(t) + 2*sum(conj(v.vh).*p.etah).re / (g.nx^2.0*g.ny^2.0)
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

  ldiv!(v.zeta, g.rfftplan, zetah1)
  ldiv!(v.psi, g.rfftplan, v.psih)
  ldiv!(v.u, g.rfftplan, uh1)
  ldiv!(v.v, g.rfftplan, vh1)

  @. v.q = v.zeta + p.eta
  nothing
end

updatevars!(prob) = updatevars!(prob.state, prob.vars, prob.params, prob.grid)

"""
    set_zeta!(prob, zeta)
    set_zeta!(s, v, g, zeta)

Set the solution s.sol as the transform of zeta and update variables v
on the grid g.
"""
function set_zeta!(s, v::Vars, p, g, zeta)
  mul!(v.zetah, g.rfftplan, zeta)
  v.zetah[1, 1] = 0.0
  @. s.sol = v.zetah

  updatevars!(s, v, p, g)
  nothing
end

function set_zeta!(s, v::ForcedVars, p, g, zeta)
  v.U = deepcopy(s.sol[1, 1])
  mul!(v.zetah, g.rfftplan, zeta)
  v.zetah[1, 1] = 0.0
  @. s.sol = v.zetah
  s.sol[1, 1] = v.U

  updatevars!(s, v, p, g)
  nothing
end

set_zeta!(prob::AbstractProblem, zeta) = set_zeta!(prob.state, prob.vars, prob.params, prob.grid, zeta)

"""
    set_U!(prob, U)
    set_U!(s, v, g, U)

Set the (kx,ky)=(0,0) part of solution s.sol as the domain-average zonal flow U.
"""
function set_U!(s, v, p, g, U::Float64)
  s.sol[1, 1] = U
  updatevars!(s, v, p, g)
  nothing
end

set_U!(prob::AbstractProblem, U::Float64) = set_U!(prob.state, prob.vars, prob.params, prob.grid, U)


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
  s, v, g = prob.state, prob.vars, prob.grid
  @. v.uh = s.sol
  v.uh[1, 1] = 0
  0.5*FourierFlows.parsevalsum2(v.uh, g)/(g.Lx*g.Ly)
end

"""
Returns the energy of the domain-averaged U.
"""
energy00(prob) = real(0.5*prob.state.sol[1, 1].^2)

"""
Returns the enstrophy of the domain-averaged U.
"""
enstrophy00(prob) = real(prob.params.beta*prob.state.sol[1, 1])

"""
    dissipation(prob)
    dissipation(s, v, p, g)

Returns the domain-averaged dissipation rate. nnu must be >= 1.
"""
@inline function dissipation(s, v, p, g)
  @. v.uh = g.KKrsq^(p.nnu-1) * abs2(s.sol)
  v.uh[1, 1] = 0
  p.nu/(g.Lx*g.Ly)*FourierFlows.parsevalsum(v.uh, g)
end

@inline dissipation(prob::AbstractProblem) = dissipation(prob.state, prob.vars, prob.params, prob.grid)

"""
    work(prob)
    work(s, v, p, g)

Returns the domain-averaged rate of work of energy by the forcing Fqh.
"""
@inline function work(s, v::ForcedVars, g)
  @. v.uh = g.invKKrsq * (v.prevsol + s.sol)/2.0 * conj(v.Fqh) # Stratonovich
  # @. v.uh = g.invKKrsq * v.prevsol * conj(v.Fqh)             # Ito
  1/(g.Lx*g.Ly)*FourierFlows.parsevalsum(v.uh, g)
end

@inline work(prob::AbstractProblem) = work(prob.state, prob.vars, prob.grid)

"""
    drag(prob)
    drag(s, v, p, g)

Returns the extraction of domain-averaged energy by drag mu.
"""
@inline function drag(s, v, p, g)
  @. v.uh = g.KKrsq^(-1) * abs2(s.sol)
  v.uh[1, 1] = 0
  p.mu/(g.Lx*g.Ly)*FourierFlows.parsevalsum(v.uh, g)
end

@inline drag(prob::AbstractProblem) = drag(prob.state, prob.vars, prob.params, prob.grid)


end # module
