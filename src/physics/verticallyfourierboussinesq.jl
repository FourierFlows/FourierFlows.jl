module VerticallyFourierBoussinesq
using FourierFlows, FFTW

import LinearAlgebra: mul!, ldiv!
import FourierFlows: getfieldspecs, structvarsexpr, parsevalsum, parsevalsum2, autoconstructtimestepper

# -------
# Problem
# -------

"""
    Problem(; parameters...)

Construct a VerticallyFourierBoussinesq initial value problem.
"""
function Problem(;
  # Numerical parameters
    nx = 128,
    Lx = 2π,
    ny = nx,
    Ly = Lx,
    dt = 0.01,
  # Drag and/or hyper-/hypo-viscosity
   nu0 = 0,
  nnu0 = 1,
   nu1 = 0,
  nnu1 = 1,
   mu0 = 0,
  nmu0 = 0,
   mu1 = 0,
  nmu1 = 0,
  # Physical parameters
     f = 1,
     N = 10,
     m = 40,
  # Optional uniform and steady background flow
    Ub = 0,
    Vb = 0,
  # Timestepper and eqn options
  stepper = "RK4",
  calcF = nothing,
  nthreads = Sys.CPU_THREADS,
  T = Float64
  )

  g  = TwoDGrid(nx, Lx, ny, Ly)
  pr = Params{T}(nu0, nnu0, nu1, nnu1, mu0, nmu0, mu1, nmu1, f, N, m, Ub, Vb)
  vs = Vars(g)
  eq = Equation(pr, g)
  ts = autoconstructtimestepper(stepper, dt, eq.LCc, eq.LCr, g)

  FourierFlows.Problem(g, vs, pr, eq, ts)
end


# ------
# Params
# ------

abstract type TwoModeParams <: AbstractParams end

"""
    Params(nu0, nnu0, nu1, nnu1, f, N, m, Ub, Vb)

Construct parameters for the Two-Fourier-mode Boussinesq problem. Suffix 0
refers to zeroth mode; 1 to first mode. f, N, m are Coriolis frequency,
buoyancy frequency, and vertical wavenumber of the first mode, respectively.
The optional constant background velocity (Ub,Vb) is set to zero by default.
The viscosity is applied only to the first-mode horizontal velocities.
"""
struct Params{T} <: TwoModeParams
  nu0::T     # Mode-0 viscosity
  nnu0::Int  # Mode-0 hyperviscous order
  nu1::T     # Mode-1 viscosity
  nnu1::Int  # Mode-1 hyperviscous order
  mu0::T     # Mode-0 drag / hypoviscosity
  nmu0::Int  # Mode-0 drag / hypoviscous order
  mu1::T     # Mode-1 drag / hypoviscosity
  nmu1::Int  # Mode-1 drag / hypoviscous order
  f::T       # Planetary vorticity
  N::T       # Buoyancy frequency
  m::T       # Mode-one wavenumber
  Ub::T      # Steady mode-0 mean x-velocity
  Vb::T      # Steady mode-0 mean y-velocity
end

Params(nu0, nnu0, nu1, nnu1, f, N, m) = Params(nu0, nnu0, nu1, nnu1, f, N, m, Ub, Vb)


# ---------
# Equations
# ---------

function Equation(p::TwoModeParams, g::TwoDGrid)
  LCc, LCr = getlinearcoefficients(p, g)
  DualEquation(LCc, LCr, calcN!)
end

function getlinearcoefficients(p::TwoModeParams, g::TwoDGrid)
  LCr = @. -p.nu0*g.KKrsq^p.nnu0 - p.mu0*g.KKrsq^p.nmu0

  LCc = zeros(g.nk, g.nl, 3)
  LCc[:, :, 1] = @. -p.nu1*g.KKsq^p.nnu1 - p.mu1*g.KKsq^p.nmu1
  LCc[:, :, 2] = @. -p.nu1*g.KKsq^p.nnu1 - p.mu1*g.KKsq^p.nmu1

  LCr[1, 1] = 0
  LCc[1 ,1, 1] = 0
  LCc[1 ,1, 2] = 0
  LCc, LCr
end


# ----
# Vars
# ----

abstract type VerticallyFourierVars <: AbstractVars end

physicalvarsr = [:Z, :U, :V, :UZuzvw, :VZvzuw, :Ux, :Uy, :Vx, :Vy, :Psi]
physicalvarsc = [:u, :v, :w, :p, :zeta, :Uu, :Uv, :Up, :Vu, :Vv, :Vp, :uUxvUy, :uVxvVy]
transformvarsr = [ Symbol(var, :h) for var in physicalvarsr ]
transformvarsc = [ Symbol(var, :h) for var in physicalvarsc ]

varspecs = cat(
  getfieldspecs(physicalvarsr, :(Array{T,2})),
  getfieldspecs(physicalvarsc, :(Array{Complex{T},2})),
  getfieldspecs(transformvarsr, :(Array{Complex{T},2})),
  getfieldspecs(transformvarsc, :(Array{Complex{T},2})),
  dims=1)

eval(structvarsexpr(:Vars, varspecs; parent=:VerticallyFourierVars))

function Vars(g)
  @createarrays typeof(g.Lx) (g.nx, g.ny) Z U V UZuzvw VZvzuw Ux Uy Vx Vy Psi
  @createarrays Complex{typeof(g.Lx)} (g.nx, g.ny) u v w p zeta Uu Uv Up Vu Vv Vp uUxvUy uVxvVy
  @createarrays Complex{typeof(g.Lx)} (g.nkr, g.nl) Zh Uh Vh UZuzvwh VZvzuwh Uxh Uyh Vxh Vyh Psih
  @createarrays Complex{typeof(g.Lx)} (g.nk, g.nl) uh vh wh ph zetah Uuh Uvh Uph Vuh Vvh Vph uUxvUyh uVxvVyh
  Vars(Z, U, V, UZuzvw, VZvzuw, Ux, Uy, Vx, Vy, Psi, u, v, w, p, zeta, Uu, Uv, Up, Vu, Vv, Vp, uUxvUy, uVxvVy,
       Zh, Uh, Vh, UZuzvwh, VZvzuwh, Uxh, Uyh, Vxh, Vyh, Psih,
       uh, vh, wh, ph, zetah, Uuh, Uvh, Uph, Vuh, Vvh, Vph, uUxvUyh, uVxvVyh)
end


# -------
# Solvers
# -------

function calcN_linearterms!(Nc, Nr, solc, solr, t, s, v, p, g)
  @views @. Nc[:, :, 1] =  p.f*solc[:, :, 2] - im*g.k*solc[:, :, 3]               # u
  @views @. Nc[:, :, 2] = -p.f*solc[:, :, 1] - im*g.l*solc[:, :, 3]               # v
  @views @. Nc[:, :, 3] = -im*p.N^2/p.m^2*(g.k*solc[:, :, 1] + g.l*solc[:, :, 2]) # p
  nothing
end

function calcN!(Nc, Nr, solc, solr, t, s, v, p, g)
  v.Zh .= solr
  @. v.Psih = -g.invKKrsq*v.Zh

  @. v.Uh  = -im*g.l  * v.Psih
  @. v.Vh  =  im*g.kr * v.Psih
  @. v.Uxh =  im*g.kr * v.Uh
  @. v.Uyh =  im*g.l  * v.Uh
  @. v.Vxh =  im*g.kr * v.Vh
  @. v.Vyh =  im*g.l  * v.Vh

  v.Uh[1, 1] += p.Ub*g.nx*g.ny
  v.Vh[1, 1] += p.Vb*g.nx*g.ny

  # Inverse transforms
  ldiv!(v.Z,  g.rfftplan, v.Zh)
  ldiv!(v.U,  g.rfftplan, v.Uh)
  ldiv!(v.V,  g.rfftplan, v.Vh)
  ldiv!(v.Ux, g.rfftplan, v.Uxh)
  ldiv!(v.Uy, g.rfftplan, v.Uyh)
  ldiv!(v.Vx, g.rfftplan, v.Vxh)
  ldiv!(v.Vy, g.rfftplan, v.Vyh)

  @views ldiv!(v.u, g.fftplan, solc[:, :, 1])
  @views ldiv!(v.v, g.fftplan, solc[:, :, 2])
  @views ldiv!(v.p, g.fftplan, solc[:, :, 3])

  @views @. v.zetah = im*g.k*solc[:, :, 2] - im*g.l*solc[:, :, 1]
  @views @. v.wh = -(g.k*solc[:, :, 1] + g.l*solc[:, :, 2]) / p.m

  ldiv!(v.w, g.fftplan, v.wh)
  ldiv!(v.zeta, g.fftplan, v.zetah)

  # Multiplies
  @. v.UZuzvw = (v.U * v.Z
    + real(   v.u*conj(v.zeta)  +  im*p.m*v.v*conj(v.w)
            + conj(v.u)*v.zeta  -  im*p.m*conj(v.v)*v.w   ))

  @. v.VZvzuw = (v.V * v.Z
    + real(   v.v*conj(v.zeta)  -  im*p.m*v.u*conj(v.w)
            + conj(v.v)*v.zeta  +  im*p.m*conj(v.u)*v.w   ))

  @. v.uUxvUy = v.u*v.Ux + v.v*v.Uy
  @. v.uVxvVy = v.u*v.Vx + v.v*v.Vy

  @. v.Uu = v.U * v.u
  @. v.Vu = v.V * v.u
  @. v.Uv = v.U * v.v
  @. v.Vv = v.V * v.v
  @. v.Up = v.U * v.p
  @. v.Vp = v.V * v.p

  # Forward transforms
  mul!(v.UZuzvwh, g.rfftplan, v.UZuzvw)
  mul!(v.VZvzuwh, g.rfftplan, v.VZvzuw)
  mul!(v.uUxvUyh, g.fftplan, v.uUxvUy)
  mul!(v.uVxvVyh, g.fftplan, v.uVxvVy)

  mul!(v.Uuh, g.fftplan, v.Uu)
  mul!(v.Uvh, g.fftplan, v.Uv)
  mul!(v.Vuh, g.fftplan, v.Vu)
  mul!(v.Vvh, g.fftplan, v.Vv)
  mul!(v.Uph, g.fftplan, v.Up)
  mul!(v.Vph, g.fftplan, v.Vp)

  # Zeroth-mode nonlinear term
  calcN_linearterms!(Nc, Nr, solc, solr, t, s, v, p, g)
  @. Nr = - im*g.kr*v.UZuzvwh - im*g.l*v.VZvzuwh
  @views @. Nc[:, :, 1] += -im*g.k*v.Uuh - im*g.l*v.Vuh - v.uUxvUyh # u
  @views @. Nc[:, :, 2] += -im*g.l*v.Vvh - im*g.k*v.Uvh - v.uVxvVyh # v
  @views @. Nc[:, :, 3] += -im*g.k*v.Uph - im*g.l*v.Vph             # p
  nothing
end


# ----------------
# Helper functions
# ----------------

"""
    updatevars!(prob)

Update variables to correspond to the solution in s.sol or prob.state.sol.
"""
function updatevars!(v, s, p, g)
  v.Zh .= s.solr
  @. v.Psih = -g.invKKrsq*v.Zh
  @. v.Uh   = -im*g.l*v.Psih
  @. v.Vh   =  im*g.kr*v.Psih

  Psih1 = deepcopy(v.Psih)
  Zh1 = deepcopy(v.Zh)
  Uh1 = deepcopy(v.Uh)
  Vh1 = deepcopy(v.Vh)

  ldiv!(v.Psi, g.rfftplan, Psih1)
  ldiv!(v.Z, g.rfftplan, Zh1)
  ldiv!(v.U, g.rfftplan, Uh1)
  ldiv!(v.V, g.rfftplan, Vh1)

  @views v.uh .= s.solc[:, :, 1]
  @views v.vh .= s.solc[:, :, 2]
  @views v.ph .= s.solc[:, :, 3]

  @. v.wh = -1.0/p.m*(g.k*v.uh + g.l*v.vh)

  ldiv!(v.u, g.fftplan, v.uh)
  ldiv!(v.v, g.fftplan, v.vh)
  ldiv!(v.p, g.fftplan, v.ph)
  ldiv!(v.w, g.fftplan, v.wh)
  nothing
end
updatevars!(prob) = updatevars!(prob.vars, prob.state, prob.params, prob.grid)

"""
    set_Z!(prob, Z)

Set zeroth mode vorticity and update variables.
"""
function set_Z!(s, v, p, g, Z)
  mul!(s.solr, g.rfftplan, Z)
  updatevars!(v, s, p, g)
  nothing
end
set_Z!(prob, Z) = set_Z!(prob.state, prob.vars, prob.params, prob.grid, Z)

"""
    set_uvp!(prob, u, v, p)

Set first mode u, v, and p and update variables.
"""
function set_uvp!(s, vs, pr, g, u, v, p)
  uh = fft(u)
  vh = fft(v)
  ph = fft(p)

  @. s.solc[:, :, 1] = uh
  @. s.solc[:, :, 2] = vh
  @. s.solc[:, :, 3] = ph

  updatevars!(vs, s, pr, g)
  nothing
end
set_uvp!(prob, u, v, p) = set_uvp!(prob.state, prob.vars, prob.params, prob.grid, u, v, p)


"""
    set_planewave!(prob, uw, nkw)

Set a plane wave solution with initial speed uw and non-dimensional wave
number nkw. The dimensional wavenumber will be 2π*nkw/Lx.
"""
function set_planewave!(s, vs, pr, g, uw, κ, θ=0; envelope=nothing)
  k = 2π/g.Lx*round(Int, κ*cos(θ))
  l = 2π/g.Lx*round(Int, κ*sin(θ))
  x, y = g.X, g.Y

  # Wave parameters
  f, N, m = pr.f, pr.N, pr.m
  σ = sqrt(f^2 + N^2*(k^2 + l^2)/m^2)

  u0 = uw/2
  v0 = u0 * (σ*l - im*f*k)/(σ*k + im*f*l)
  p0 = u0 * (σ^2 - f^2)/(σ*k + im*f*l)

  Φ = @. k*x + l*y
  u = u0 * exp.(im*Φ)
  v = v0 * exp.(im*Φ)
  p = p0 * exp.(im*Φ)

  if envelope != nothing
    @. u *= envelope(x, y)
    @. v *= envelope(x, y)
    @. p *= envelope(x, y)
  end

  set_uvp!(s, vs, pr, g, u, v, p)
  nothing
end
set_planewave!(prob, uw, nkw; kwargs...) = set_planewave!(
  prob.state, prob.vars, prob.params, prob.grid, uw::Real, nkw::Int; kwargs...)

"""
    set_isotropicwavefield!(prob, amplitude; KE=1.0, maxspeed=nothing)

Generate an isotropic spectrum of waves with an envelope given by the function
amplitude(k, l), and either total kinetic energy KE or maximum speed maxspeed.
"""
function set_isotropicwavefield!(s, vs, pr, g, amplitude; KE=1.0, maxspeed=nothing)
  f, N, m = pr.f, pr.N, pr.m # for clarity
  @createarrays Complex{Float64} (g.nx, g.ny) phase u0 v0 p0

  # Sum Fourier components
  for k in real.(g.k), l in real.(g.l)
    if amplitude(k, l) > 1e-15
      σ = sqrt(f^2 + N^2/m^2*(k^2 + l^2))   # dispersion relation
      phase .= k*g.X .+ l*g.Y .+ 2π*rand()  # random phases
      # Polarization relations
      u0 .+= amplitude(k, l)*exp.(im*phase)
      v0 .+= -u0*(im*f/σ - k*l*N^2/(σ*m)^2)/(1 - (l*N)^2/(σ*m)^2)
      p0 .+= N^2/(σ*m^2) * (k*u0 .+ l*v0)
    end
  end

  if maxspeed == nothing # Normalize by kinetic energy
    uh, vh = fft(u0), fft(v0)
    ke = mode1ke(uh, vh, g)
    norm = sqrt(KE)/sqrt(ke/(g.Lx*g.Ly))
  else
    norm = maxspeed / maximum(
      sqrt.(real.(u0+conj.(u0)).^2 + real.(v0+conj.(v0)).^2))
  end

  u0 .*= norm
  v0 .*= norm
  p0 .*= norm

  set_uvp!(s, vs, pr, g, u0, v0, p0)
  nothing
end
set_isotropicwavefield!(prob::AbstractProblem, amplitude; kwargs...) = set_isotropicwavefield!(prob.state, prob.vars,
  prob.params, prob.grid, amplitude; kwargs...)


# -----------
# Diagnostics
# -----------

"""
    mode0energy(prob)

Returns the domain-averaged energy in the zeroth mode.
"""
@inline function mode0energy(s, v, g)
  @. v.Uh = g.invKKrsq * abs2(s.solr) # qh*Psih
  1/(2*g.Lx*g.Ly)*parsevalsum(v.Uh, g)
end
@inline mode0energy(prob) = mode0energy(prob.state, prob.vars, prob.grid)

"""
    mode1ke(prob)

Returns the domain-averaged kinetic energy in the first mode.
"""
@inline mode1ke(uh, vh, g) = parsevalsum2(uh, g) + parsevalsum2(vh, g)/(g.Lx*g.Ly)
@inline mode1ke(s, g) = @views mode1ke(s.solc[:, :, 1], s.solc[:, :, 2], g)
@inline mode1ke(prob) = mode1ke(prob.state, prob.grid)

"""
    mode1pe(prob)

Returns the domain-averaged potential energy in the first mode.
"""
@inline mode1pe(s, p, g) = @views p.m^2/(g.Lx*g.Ly*p.N^2)*parsevalsum2(s.solc[:, :, 3], g)
@inline mode1pe(prob) = mode1pe(prob.state, prob.params, prob.grid)

"""
    mode1energy(prob)

Returns the domain-averaged total energy in the first mode.
"""
mode1energy(s, p, g) = mode1ke(s, g) + mode1pe(s, p, g)
mode1energy(prob) = mode1energy(prob.state, prob.params, prob.grid)

"""
    mode0dissipation(prob)

Returns the domain-averaged kinetic energy dissipation of the zeroth mode.
"""
@inline function mode0dissipation(s, v, p, g)
  @. v.Uh = g.KKrsq^(p.nnu0-1) * abs2(s.solr)
  p.nu0/(g.Lx*g.Ly)*parsevalsum(v.Uh, g)
end
@inline mode0dissipation(prob) = mode0dissipation(prob.state, prob.vars, prob.params, prob.grid)

"""
    mode0drag(prob)

Returns the extraction of domain-averaged energy extraction by the drag μ.
"""
@inline function mode0drag(s, v, p, g)
  @. v.Uh = g.KKrsq^(p.nmu0-1) * abs2(s.solr)
  @. v.Uh[1, 1] = 0
  p.mu0/(g.Lx*g.Ly)*parsevalsum(v.Uh, g)
end
@inline mode0drag(prob) = mode0drag(prob.state, prob.vars, prob.params, prob.grid)

"""
    mode1dissipation(prob)

Returns the domain-averaged kinetic energy dissipation of the first mode
by horizontal viscosity.
"""
@inline function mode1dissipation(s, v, p, g)
  @views @. v.Uuh = g.k^p.nnu1*s.solc[:, :, 1]
  @views @. v.Vuh = g.l^p.nnu1*s.solc[:, :, 2]
  2*p.nu1/(g.Lx*g.Ly)*(parsevalsum2(v.Uuh, g) + parsevalsum2(v.Vuh, g))
end
@inline mode1dissipation(prob) = mode1dissipation(prob.state, prob.vars, prob.params, prob.grid)

"""
    mode1drag(prob)

Returns the domain-averaged kinetic energy dissipation of the first mode
by horizontal viscosity.
"""
@inline function mode1drag(s, v, p, g)
  @views @. v.Uuh = g.k^p.nmu1*s.solc[:, :, 1]
  @views @. v.Vuh = g.l^p.nmu1*s.solc[:, :, 2]
  if p.nmu1 != 0 # zero out zeroth mode
    @views @. v.Uuh[1, :] = 0
    @views @. v.Vuh[:, 1] = 0
  end
  2*p.mu1/(g.Lx*g.Ly)*(parsevalsum2(v.Uuh, g) + parsevalsum2(v.Vuh, g))
end
@inline mode1drag(prob) = mode1drag(prob.state, prob.vars, prob.params,
                                    prob.grid)

"""
    totalenergy(prob)

Returns the total energy projected onto the zeroth mode.
"""
@inline totalenergy(s, v, p, g) = mode0energy(s, v, g) + mode1energy(s, p, g)
@inline totalenergy(prob) = totalenergy(prob.state, prob.vars, prob.params, prob.grid)

"""
    shearp(prob)

Returns the domain-integrated shear production.
"""
function shearp(s, v, p, g)
  v.Zh .= s.solr
  @. v.Psih = -g.invKKrsq*v.Zh
  @. v.Uh  = -im*g.l  * v.Psih
  @. v.Vh  =  im*g.kr * v.Psih
  @. v.Uxh =  im*g.kr * v.Uh
  @. v.Uyh =  im*g.l  * v.Uh
  @. v.Vxh =  im*g.kr * v.Vh
  @. v.Vyh =  im*g.l  * v.Vh

  ldiv!(v.Ux, g.rfftplan, v.Uxh)
  ldiv!(v.Uy, g.rfftplan, v.Uyh)
  ldiv!(v.Vx, g.rfftplan, v.Vxh)
  ldiv!(v.Vy, g.rfftplan, v.Vyh)

  @views ldiv!(v.u, g.fftplan, s.solc[:, :, 1])
  @views ldiv!(v.v, g.fftplan, s.solc[:, :, 2])

  @. v.UZuzvw = real( 2.0*abs2(v.u)*v.Ux + 2.0*abs2(v.v)*v.Vy
                       + (conj(v.u)*v.v + v.u*conj(v.v))*(v.Uy + v.Vx) )

  g.dx*g.dy*sum(v.UZuzvw)
end
shearp(prob) = shearp(prob.state, prob.vars, prob.params, prob.grid)


"""
    conversion(prob)

Return the domain-integrated conversion from potential to kinetic energy.
"""
function conversion(s, v, p, g)
  @views @. v.wh = -(g.k*s.solc[:, :, 1] + g.l*s.solc[:, :, 2]) / p.m
  @views ldiv!(v.p, g.fftplan, s.solc[:, :, 3])
  ldiv!(v.w, g.fftplan, v.wh)
  # b = i*m*p
  @. v.UZuzvw = real(im*p.m*conj(v.w)*v.p - im*p.m*v.w*conj(v.p))
  g.dx*g.dy*sum(v.UZuzvw)
end
conversion(prob) = conversion(prob.state, prob.vars, prob.params, prob.grid)


"""
    mode0apv(prob)

Returns the projection of available potential vorticity onto the
zeroth mode.
"""
function mode0apv(Z, u, v, p, pr::TwoModeParams, g::TwoDGrid)
  Z .+ irfft( im*pr.m^2/pr.N^2 * (
      g.l .*rfft( @. real(u*conj(p) + conj(u)*p) )
    - g.kr.*rfft( @. real(v*conj(p) + conj(v)*p) )), g.nx)
end

function mode0apv(s, v::Vars, p::TwoModeParams, g::TwoDGrid)
  v.Z = irfft(s.solr, g.nx)
  @views ldiv!(v.u, g.fftplan, s.solc[:, :, 1])
  @views ldiv!(v.v, g.fftplan, s.solc[:, :, 2])
  @views ldiv!(v.p, g.fftplan, s.solc[:, :, 3])
  mode0apv(v.Z, v.u, v.v, v.p, p, g)
end
mode0apv(prob) = mode0apv(prob.state, prob.vars, prob.params, prob.grid)


"""
    mode1apv(prob)

Returns the projection of available potential energy onto the first mode.
"""
mode1apv(Z, zeta, p, pr, g) = @. zeta - pr.m^2/pr.N^2*(pr.f + Z)*p

function mode1apv(Z, s::DualState, v, p, g)
  @views @. v.ph = s.solc[:, :, 3]
  @views @. v.zetah = im*g.k*s.solc[:, :, 2] - im*g.l*s.solc[:, :, 1]

  ldiv!(v.p,  g.fftplan, v.ph)
  ldiv!(v.zeta,  g.fftplan, v.zetah)

  mode1apv(Z, v.zeta, v.p, p, g)
end

function mode1apv(s::DualState, v, p, g)
  v.Zh .= s.solr
  ldiv!(v.Z, g.rfftplan, v.Zh)
  mode1apv(v.Z, v, p, g)
end
mode1apv(prob) = mode1apv(prob.state, prob.vars, prob.params, prob.grid)


"""
    mode1u(prob)

Return the x-velocity associated with mode-1 at z=0.
"""
mode1u(v) = @. real(v.u + conj.(v.u))
mode1u(prob::AbstractProblem) = mode1u(prob.vars)

"""
    mode1v(prob)

Return the y-velocity associated with mode-1 at z=0.
"""
mode1v(v) = @. real(v.v + conj(v.v))
mode1v(prob::AbstractProblem) = mode1v(prob.vars)

"""
    mode1w(prob)

Return the z-velocity associated with mode-1 at z=0.
"""
mode1w(v) = @. real(v.w + conj(v.w))
mode1w(prob::AbstractProblem) = mode1w(prob.vars)

"""
    mode1p(prob)

Return the pressure associated with mode-1 at z=0.
"""
mode1p(v) = @. real(v.p + conj(v.p))
mode1p(prob::AbstractProblem) = mode1p(prob.vars)

"""
    mode1buoyancy(prob)
Return the buoyancy associated with mode-1 at z=0.
"""
mode1buoyancy(v, p) = @. real(im*p.m*v.p - im*p.m*conj(v.p))
mode1buoyancy(prob) = mode1buoyancy(prob.vars, prob.params)

"""
    mode1speed(prob)

Return the speed associated with mode-1 at z=0.
"""
mode1speed(v) = @. sqrt($mode1u(v)^2 + $mode1v(v)^2)
mode1speed(prob::AbstractProblem) = mode1speed(prob.vars)

"""
    mode0speed(prob)

Return the speed associated with mode-0.
"""
mode0speed(prob) = @. sqrt(prob.vars.U^2 + prob.vars.V^2)

end # module
