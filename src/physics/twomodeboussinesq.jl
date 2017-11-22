__precompile__()

# A module for solving a two-vertical-mode truncation of the Boussinesq 
# equations
module TwoModeBoussinesq

using FourierFlows




# Problem --------------------------------------------------------------------- 
""" 
Construct a TwoModeBoussinesq initial value problem.
"""
function InitialValueProblem(;
  nx   = 128, 
  Lx   = 2π, 
  ny   = nothing,
  Ly   = nothing,
  ν0   = nothing, 
  nν0  = 2, 
  ν1   = nothing,
  nν1  = 2, 
  f    = 1.0,
  N    = 10.0,
  m    = 40.0,
  Us   = 0.0,
  Vs   = 0.0,
  dt   = 0.01
  )

  if Ly == nothing; Ly = Lx; end
  if ny == nothing; ny = nx; end
  if ν0 == nothing; ν0 = 1e-1/(dt*(0.65π*nx/Lx)^nν0); end
  if ν1 == nothing; ν1 = 1e-1/(dt*(0.65π*nx/Lx)^nν1); end

  g  = TwoDGrid(nx, Lx, ny, Ly)
  pr = TwoModeBoussinesq.Params(ν0, nν0, ν1, nν1, f, N, m)
  vs = TwoModeBoussinesq.Vars(g)
  eq = TwoModeBoussinesq.Equation(pr, g)
  ts = ETDRK4TimeStepper(dt, eq.LCc, eq.LCr)

  FourierFlows.Problem(g, vs, pr, eq, ts)
end




# Params ---------------------------------------------------------------------- 
abstract type TwoModeParams <: AbstractParams end

type Params <: TwoModeParams
  ν0::Float64                 # Mode-0 viscosity
  nν0::Int                    # Mode-0 hyperviscous order
  ν1::Float64                 # Mode-1 viscosity
  nν1::Int                    # Mode-1 hyperviscous order
  f::Float64                  # Planetary vorticity
  N::Float64                  # Buoyancy frequency
  m::Float64                  # Mode-one wavenumber
  Us::Float64                 # Steady mode-0 mean x-velocity
  Vs::Float64                 # Steady mode-0 mean y-velocity
end

function Params(ν0, nν0::Int, ν1, nν1::Int, f, N, m)
  Params(ν0, nν0, ν1, nν1, f, N, m, 0.0, 0.0)
end




# Equations ------------------------------------------------------------------- 
type Equation <: AbstractEquation
  LCc::Array{Complex{Float64}, 3} # Element-wise coeff of the eqn's linear part
  LCr::Array{Complex{Float64}, 2} # Element-wise coeff of the eqn's linear part
  calcNL!::Function               # Function to calculate eqn's nonlinear part
end

function Equation(p::TwoModeParams, g::TwoDGrid)
  LCc, LCr = getlinearcoefficients(p, g)
  Equation(LCc, LCr, calcNL!)
end

function getlinearcoefficients(p::TwoModeParams, g::TwoDGrid)
  LCr = -p.ν0 * g.KKrsq.^(0.5*p.nν0)

  LCc = zeros(g.nk, g.nl, 3)
  LCc[:, :, 1] = -p.ν1 * g.KKsq.^(0.5*p.nν1)
  LCc[:, :, 2] = -p.ν1 * g.KKsq.^(0.5*p.nν1)

  LCc, LCr
end




# Vars ------------------------------------------------------------------------ 
abstract type TwoModeVars <: AbstractVars end

type Vars <: TwoModeVars
  t::Float64
  solr::Array{Complex128, 2}
  solc::Array{Complex128, 3}

  # Auxiliary zeroth-mode vars
  Z::Array{Float64, 2}
  U::Array{Float64, 2}
  V::Array{Float64, 2}
  UZuz::Array{Float64, 2}
  VZvz::Array{Float64, 2}
  Ux::Array{Float64, 2}
  Uy::Array{Float64, 2}
  Vx::Array{Float64, 2}
  Vy::Array{Float64, 2}
  Psi::Array{Float64, 2}

  # Auxiliary first-mode vars
  u::Array{Complex{Float64}, 2}
  v::Array{Complex{Float64}, 2}
  w::Array{Complex{Float64}, 2}
  p::Array{Complex{Float64}, 2}
  zeta::Array{Complex{Float64}, 2}

  # Multiplies
  Uu::Array{Complex{Float64}, 2}
  Uv::Array{Complex{Float64}, 2}
  Up::Array{Complex{Float64}, 2}
  Vu::Array{Complex{Float64}, 2}
  Vv::Array{Complex{Float64}, 2}
  Vp::Array{Complex{Float64}, 2}
  uUxvUy::Array{Complex{Float64}, 2}
  uVxvVy::Array{Complex{Float64}, 2}

  # Zeroth-mode transforms
  Zh::Array{Complex{Float64}, 2}
  Uh::Array{Complex{Float64}, 2}
  Vh::Array{Complex{Float64}, 2}
  UZuzh::Array{Complex{Float64}, 2}
  VZvzh::Array{Complex{Float64}, 2}
  Uxh::Array{Complex{Float64}, 2}
  Uyh::Array{Complex{Float64}, 2}
  Vxh::Array{Complex{Float64}, 2}
  Vyh::Array{Complex{Float64}, 2}
  Psih::Array{Complex{Float64}, 2}

  # First-mode transforms
  uh::Array{Complex{Float64}, 2}
  vh::Array{Complex{Float64}, 2}
  wh::Array{Complex{Float64}, 2}
  ph::Array{Complex{Float64}, 2}
  zetah::Array{Complex{Float64}, 2}

  # Multiply transforms
  Uuh::Array{Complex{Float64}, 2}
  Uvh::Array{Complex{Float64}, 2}
  Uph::Array{Complex{Float64}, 2}
  Vuh::Array{Complex{Float64}, 2}
  Vvh::Array{Complex{Float64}, 2}
  Vph::Array{Complex{Float64}, 2}
  uUxvUyh::Array{Complex{Float64}, 2}
  uVxvVyh::Array{Complex{Float64}, 2}
end

function Vars(g::TwoDGrid)
  # Initialize with t=0
  t = 0.0
  solc = zeros(Complex{Float64}, g.nk, g.nl, 3)
  solr = zeros(Complex{Float64}, g.nkr, g.nl)

  # Auxiliary zeroth-mode vars
  Z      = zeros(Float64, g.nx, g.ny)
  U      = zeros(Float64, g.nx, g.ny)
  V      = zeros(Float64, g.nx, g.ny)
  UZuz   = zeros(Float64, g.nx, g.ny)
  VZvz   = zeros(Float64, g.nx, g.ny)
  Ux     = zeros(Float64, g.nx, g.ny)
  Uy     = zeros(Float64, g.nx, g.ny)
  Vx     = zeros(Float64, g.nx, g.ny)
  Vy     = zeros(Float64, g.nx, g.ny)
  Psi    = zeros(Float64, g.nx, g.ny)
  
  # Auxiliary first-mode vars
  u      = zeros(Complex{Float64}, g.nx, g.ny)
  v      = zeros(Complex{Float64}, g.nx, g.ny)
  w      = zeros(Complex{Float64}, g.nx, g.ny)
  p      = zeros(Complex{Float64}, g.nx, g.ny)
  zeta   = zeros(Complex{Float64}, g.nx, g.ny)

  Uu     = zeros(Complex{Float64}, g.nx, g.ny)
  Uv     = zeros(Complex{Float64}, g.nx, g.ny)
  Up     = zeros(Complex{Float64}, g.nx, g.ny)
  Vu     = zeros(Complex{Float64}, g.nx, g.ny)
  Vv     = zeros(Complex{Float64}, g.nx, g.ny)
  Vp     = zeros(Complex{Float64}, g.nx, g.ny)

  uUxvUy = zeros(Complex{Float64}, g.nx, g.ny)
  uVxvVy = zeros(Complex{Float64}, g.nx, g.ny)

  # Transforms
  Zh     = zeros(Complex{Float64}, g.nkr, g.nl)
  Uh     = zeros(Complex{Float64}, g.nkr, g.nl)
  Vh     = zeros(Complex{Float64}, g.nkr, g.nl)
  UZuzh  = zeros(Complex{Float64}, g.nkr, g.nl)
  VZvzh  = zeros(Complex{Float64}, g.nkr, g.nl)
  Uxh    = zeros(Complex{Float64}, g.nkr, g.nl)
  Uyh    = zeros(Complex{Float64}, g.nkr, g.nl)
  Vxh    = zeros(Complex{Float64}, g.nkr, g.nl)
  Vyh    = zeros(Complex{Float64}, g.nkr, g.nl)
  Psih   = zeros(Complex{Float64}, g.nkr, g.nl)

  uh     = zeros(Complex{Float64}, g.nk, g.nl)
  vh     = zeros(Complex{Float64}, g.nk, g.nl)
  wh     = zeros(Complex{Float64}, g.nk, g.nl)
  ph     = zeros(Complex{Float64}, g.nk, g.nl)
  zetah  = zeros(Complex{Float64}, g.nk, g.nl)
  
  Uuh    = zeros(Complex{Float64}, g.nk, g.nl)
  Uvh    = zeros(Complex{Float64}, g.nk, g.nl)
  Uph    = zeros(Complex{Float64}, g.nk, g.nl)
  Vuh    = zeros(Complex{Float64}, g.nk, g.nl)
  Vvh    = zeros(Complex{Float64}, g.nk, g.nl)
  Vph    = zeros(Complex{Float64}, g.nk, g.nl)

  uUxvUyh= zeros(Complex{Float64}, g.nk, g.nl)
  uVxvVyh= zeros(Complex{Float64}, g.nk, g.nl)

  return Vars(t, solr, solc, 
    Z, U, V, UZuz, VZvz, Ux, Uy, Vx, Vy, Psi, 
    u, v, w, p, zeta, Uu, Uv, Up, Vu, Vv, Vp, uUxvUy, uVxvVy,
    Zh, Uh, Vh, UZuzh, VZvzh, Uxh, Uyh, Vxh, Vyh, Psih, 
    uh, vh, wh, ph, zetah, Uuh, Uvh, Uph, Vuh, Vvh, Vph, uUxvUyh, uVxvVyh,
    )
end




# Solvers --------------------------------------------------------------------- 
function calcNL!(
  NLc::Array{Complex{Float64}, 3},  NLr::Array{Complex{Float64}, 2}, 
  solc::Array{Complex{Float64}, 3}, solr::Array{Complex{Float64}, 2}, 
  t::Float64, v::Vars, p::TwoModeParams, g::TwoDGrid)
  
  # Spectral-space calculations
  @views @. v.wh = -1.0/p.m*(g.K*solc[:, :, 1] + g.L*solc[:, :, 2])

  # This copy is necessary because calling A_mul_B(v.Z, g.irfftplan, sol) 
  # a few lines below destroys sol when using Julia's FFTW.
  v.Zh .= solr

  @. v.Psih = -g.invKKrsq*v.Zh

  @. v.Uh = -im*g.Lr*v.Psih
  @. v.Vh =  im*g.Kr*v.Psih

  @. v.Uxh = im*g.Kr*v.Uh
  @. v.Vxh = im*g.Kr*v.Vh

  @. v.Uyh = im*g.Lr*v.Uh
  @. v.Vyh = im*g.Lr*v.Vh

  v.Uh[1, 1] += p.Us*g.nx*g.ny
  v.Vh[1, 1] += p.Vs*g.nx*g.ny

  @views @. v.zetah = im*g.K*solc[:, :, 2] - im*g.L*solc[:, :, 1]
 
  # Inverse transforms
  A_mul_B!(v.Z, g.irfftplan, v.Zh)
  A_mul_B!(v.U, g.irfftplan, v.Uh)
  A_mul_B!(v.V, g.irfftplan, v.Vh)

  A_mul_B!(v.Ux, g.irfftplan, v.Uxh)
  A_mul_B!(v.Uy, g.irfftplan, v.Uyh)
  A_mul_B!(v.Vx, g.irfftplan, v.Vxh)
  A_mul_B!(v.Vy, g.irfftplan, v.Vyh)

  @views A_mul_B!(v.u, g.ifftplan, solc[:, :, 1])
  @views A_mul_B!(v.v, g.ifftplan, solc[:, :, 2])
  @views A_mul_B!(v.p, g.ifftplan, solc[:, :, 3])

  A_mul_B!(v.zeta,  g.ifftplan, v.zetah)
  A_mul_B!(v.w,  g.ifftplan, v.wh)

  # Multiplies
  @. v.UZuz = (v.U * v.Z
    + real(v.u*conj(v.zeta) + conj(v.u)*v.zeta)
    + real(im*p.m*v.v*conj(v.w) - im*p.m*conj(v.v)*v.w)
  )

  @. v.VZvz = (v.V * v.Z
    + real(v.v*conj(v.zeta) + conj(v.v)*v.zeta)
    - real(im*p.m*v.u*conj(v.w) - im*p.m*conj(v.u)*v.w )
  )

  @. v.Uu = v.U * v.u
  @. v.Vu = v.V * v.u
  @. v.Uv = v.U * v.v
  @. v.Vv = v.V * v.v
  @. v.Up = v.U * v.p
  @. v.Vp = v.V * v.p

  @. v.uUxvUy = v.u*v.Ux + v.v*v.Uy
  @. v.uVxvVy = v.u*v.Vx + v.v*v.Vy

  # Forward transforms
  A_mul_B!(v.UZuzh, g.rfftplan, v.UZuz)
  A_mul_B!(v.VZvzh, g.rfftplan, v.VZvz)

  A_mul_B!(v.Uuh, g.fftplan, v.Uu)
  A_mul_B!(v.Uvh, g.fftplan, v.Uv)
  A_mul_B!(v.Vuh, g.fftplan, v.Vu)
  A_mul_B!(v.Vvh, g.fftplan, v.Vv)
  A_mul_B!(v.Uph, g.fftplan, v.Up)
  A_mul_B!(v.Vph, g.fftplan, v.Vp)

  A_mul_B!(v.uUxvUyh, g.fftplan, v.uUxvUy)
  A_mul_B!(v.uVxvVyh, g.fftplan, v.uVxvVy)


  # Zeroth-mode nonlinear term
  @. NLr = - im*g.Kr*v.UZuzh - im*g.Lr*v.VZvzh


  # First-mode nonlinear terms:
  # u
  @views @. NLc[:, :, 1] = ( p.f*solc[:, :, 2] - im*g.K*solc[:, :, 3]
    - im*g.K*v.Uuh - im*g.L*v.Vuh - v.uUxvUyh
  )

  # v
  @views @. NLc[:, :, 2] = ( -p.f*solc[:, :, 1] - im*g.L*solc[:, :, 3]
    - im*g.K*v.Uvh - im*g.L*v.Vvh - v.uVxvVyh
  )

  # p
  @views @. NLc[:, :, 3] = ( im*p.N^2.0/p.m*v.wh
    - im*g.K*v.Uph - im*g.L*v.Vph
  )

  dealias!(NLr, g)
  dealias!(NLc, g)

  nothing
end




# Helper functions ------------------------------------------------------------ 
function updatevars!(v::TwoModeVars, p::TwoModeParams, g::TwoDGrid, 
  Zh::AbstractArray)
  # We don't use A_mul_B here because irfft destroys its input.
  v.Z = irfft(Zh, g.nx)

  @. v.Psih = -g.invKKrsq*v.Zh
  @. v.Uh   = -im*g.Lr*v.Psih
  @. v.Vh   =  im*g.Kr*v.Psih
 
  # We don't use A_mul_B here because irfft destroys its input.
  v.Psi = irfft(v.Psih, g.nx)
  v.U = irfft(v.Uh, g.nx)
  v.V = irfft(v.Vh, g.nx)

  @views v.uh .= v.solc[:, :, 1]
  @views v.vh .= v.solc[:, :, 2]
  @views v.ph .= v.solc[:, :, 3]

  @. v.wh = -1.0/p.m*(g.K*v.uh + g.L*v.vh)

  A_mul_B!(v.u, g.ifftplan, v.uh)
  A_mul_B!(v.v, g.ifftplan, v.vh)
  A_mul_B!(v.p, g.ifftplan, v.ph)
  A_mul_B!(v.w, g.ifftplan, v.wh)

  nothing
end

function updatevars!(v::Vars, p::TwoModeParams, g::TwoDGrid)
  v.Zh .= v.solr
  updatevars!(v, p, g, v.Zh)
end

function updatevars!(prob::AbstractProblem)
  updatevars!(prob.vars, prob.params, prob.grid)
end




"""
Set zeroth mode vorticity and update vars. 
"""
function set_Z!(v::Vars, p::TwoModeParams, g::TwoDGrid, Z)
  A_mul_B!(v.solr, g.rfftplan, Z)
  updatevars!(v, p, g)
  nothing
end

function set_Z!(prob::AbstractProblem, Z)
  set_Z!(prob.vars, prob.params, prob.grid, Z)
end




""" 
Set first mode u, v, and p and update vars.
"""
function set_uvp!(vs::TwoModeVars, pr::TwoModeParams, g::TwoDGrid, u, v, p)
  @views A_mul_B!(vs.solc[:, :, 1], g.fftplan, u)
  @views A_mul_B!(vs.solc[:, :, 2], g.fftplan, v)
  @views A_mul_B!(vs.solc[:, :, 3], g.fftplan, p)
  updatevars!(vs, pr, g)
  nothing
end

function set_uvp!(prob::AbstractProblem, u, v, p)
  set_uvp!(prob.vars, prob.params, prob.grid, u, v, p)
end




""" 
Set a plane wave solution with initial speed uw and non-dimensional wave
number nkw. The dimensional wavenumber will be 2π*nkw/Lx. 
"""
function set_planewave!(vs::TwoModeVars, pr::TwoModeParams, g::TwoDGrid,
  uw::Real, nkw::Int)

  x, y = g.X, g.Y

  # Wave parameters
  kw = 2π*nkw/g.Lx
  σ = sqrt(pr.f^2 + pr.N^2*kw^2/pr.m^2)
  alpha = pr.N^2*kw^2/(pr.f^2*pr.m^2) # also (sig^2-f^2)/f^2

  u0 = uw/2
  v0 = -uw * im*pr.f/2σ
  p0 = uw * kw*pr.N^2/(2σ*pr.m^2)

  u = u0 * exp.(im*kw*x)
  v = v0 * exp.(im*kw*x)
  p = p0 * exp.(im*kw*x)
  
  set_uvp!(vs, pr, g, u, v, p)
  nothing
end

function set_planewave!(prob::AbstractProblem, uw::Real, nkw::Int)
  set_planewave!(prob.vars, prob.params, prob.grid, uw::Real, nkw::Int)
end




""" 
Generate an isotropic spectrum of waves.
"""
function set_isotropicwavefield!(
  vs::TwoModeVars, pr::TwoModeParams, g::TwoDGrid, amplitude::Function; KE=1.0,
  maxspeed=nothing)

  # For clarity
  f, N, m = pr.f, pr.N, pr.m

  # Initialize
  phase = zeros(Complex{Float64}, g.nx, g.ny)
  u0 = zeros(Complex{Float64}, g.nx, g.ny)
  v0 = zeros(Complex{Float64}, g.nx, g.ny)
  p0 = zeros(Complex{Float64}, g.nx, g.ny)

  # Sum Fourier components
  for k in real.(g.k), l in real.(g.l)
    if amplitude(k, l) > 1e-15
     
      # Dispersion relation
      σ = sqrt(f^2 + N^2/m^2*(k^2 + l^2))

      # Random phases
      phase .= k*g.X .+ l*g.Y .+ 2π*rand()

      # Polarization
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
   
  set_uvp!(vs, pr, g, u0, v0, p0)

  nothing
end

function set_isotropicwavefield!(
  vs::TwoModeVars, pr::TwoModeParams, g::TwoDGrid; kwargs...)
  amplitude(k, l) = 1.0
  set_isotropicwavefield!(vs, pr, g, amplitude; kwargs...)
end

function set_isotropicwavefield!(prob::AbstractProblem, amplitude::Function;
  kwargs...)
  set_isotropicwavefield!(prob.vars, prob.params, prob.grid, amplitude; 
    kwargs...)
end


 

""" 
Returns the integrated energy in the zeroth mode energy. 
"""
function mode0energy(v::Vars, p::TwoModeParams, g::TwoDGrid)
  0.5*FourierFlows.parsevalsum(real.(g.invKKrsq).*abs2.(v.solr), g)
end

function mode0energy(prob::AbstractProblem)
  mode0energy(prob.vars, prob.params, prob.grid)
end




""" 
Returns the projection of the integrated first mode kinetic energy
onto the zeroth mode.
"""
function mode1ke(uh, vh, g)
  FourierFlows.parsevalsum2(uh, g) + FourierFlows.parsevalsum2(vh, g)
end

function mode1ke(v::TwoModeVars, p::TwoModeParams, g::TwoDGrid)
  @views mode1ke(v.solc[:, :, 1], v.solc[:, :, 2], g)
end

function mode1ke(prob::AbstractProblem)
  mode1ke(prob.vars, prob.params, prob.grid)
end




""" 
Returns the projection of the integrated first mode potential energy onto the
zeroth mode. 
"""
function mode1pe(v::TwoModeVars, p::TwoModeParams, g::TwoDGrid)
  p.m^2/p.N^2*FourierFlows.parsevalsum2(v.solc[:, :, 3], g)
end

function mode1pe(prob::AbstractProblem)
  mode1pe(prob.vars, prob.params, prob.grid)
end




""" 
Returns the projection of the total integrated first mode energy onto the
zeroth mode.
"""
function mode1energy(v::TwoModeVars, p::TwoModeParams, g::TwoDGrid)
  mode1ke(v, p, g) + mode1pe(v, p, g)
end

function mode1energy(prob::AbstractProblem)
  mode1energy(prob.vars, prob.params, prob.grid)
end




""" 
Returns the total energy projected onto the zeroth mode.
"""
function totalenergy(v::TwoModeVars, p::TwoModeParams, g::TwoDGrid)
  mode0energy(v, p, g) + mode1energy(v, p, g)
end

function totalenergy(prob::AbstractProblem)
  totalenergy(prob.vars, prob.params, prob.grid)
end




""" 
Returns kinetic energy dissipation of the zeroth mode. 
"""
function mode0dissipation(v::TwoModeVars, p::TwoModeParams, g::TwoDGrid)
  delzeta = irfft(
    (-1.0)^(p.nν0/2) .* g.KKrsq.^(p.nν0/2) .* vs.solr, g.nx)
  -p.nu*g.dx*g.dy*sum(vs.Psi.*delzeta)
end




""" 
Returns the projection of available potential vorticity onto the 
zeroth mode.
"""
function mode0apv(Z, u, v, p, pr::TwoModeParams, g::TwoDGrid)
  (Z .+ irfft( pr.m^2.0./pr.N^2.0 .* (
      im.*g.Lr.*rfft(real.(u.*conj.(p) .+ conj.(u).*p)) 
    - im.*g.Kr.*rfft(real.(v.*conj.(p) .+ conj.(v).*p))
  ), g.nx))
end

function mode0apv(v::Vars, p::TwoModeParams, g::TwoDGrid)
  v.Z = irfft(v.solr, g.nx)
  @views A_mul_B!(v.u, g.ifftplan, v.solc[:, :, 1])
  @views A_mul_B!(v.v, g.ifftplan, v.solc[:, :, 2])
  @views A_mul_B!(v.p, g.ifftplan, v.solc[:, :, 3])
  mode0apv(v.Z, v.u, v.v, v.p, p, g)
end

function mode0apv(prob::AbstractProblem)
  mode0apv(prob.vars, prob.params, prob.grid)
end




""" 
Returns the projection of available potential energy onto the first mode.
"""
function mode1apv(Z, zeta, p, pr::TwoModeParams, g::TwoDGrid)
  zeta .- pr.m.^2.0./pr.N.^2.0 .* (pr.f .+ Z) .* p
end

function mode1apv(Z, v::TwoModeVars, p::TwoModeParams, g::TwoDGrid)
  @views @. v.ph = v.solc[:, :, 3]
  @views @. v.zetah = im*g.K*v.solc[:, :, 2] - im*g.L*v.solc[:, :, 1]

  A_mul_B!(v.p,  g.ifftplan, v.ph)
  A_mul_B!(v.zeta,  g.ifftplan, v.zetah)

  mode1apv(Z, v.zeta, v.p, p, g)
end




"""
Return the apv associated with mode-1.
"""

function mode1apv(v::Vars, p::TwoModeParams, g::TwoDGrid)
  v.Z = irfft(v.solr, g.nx)
  mode1apv(v.Z, v, p, g)
end

mode1apv(prob::AbstractProblem) = mode1apv(prob.vars, prob.params, prob.grid) 




"""
Return the x-velocity associated with mode-1 at z=0.
"""
mode1u(v::AbstractVars) = real.(v.u .+ conj.(v.u))
mode1u(prob::AbstractProblem) = mode1u(prob.vars)




"""
Return the y-velocity associated with mode-1 at z=0.
"""
mode1v(v::AbstractVars) = real.(v.v .+ conj.(v.v))
mode1v(prob::AbstractProblem) = mode1v(prob.vars)




"""
Return the z-velocity associated with mode-1 at z=0.
"""
mode1w(v::AbstractVars) = real.(v.w .+ conj.(v.w))
mode1w(prob::AbstractProblem) = mode1w(prob.vars)




"""
Return the pressure associated with mode-1 at z=0.
"""
mode1p(v::AbstractVars) = real.(v.p .+ conj.(v.p))
mode1p(prob::AbstractProblem) = mode1p(prob.vars)




"""
Return the buoyancy associated with mode-1 at z=0.
"""
mode1buoyancy(v, p) = real.(im.*p.m.*v.p .- im.*p.m.*conj.(v.p))
mode1buoyancy(prob::AbstractProblem) = mode1buoyancy(prob.vars, prob.params)




"""
Return the speed associated with mode-1 at z=0.
"""
mode1speed(v::AbstractVars) = sqrt.(mode1u(v).^2.0 .+ mode1v(v).^2.0)
mode1speed(prob::AbstractProblem) = mode1speed(prob.vars)




"""
Return the speed associated with mode-0.
"""
mode0speed(prob) = sqrt.(prob.vars.U.^2.0 .+ prob.vars.V.^2.0)




""" 
Returns the Courant-Freidrichs-Lewy number. 
"""
function CFL(prob, dt)
  dx = minimum([prob.grid.dx, prob.grid.dy])
  U = maximum([
    maximum(abs.(prob.vars.U)), 
    maximum(abs.(prob.vars.V)),
    maximum(abs.(prob.vars.u)),
    maximum(abs.(prob.vars.v))])

  U*dt/dx
end

 


# End module
end
