__precompile__()

# A module for solving a two-vertical-mode truncation of the Boussinesq 
# equations
module TwoModeBoussinesq

using FourierFlows

export Params,
       Vars,
       Equation

export set_zeta!, updatevars!


# Problem --------------------------------------------------------------------- 
""" Construct a FourierFlows Problem. """
function InitialValueProblem(;
  nx   = 128, 
  Lx   = 2pi, 
  ny   = nothing,
  Ly   = nothing,
  nu0  = nothing, 
  nnu0 = 2, 
  nu1  = nothing,
  nnu1 = 2, 
  f    = 1.0,
  N    = 10.0,
  m    = 40.0,
  Us   = 0.0,
  Vs   = 0.0,
  dt   = 0.01
  )

  if Ly == nothing; Ly = Lx; end
  if ny == nothing; ny = nx; end
  if nu0 == nothing; nu0 = 1e-1/(dt*(0.65*pi*nx/Lx)^nnu0); end
  if nu1 == nothing; nu1 = 1e-1/(dt*(0.65*pi*nx/Lx)^nnu1); end

  g  = TwoDGrid(nx, Lx)
  pr = TwoModeBoussinesq.Params(nu0, nnu0, nu1, nnu1, f, N, m)
  vs = TwoModeBoussinesq.Vars(g)
  eq = TwoModeBoussinesq.Equation(pr, g)
  ts = ETDRK4TimeStepper(dt, eq.LCc, eq.LCr)

  FourierFlows.Problem(g, vs, pr, eq, ts)
end



# Params ---------------------------------------------------------------------- 
abstract type TwoModeParams <: AbstractParams end

type Params <: TwoModeParams
  nu0::Float64                    # Mode-0 viscosity
  nnu0::Int                       # Mode-0 hyperviscous order
  nu1::Float64                    # Mode-1 viscosity
  nnu1::Int                       # Mode-1 hyperviscous order
  f::Float64                      # Planetary vorticity
  N::Float64                      # Buoyancy frequency
  m::Float64                      # Mode-one wavenumber
  Us::Float64                     # Steady mode-0 mean x-velocity
  Vs::Float64                     # Steady mode-0 mean y-velocity
end

function Params(nu0, nnu0::Int, nu1, nnu1::Int, f, N, m)
  Params(nu0, nnu0, nu1, nnu1, f, N, m, 0.0, 0.0)
end




# Equations ------------------------------------------------------------------- 
type Equation <: AbstractEquation
  LCc::Array{Complex{Float64}, 3}  # Element-wise coeff of the eqn's linear part
  LCr::Array{Complex{Float64}, 2}  # Element-wise coeff of the eqn's linear part
  calcNL!::Function               # Function to calculate eqn's nonlinear part
end

function Equation(p::TwoModeParams, g::TwoDGrid)
  LCr = -p.nu0 * g.KKrsq.^(0.5*p.nnu0)

  LCc = zeros(g.nx, g.ny, 3)
  LCc[:, :, 1] = -p.nu1 * g.KKsq.^(0.5*p.nnu1)
  LCc[:, :, 2] = -p.nu1 * g.KKsq.^(0.5*p.nnu1)
  LCc[:, :, 3] = -p.nu1 * g.KKsq.^(0.5*p.nnu1)

  # Function calcNL! is defined below.
  Equation(LCc, LCr, calcNL!)
end




# Vars ------------------------------------------------------------------------ 
type Vars <: AbstractVars

  # Z : zeta_0
  # u0, v0, psi0 : U, V, psi
  # u1, v1, w1, p1 : u, v, w, p

  t::Float64
  solr::Array{Complex128, 2}
  solc::Array{Complex128, 3}

  # Auxiliary zeroth-mode vars
  Z::Array{Float64, 2}
  U::Array{Float64, 2}
  V::Array{Float64, 2}
  UZ::Array{Float64, 2}
  VZ::Array{Float64, 2}
  Ux::Array{Float64, 2}
  Uy::Array{Float64, 2}
  Vx::Array{Float64, 2}
  Vy::Array{Float64, 2}
  uw::Array{Float64, 2}
  vw::Array{Float64, 2}
  uzeta::Array{Float64, 2}
  vzeta::Array{Float64, 2}
  psi::Array{Float64, 2}

  # Auxiliary first-mode vars
  u::Array{Complex{Float64}, 2}
  v::Array{Complex{Float64}, 2}
  w::Array{Complex{Float64}, 2}
  p::Array{Complex{Float64}, 2}
  vx::Array{Complex{Float64}, 2}
  uy::Array{Complex{Float64}, 2}
  zeta::Array{Complex{Float64}, 2}

  # Multiplies
  Uu::Array{Complex{Float64}, 2}
  Uv::Array{Complex{Float64}, 2}
  Up::Array{Complex{Float64}, 2}
  Vu::Array{Complex{Float64}, 2}
  Vv::Array{Complex{Float64}, 2}
  Vp::Array{Complex{Float64}, 2}
  uUx::Array{Complex{Float64}, 2}
  uVx::Array{Complex{Float64}, 2}
  vUy::Array{Complex{Float64}, 2}
  vVy::Array{Complex{Float64}, 2}

  # Zeroth-mode transforms
  Zh::Array{Complex{Float64}, 2}
  Uh::Array{Complex{Float64}, 2}
  Vh::Array{Complex{Float64}, 2}
  UZh::Array{Complex{Float64}, 2}
  VZh::Array{Complex{Float64}, 2}
  Uxh::Array{Complex{Float64}, 2}
  Uyh::Array{Complex{Float64}, 2}
  Vxh::Array{Complex{Float64}, 2}
  Vyh::Array{Complex{Float64}, 2}
  uwh::Array{Complex{Float64}, 2}
  vwh::Array{Complex{Float64}, 2}
  uzetah::Array{Complex{Float64}, 2}
  vzetah::Array{Complex{Float64}, 2}
  psih::Array{Complex{Float64}, 2}

  # First-mode transforms
  uh::Array{Complex{Float64}, 2}
  vh::Array{Complex{Float64}, 2}
  wh::Array{Complex{Float64}, 2}
  ph::Array{Complex{Float64}, 2}
  vxh::Array{Complex{Float64}, 2}
  uyh::Array{Complex{Float64}, 2}

  # Multiply transforms
  Uuh::Array{Complex{Float64}, 2}
  Uvh::Array{Complex{Float64}, 2}
  Uph::Array{Complex{Float64}, 2}
  Vuh::Array{Complex{Float64}, 2}
  Vvh::Array{Complex{Float64}, 2}
  Vph::Array{Complex{Float64}, 2}
  uUxh::Array{Complex{Float64}, 2}
  uVxh::Array{Complex{Float64}, 2}
  vUyh::Array{Complex{Float64}, 2}
  vVyh::Array{Complex{Float64}, 2}


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
  UZ     = zeros(Float64, g.nx, g.ny)
  VZ     = zeros(Float64, g.nx, g.ny)
  Ux     = zeros(Float64, g.nx, g.ny)
  Uy     = zeros(Float64, g.nx, g.ny)
  Vx     = zeros(Float64, g.nx, g.ny)
  Vy     = zeros(Float64, g.nx, g.ny)
  uw     = zeros(Float64, g.nx, g.ny)
  vw     = zeros(Float64, g.nx, g.ny)
  uzeta  = zeros(Float64, g.nx, g.ny)
  vzeta  = zeros(Float64, g.nx, g.ny)
  psi    = zeros(Float64, g.nx, g.ny)
  
  # Auxiliary first-mode vars
  u      = zeros(Complex{Float64}, g.nx, g.ny)
  v      = zeros(Complex{Float64}, g.nx, g.ny)
  w      = zeros(Complex{Float64}, g.nx, g.ny)
  p      = zeros(Complex{Float64}, g.nx, g.ny)
  vx     = zeros(Complex{Float64}, g.nx, g.ny)
  uy     = zeros(Complex{Float64}, g.nx, g.ny)
  zeta   = zeros(Complex{Float64}, g.nx, g.ny)

  Uu     = zeros(Complex{Float64}, g.nx, g.ny)
  Uv     = zeros(Complex{Float64}, g.nx, g.ny)
  Up     = zeros(Complex{Float64}, g.nx, g.ny)
  Vu     = zeros(Complex{Float64}, g.nx, g.ny)
  Vv     = zeros(Complex{Float64}, g.nx, g.ny)
  Vp     = zeros(Complex{Float64}, g.nx, g.ny)

  uUx    = zeros(Complex{Float64}, g.nx, g.ny)
  uVx    = zeros(Complex{Float64}, g.nx, g.ny)
  vUy    = zeros(Complex{Float64}, g.nx, g.ny)
  vVy    = zeros(Complex{Float64}, g.nx, g.ny)

  # Transforms
  Zh     = zeros(Complex{Float64}, g.nkr, g.nl)
  Uh     = zeros(Complex{Float64}, g.nkr, g.nl)
  Vh     = zeros(Complex{Float64}, g.nkr, g.nl)
  UZh    = zeros(Complex{Float64}, g.nkr, g.nl)
  VZh    = zeros(Complex{Float64}, g.nkr, g.nl)
  Uxh    = zeros(Complex{Float64}, g.nkr, g.nl)
  Uyh    = zeros(Complex{Float64}, g.nkr, g.nl)
  Vxh    = zeros(Complex{Float64}, g.nkr, g.nl)
  Vyh    = zeros(Complex{Float64}, g.nkr, g.nl)
  uwh    = zeros(Complex{Float64}, g.nkr, g.ny)
  vwh    = zeros(Complex{Float64}, g.nkr, g.ny)
  uzetah = zeros(Complex{Float64}, g.nkr, g.ny)
  vzetah = zeros(Complex{Float64}, g.nkr, g.ny)
  psih   = zeros(Complex{Float64}, g.nkr, g.nl)

  uh     = zeros(Complex{Float64}, g.nk, g.nl)
  vh     = zeros(Complex{Float64}, g.nk, g.nl)
  wh     = zeros(Complex{Float64}, g.nk, g.nl)
  ph     = zeros(Complex{Float64}, g.nk, g.nl)
  vxh    = zeros(Complex{Float64}, g.nk, g.nl)
  uyh    = zeros(Complex{Float64}, g.nk, g.nl)
  
  Uuh    = zeros(Complex{Float64}, g.nk, g.nl)
  Uvh    = zeros(Complex{Float64}, g.nk, g.nl)
  Uph    = zeros(Complex{Float64}, g.nk, g.nl)
  Vuh    = zeros(Complex{Float64}, g.nk, g.nl)
  Vvh    = zeros(Complex{Float64}, g.nk, g.nl)
  Vph    = zeros(Complex{Float64}, g.nk, g.nl)

  uUxh   = zeros(Complex{Float64}, g.nk, g.nl)
  uVxh   = zeros(Complex{Float64}, g.nk, g.nl)
  vUyh   = zeros(Complex{Float64}, g.nk, g.nl)
  vVyh   = zeros(Complex{Float64}, g.nk, g.nl)

  return Vars(t, solr, solc, 
    Z, U, V, UZ, VZ, Ux, Uy, Vx, Vy, uw, vw, uzeta, vzeta, psi, 
    u, v, w, p, vx, uy, zeta, Uu, Uv, Up, Vu, Vv, Vp, uUx, uVx, vUy, vVy,
    Zh, Uh, Vh, UZh, VZh, Uxh, Uyh, Vxh, Vyh, uwh, vwh, uzetah, vzetah, psih, 
    uh, vh, wh, ph, vxh, uyh, Uuh, Uvh, Uph, Vuh, Vvh, Vph, uUxh, uVxh, vUyh, 
    vVyh,
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

  @. v.psih = -g.invKKrsq*v.Zh

  @. v.Uh = -im*g.Lr*v.psih
  @. v.Vh =  im*g.Kr*v.psih

  @. v.Uxh = im*g.Kr*v.Uh
  @. v.Vxh = im*g.Kr*v.Vh

  @. v.Uyh = im*g.Lr*v.Uh
  @. v.Vyh = im*g.Lr*v.Vh

  v.Uh[1, 1] += p.Us*g.nx*g.ny
  v.Vh[1, 1] += p.Vs*g.nx*g.ny

  @views @. v.vxh = im*g.K*solc[:, :, 2]
  @views @. v.uyh = im*g.L*solc[:, :, 1]
 
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

  A_mul_B!(v.vx, g.ifftplan, v.vxh)
  A_mul_B!(v.uy, g.ifftplan, v.uyh)
  A_mul_B!(v.w,  g.ifftplan, v.wh)


  # Multiplies
  @. v.UZ = v.U * v.Z
  @. v.VZ = v.V * v.Z

  @. v.uw = real(im*p.m*v.u*conj(v.w) - im*p.m*conj(v.u)*v.w )
  @. v.vw = real(im*p.m*v.v*conj(v.w) - im*p.m*conj(v.v)*v.w )

  @. v.zeta = v.vx - v.uy
  @. v.uzeta = real(v.u*conj(v.zeta) + conj(v.u)*v.zeta)
  @. v.vzeta = real(v.v*conj(v.zeta) + conj(v.v)*v.zeta)
               
  @. v.Uu = v.U * v.u
  @. v.Vu = v.V * v.u
  @. v.Uv = v.U * v.v
  @. v.Vv = v.V * v.v
  @. v.Up = v.U * v.p
  @. v.Vp = v.V * v.p

  @. v.uUx = v.u * v.Ux
  @. v.uVx = v.u * v.Vx
  @. v.vUy = v.v * v.Uy
  @. v.vVy = v.v * v.Vy


  # Forward transforms
  A_mul_B!(v.UZh, g.rfftplan, v.UZ)
  A_mul_B!(v.VZh, g.rfftplan, v.VZ)

  A_mul_B!(v.uwh, g.rfftplan, v.uw)
  A_mul_B!(v.vwh, g.rfftplan, v.vw)

  A_mul_B!(v.uzetah, g.rfftplan, v.uzeta)
  A_mul_B!(v.vzetah, g.rfftplan, v.vzeta)

  A_mul_B!(v.Uuh, g.fftplan, v.Uu)
  A_mul_B!(v.Uvh, g.fftplan, v.Uv)
  A_mul_B!(v.Vuh, g.fftplan, v.Vu)
  A_mul_B!(v.Vvh, g.fftplan, v.Vv)
  A_mul_B!(v.Uph, g.fftplan, v.Up)
  A_mul_B!(v.Vph, g.fftplan, v.Vp)

  A_mul_B!(v.uUxh, g.fftplan, v.uUx)
  A_mul_B!(v.uVxh, g.fftplan, v.uVx)
  A_mul_B!(v.vUyh, g.fftplan, v.vUy)
  A_mul_B!(v.vVyh, g.fftplan, v.vVy)


  # ---------------------------------------------------------------------------   
  # Zeroth-mode nonlinear term
  @. NLr = ( - im*g.Kr*v.UZh    - im*g.Lr*v.VZh
             - im*g.Kr*v.uzetah - im*g.Lr*v.vzetah
             + im*g.Lr*v.uwh - im*g.Kr*v.vwh
  )

  # First-mode nonlinear terms:
  # u
  @views @. NLc[:, :, 1] = ( p.f*solc[:, :, 2] - im*g.K*solc[:, :, 3]
    - im*g.K*v.Uuh - im*g.L*v.Vuh - v.uUxh - v.vUyh
  )

  # v
  @views @. NLc[:, :, 2] = ( -p.f*solc[:, :, 1] - im*g.L*solc[:, :, 3]
    - im*g.K*v.Uvh - im*g.L*v.Vvh - v.uVxh - v.vVyh
  )

  # p
  @views @. NLc[:, :, 3] = ( im*p.N^2.0/p.m*v.wh
    - im*g.K*v.Uph - im*g.L*v.Vph
  )

  nothing
end







# Helper functions ------------------------------------------------------------ 
function updatevars!(v::Vars, p::TwoModeParams, g::TwoDGrid)

  v.Zh .= v.solr

  # We don't use A_mul_B here because irfft destroys its input.
  # A_mul_B!(v.Z, g.irfftplan, v.Zh)
  v.Z = irfft(v.Zh, g.nx)

  @. v.psih =         -g.invKKrsq*v.Zh
  @. v.Uh   =  im*g.Lr*g.invKKrsq*v.Zh
  @. v.Vh   = -im*g.Kr*g.invKKrsq*v.Zh
 
  # We don't use A_mul_B here because irfft destroys its input.
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

function updatevars!(prob::AbstractProblem)
  updatevars!(prob.vars, prob.params, prob.grid)
end




""" Set zeroth mode vorticity, zeta, and update vars. """
function set_zeta!(v::Vars, p::TwoModeParams, g::TwoDGrid, Z)
  A_mul_B!(v.solr, g.rfftplan, Z)
  updatevars!(v, p, g)
  nothing
end

function set_zeta!(prob::AbstractProblem, Z)
  set_zeta!(prob.vars, prob.params, prob.grid, Z)
end




""" Set first mode u, v, and p and update vars."""
function set_uvp!(vs::Vars, pr::TwoModeParams, g::TwoDGrid, u, v, p)
  uh = fft(u)
  vh = fft(v)
  ph = fft(p)

  vs.solc[:, :, 1] .= uh
  vs.solc[:, :, 2] .= vh
  vs.solc[:, :, 3] .= ph

  updatevars!(vs, pr, g)
  nothing
end

function set_uvp!(prob::AbstractProblem, u, v, p)
  set_uvp!(prob.vars, prob.params, prob.grid, u, v, p)
end




""" Set a plane wave solution with initial speed uw and non-dimensional wave
number nkw. The dimensional wavenumber will be 2π*nkw/Lx. """
function set_planewave!(vs::Vars, pr::TwoModeParams, g::TwoDGrid,
  uw::Real, nkw::Int)

  x, y = g.X, g.Y

  # Wave parameters
  kw = 2π*nkw/g.Lx
  sig = sqrt(pr.f^2 + pr.N^2*kw^2/pr.m^2)
  alpha = pr.N^2*kw^2/(pr.f^2*pr.m^2) # also (sig^2-f^2)/f^2

  # Component amplitudes
  u0 = uw * sig/(sqrt(2)*sqrt(alpha+2)*pr.f)
  v0 = pr.f/sig * u0
  p0 = 2*u0*alpha*pr.f^2 / (sig*kw)

  # Initial conditions
  u = u0     * exp.(im*kw*x)    # u = 2*u0*cos(phi)
  v = -im*v0 * exp.(im*kw*x)    # v = 2*v0*sin(phi)
  p = p0     * exp.(im*kw*x)    # p = 2*p0*cos(phi)
  
  set_uvp!(vs, pr, g, u, v, p)
  nothing
end

function set_planewave!(prob::AbstractProblem, uw::Real, nkw::Int; rotate=0.0)
  set_planewave!(prob.vars, prob.params, prob.grid, uw::Real, nkw::Int;
    rotate=rotate)
end




""" Calculate the zeroth mode energy. """
function mode0energy(v::Vars, p::TwoModeParams, g::TwoDGrid)
  0.5*FourierFlows.parsevalsum(real.(g.invKKrsq).*abs2.(v.solr), g)
end

function mode0energy(prob::AbstractProblem)
  mode0energy(prob.vars, prob.params, prob.grid)
end




""" Calculate the projection of the first mode kinetic energy onto the
zeroth mode. """
function mode1ke(v::Vars, p::TwoModeParams, g::TwoDGrid)
  (FourierFlows.parsevalsum2(v.solc[:, :, 1], g) 
    + FourierFlows.parsevalsum2(v.solc[:, :, 2], g))
end

function mode1ke(prob::AbstractProblem)
  mode1ke(prob.vars, prob.params, prob.grid)
end




""" Calculate the projection of the first mode potential energy onto the
zeroth mode. """
function mode1pe(v::Vars, p::TwoModeParams, g::TwoDGrid)
  p.m^2/p.N^2*FourierFlows.parsevalsum2(v.solc[:, :, 3], g)
end

function mode1pe(prob::AbstractProblem)
  mode1pe(prob.vars, prob.params, prob.grid)
end




""" Calculate the first mode energy. """
function mode1energy(v::Vars, p::TwoModeParams, g::TwoDGrid)
  mode1ke(v, p, g) + mode1pe(v, p, g)
end

function mode1energy(prob::AbstractProblem)
  mode1energy(prob.vars, prob.params, prob.grid)
end




""" 
Calculate the total energy projected onto the zeroth mode. 
"""
function totalenergy(v::Vars, p::TwoModeParams, g::TwoDGrid)
  mode0energy(v, p, g) + mode1energy(v, p, g)
end

function totalenergy(prob::AbstractProblem)
  totalenergy(prob.vars, prob.params, prob.grid)
end




""" 
Return the zeroth and first mode energy as a tuple. 
"""
function twoenergies(v::Vars, p::TwoModeParams, g::TwoDGrid)
  mode0energy(v, p, g), mode1energy(v, p, g)
end

function twoenergies(prob::AbstractProblem)
  twoenergies(prob.vars, prob.params, prob.grid)
end




""" Return kinetic energy dissipation of the zeroth mode. """
function mode0dissipation(v::Vars, p::TwoModeParams, g::TwoDGrid)
  delzeta = irfft(
    (-1.0)^(p.nnu0/2) .* g.KKrsq.^(p.nnu0/2) .* vs.solr, g.nx)
  -p.nu*g.dx*g.dy*sum(vs.psi.*delzeta)
end




""" Calculate the projection of APV onto the zeroth mode. """
function calc_apv(v::Vars, p::TwoModeParams, g::TwoDGrid)
  v.Z = irfft(v.solr, g.nx)
  @views A_mul_B!(v.u, g.ifftplan, v.solc[:, :, 1])
  @views A_mul_B!(v.v, g.ifftplan, v.solc[:, :, 2])
  @views A_mul_B!(v.p, g.ifftplan, v.solc[:, :, 3])

  (v.Z .+ irfft( p.m^2.0./p.N^2.0 .* (
      im.*g.Lr.*rfft(real.(v.u.*conj.(v.p) .+ conj.(v.u).*v.p)) 
    - im.*g.Kr.*rfft(real.(v.v.*conj.(v.p) .+ conj.(v.v).*v.p))
  ), g.nx))
end

function apv(prob::AbstractProblem)
  calc_apv(prob.vars, prob.params, prob.grid)
end


mode0speed(prob) = sqrt.(prob.vars.U.^2.0 + prob.vars.V.^2.0)

mode1u(v::AbstractVars) = real.(v.u + conj.(v.u))
mode1v(v::AbstractVars) = real.(v.v + conj.(v.v))
mode1w(v::AbstractVars) = real.(v.w + conj.(v.w))
mode1p(v::AbstractVars) = real.(v.p + conj.(v.p))

mode1u(prob::AbstractProblem) = mode1u(prob.vars)
mode1v(prob::AbstractProblem) = mode1v(prob.vars)
mode1w(prob::AbstractProblem) = mode1w(prob.vars)
mode1p(prob::AbstractProblem) = mode1p(prob.vars)

mode1speed(v::AbstractVars) = sqrt.(mode1u(v).^2.0 + mode1v(v).^2.0)
mode1speed(prob::AbstractProblem) = mode1speed(prob.vars)

mode1buoyancy(v, p) = real.(im*p.m*v.p - im*p.m*conj.(v.p))
mode1buoyancy(prob::AbstractProblem) = mode1buoyancy(prob.vars, prob.params)


""" Compute the transform of the Jacobian of two fields a, b on a grid g. """
function jacobianh(a, b, g::TwoDGrid)
  # J(a, b) = dx(a b_y) - dy(a b_x)
  bh = fft(b)
  bx = ifft(im*g.K.*bh)
  by = ifft(im*g.L.*bh)
  im*g.K.*fft(a.*bx) - im*g.L.*fft(a.*by)
end

""" Compute the Jacobian of two fields a, b on a grid g. """
function jacobian(a, b, g::TwoDGrid)
  ifft(jacobianh(a, b, g))
end









# Wave-induced flow and potential vorticity ----------------------------------- 
""" Calculate the wave-induced streamfunction and velocity fields. """
function wave_induced_uv(qw, g::TwoDGrid)

  qwh = rfft(qw)

  psiwh = g.invKKrsq.*qwh
  uwh   = -im*g.Lr.*psiwh
  vwh   =  im*g.Kr.*psiwh

  psiw = irfft(psiwh, g.nx)
  uw   = irfft(uwh, g.nx)
  vw   = irfft(vwh, g.nx)

  return uw, vw
end

""" Calculate the wave-induced streamfunction and velocity fields. """
function wave_induced_uv(sig::Real, v::Vars, p::TwoModeParams, g::TwoDGrid)
  wave_induced_uv(calc_qw(sig::Real, v::Vars, p::TwoModeParams, g::TwoDGrid), 
    g::TwoDGrid)
end

function wave_induced_uv(sig, prob::AbstractProblem)
  wave_induced_uv(sig, prob.vars, prob.params, prob.grid)
end




""" Calculate the wave contribution to PV, qw. """
function calc_qw(sig::Real, v::Vars, p::TwoModeParams, g::TwoDGrid)

  usig, vsig = calc_usigvsig(sig, v, p, g)

  # Non-Jacobian terms
  usig2xx = ifft(-g.K.^2.0.*fft(abs2.(usig)))
  vsig2yy = ifft(-g.L.^2.0.*fft(abs2.(vsig)))
  usigvsigxy = ifft(-g.K.*g.L.*fft(usig.*conj.(vsig) + conj.(usig).*vsig))

  # Assemble contributions
  qw = -real.( 
       2.0*im/sig * (
        jacobian(conj.(usig), usig, g) + jacobian(conj.(vsig), vsig, g))
    + p.f/sig^2.0 * (
        jacobian(conj.(vsig), usig, g) + jacobian(vsig, conj.(usig), g))
    + p.f/sig^2.0 * (usig2xx + vsig2yy + usigvsigxy)
  )

  return qw
end

""" Calculate the wave contribution to PV, qw. """
function calc_qw(usig::AbstractArray, vsig::AbstractArray, sig::Real, 
  p::TwoModeParams, g::TwoDGrid)

  usigh = fft(usig)
  vsigh = fft(vsig)

  # Non-Jacobian terms
  usig2xx = ifft(-g.K.^2.0.*fft(abs2.(usig)))
  vsig2yy = ifft(-g.L.^2.0.*fft(abs2.(vsig)))
  usigvsigxy = ifft(-g.K.*g.L.*fft(usig.*conj.(vsig) + conj.(usig).*vsig))

  # Assemble contributionsCalulate
  qw = -real.( 
       2.0*im/sig * (
        jacobian(conj.(usig), usig, g) + jacobian(conj.(vsig), vsig, g))
    + p.f/sig^2.0 * (
        jacobian(conj.(vsig), usig, g) + jacobian(vsig, conj.(usig), g))
    + p.f/sig^2.0 * (usig2xx + vsig2yy + usigvsigxy)
  )

  return qw
end

function calc_qw(sig::Real, prob::AbstractProblem)
  calc_qw(sig, prob.vars, prob.params, prob.grid)
end


""" Calculate usig and vsig, the complex, sigm-ified amplitudes 
of u_1 and v_1. """
function calc_usigvsig(sig, v::Vars, p::TwoModeParams, g::TwoDGrid)

  # Calculate u_t and v_t, and use them to find usig and vsig

  @views @. v.uh = v.solc[:, :, 1]
  @views @. v.vh = v.solc[:, :, 2]
  @views @. v.ph = v.solc[:, :, 3]

  # This copy is necessary because calling A_mul_B(v.Z, g.irfftplan, v.Zh) 
  # a few lines below destroys v.Zh
  v.Zh .= v.solr

  @. v.psih = -g.invKKrsq*v.Zh
  @. v.Uh   = -im*g.Lr*v.psih
  @. v.Vh   =  im*g.Kr*v.psih
  @. v.Uxh  = im*g.Kr*v.Uh
  @. v.Vxh  = im*g.Kr*v.Vh
  @. v.Uyh  = im*g.Lr*v.Uh
  @. v.Vyh  = im*g.Lr*v.Vh

  v.Uh[1, 1] += p.Us*g.nx*g.ny
  v.Vh[1, 1] += p.Vs*g.nx*g.ny

  # Inverse transforms
  A_mul_B!(v.U,  g.irfftplan, v.Uh)
  A_mul_B!(v.V,  g.irfftplan, v.Vh)
  A_mul_B!(v.Ux, g.irfftplan, v.Uxh)
  A_mul_B!(v.Uy, g.irfftplan, v.Uyh)
  A_mul_B!(v.Vx, g.irfftplan, v.Vxh)
  A_mul_B!(v.Vy, g.irfftplan, v.Vyh)

  A_mul_B!(v.u,  g.ifftplan, v.uh)
  A_mul_B!(v.v,  g.ifftplan, v.vh)

  @. v.Uu =  v.U * v.u
  @. v.Vu =  v.V * v.u
  @. v.Uv =  v.U * v.v
  @. v.Vv =  v.V * v.v
  @. v.uUx = v.u * v.Ux
  @. v.uVx = v.u * v.Vx
  @. v.vUy = v.v * v.Uy
  @. v.vVy = v.v * v.Vy

  A_mul_B!(v.Uuh,  g.fftplan, v.Uu)
  A_mul_B!(v.Uvh,  g.fftplan, v.Uv)
  A_mul_B!(v.Vuh,  g.fftplan, v.Vu)
  A_mul_B!(v.Vvh,  g.fftplan, v.Vv)
  A_mul_B!(v.uUxh, g.fftplan, v.uUx)
  A_mul_B!(v.uVxh, g.fftplan, v.uVx)
  A_mul_B!(v.vUyh, g.fftplan, v.vUy)
  A_mul_B!(v.vVyh, g.fftplan, v.vVy)

  # First-mode nonlinear terms:
  # u
  uth = ( -p.nu1*g.KKsq.^(0.5*p.nnu1) .* v.uh
    + p.f*v.vh - im*g.K.*v.ph
    - im*g.K.*v.Uuh - im*g.L.*v.Vuh - v.uUxh - v.vUyh
  )

  # v
  vth = ( -p.nu1*g.KKsq.^(0.5*p.nnu1) .* v.vh
    - p.f*v.uh - im*g.L.*v.ph
    - im*g.K.*v.Uvh - im*g.L.*v.Vvh - v.uVxh - v.vVyh
  )

  ut = ifft(uth)
  vt = ifft(vth)

  # Calculate amplitudes
  usig = exp(im*sig*v.t) * (v.u + im/sig*ut)
  vsig = exp(im*sig*v.t) * (v.v + im/sig*vt)

  return usig, vsig
end



""" Returns the wave-induced speed. """
function wave_induced_speed(vs, pr, g)
  uw, vw = wave_induced_uv(sig, vs, pr, g)
  return sqrt.(uw.^2.0 + vw.^2.0)
end

function wave_induced_speed(prob::AbstractProblem)
  wave_induced_speed(prob.vars, prob.params, prob.grid)
end


""" Returns the wave-induced x-velocity. """
function wave_induced_u(vs, pr, g)
  uw, vw = wave_induced_uv(sig, vs, pr, g)
  return uw
end

function wave_induced_u(prob::AbstractProblem)
  wave_induced_u(prob.vars, prob.params, prob.grid)
end


""" Returns the wave-induced y-velocity. """
function wave_induced_v(vs, pr, g)
  uw, vw = wave_induced_uv(sig, vs, pr, g)
  return vw
end

function wave_induced_v(prob::AbstractProblem)
  wave_induced_v(prob.vars, prob.params, prob.grid)
end


function apvinducedflow(vs, pr, g)
  q = TwoModeBoussinesq.calc_apv(vs, pr, g)
  psiqh = -g.invKKrsq.*rfft(q)
  uq = irfft(-im*g.Lr.*psiqh, g.nx)
  vq = irfft( im*g.Kr.*psiqh, g.nx)
  return sqrt.(uq.^2.0+vq.^2.0)
end



""" Calculate the Courant-Freidrich-Levant number. """
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
