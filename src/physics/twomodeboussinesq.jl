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

  g  = TwoDGrid(nx, Lx)
  pr = TwoModeBoussinesq.Params(ν0, nν0, ν1, nν1, f, N, m)
  vs = TwoModeBoussinesq.Vars(g)
  eq = TwoModeBoussinesq.Equation(pr, g)
  ts = ETDRK4TimeStepper(dt, eq.LCc, eq.LCr)

  FourierFlows.Problem(g, vs, pr, eq, ts)
end




"""
Returns a TwoModeBoussinesq problem in which the APV is used
instead of vorticity as the barotropic prognostic variable.
"""
function PrognosticAPVInitialValueProblem(;
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
  dt   = 0.01,
  nowavecubics = false,
  nowaveinducedflow = false,
  linearwaves = false
  )

  if sum([nowavecubics, nowaveinducedflow, linearwaves]) > 1
    throw("Only one type of equation may be seleted by nowavecubics, "*
            "nowaveinducedflow, or linearwaves.")
  end

  if Ly == nothing; Ly = Lx; end
  if ny == nothing; ny = nx; end
  if ν0 == nothing; ν0 = 1e-1/(dt*(0.65π*nx/Lx)^nν0); end
  if ν1 == nothing; ν1 = 1e-1/(dt*(0.65π*nx/Lx)^nν1); end

  g  = TwoDGrid(nx, Lx)
  pr = TwoModeBoussinesq.PrognosticAPVParams(ν0, nν0, ν1, nν1, f, N, m)
  vs = TwoModeBoussinesq.PrognosticAPVVars(g)

  if nowavecubics
    eq = TwoModeBoussinesq.NoWaveCubicsEquation(pr, g)
  elseif nowaveinducedflow
    eq = TwoModeBoussinesq.NoWaveInducedFlowEquation(pr, g)
  elseif linearwaves
    eq = TwoModeBoussinesq.LinearWavesEquation(pr, g)
  else
    eq = TwoModeBoussinesq.Equation(pr, g)
  end

  ts = ETDRK4TimeStepper(dt, eq.LCc, eq.LCr)

  FourierFlows.Problem(g, vs, pr, eq, ts)
end




"""
Returns a TwoModeBoussinesq problem in which APV is advected as a passive
scalar by an diagnosed Lagrangian-mean flow simultaneous to the evoution of
the wave field and barotropic vorticity field. The Lagrangian-mean flow can
be diagnosed when one wave frequency, σ, is present and known.
"""
function PassiveAPVInitialValueProblem(;
  nx  = 128, 
  Lx  = 2π, 
  ny  = nothing,
  Ly  = nothing,
  ν0  = nothing, 
  nν0 = 2, 
  ν1  = nothing,
  nν1 = 2, 
  f   = 1.0,
  N   = 10.0,
  m   = 40.0,
  Us  = 0.0,
  Vs  = 0.0,
  dt  = 0.01,
  σ   = f
  )

  if Ly == nothing; Ly = Lx; end
  if ny == nothing; ny = nx; end
  if ν0 == nothing; ν0 = 1e-1/(dt*(0.65π*nx/Lx)^nν0); end
  if ν1 == nothing; ν1 = 1e-1/(dt*(0.65π*nx/Lx)^nν1); end

  g  = TwoDGrid(nx, Lx)
  pr = TwoModeBoussinesq.PassiveAPVParams(
    ν0, nν0, ν1, nν1, f, N, m, σ)
  vs = TwoModeBoussinesq.PassiveAPVVars(g)
  eq = TwoModeBoussinesq.PassiveAPVEquation(pr, g)
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



"""
Defines a type for problems in which the baroclinic APV is zero, allowing
for the barotropic dynamics to be solved by evolving a barotropic APV 
equation.
"""
type PrognosticAPVParams <: TwoModeParams
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

function PrognosticAPVParams(ν0, nν0::Int, ν1, nν1::Int, f, N, m)
  PrognosticAPVParams(ν0, nν0, ν1, nν1, f, N, m, 0.0, 0.0)
end




type PassiveAPVParams <: TwoModeParams
  _params::Params
  ν0::Float64                 # Mode-0 viscosity
  nν0::Int                    # Mode-0 hyperviscous order
  ν1::Float64                 # Mode-1 viscosity
  nν1::Int                    # Mode-1 hyperviscous order
  f::Float64                  # Planetary vorticity
  N::Float64                  # Buoyancy frequency
  m::Float64                  # Mode-one wavenumber
  Us::Float64                 # Steady mode-0 mean x-velocity
  Vs::Float64                 # Steady mode-0 mean y-velocity
  σ::Float64                  # Wave frequency
end


function PassiveAPVParams(ν0, nν0::Int, ν1, nν1::Int, f, N, m, σ)
  _p = Params(ν0, nν0, ν1, nν1, f, N, m, 0.0, 0.0)
  PassiveAPVParams(_p, ν0, nν0, ν1, nν1, f, N, m, 0.0, 0.0, σ)
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

function NoWaveCubicsEquation(p::PrognosticAPVParams, g::TwoDGrid)
  LCc, LCr = getlinearcoefficients(p, g)
  Equation(LCc, LCr, calcNL_nowavecubics!)
end

function NoWaveInducedFlowEquation(p::PrognosticAPVParams, g::TwoDGrid)
  LCc, LCr = getlinearcoefficients(p, g)
  Equation(LCc, LCr, calcNL_nowaveinducedflow!)
end

function LinearWavesEquation!(p::PrognosticAPVParams, g::TwoDGrid)
  LCc, LCr = getlinearcoefficients(p, g)
  Equation(LCc, LCr, calcNL_linearwaves!)
end


function getlinearcoefficients(p::TwoModeParams, g::TwoDGrid)
  LCr = -p.ν0 * g.KKrsq.^(0.5*p.nν0)

  LCc = zeros(g.nk, g.nl, 3)
  LCc[:, :, 1] = -p.ν1 * g.KKsq.^(0.5*p.nν1)
  LCc[:, :, 2] = -p.ν1 * g.KKsq.^(0.5*p.nν1)

  LCc, LCr
end






type PassiveAPVEquation <: AbstractEquation
  LCc::Array{Complex{Float64}, 3} 
  LCr::Array{Complex{Float64}, 3}
  calcNL!::Function             
end

function PassiveAPVEquation(p::PassiveAPVParams, g::TwoDGrid)
  LCr = zeros(g.nkr, g.nl, 2)
  LCr[:, :, 1] = -p.ν0 * g.KKrsq.^(0.5*p.nν0)
  LCr[:, :, 2] = -p.ν0 * g.KKrsq.^(0.5*p.nν0)

  LCc = zeros(g.nk, g.nl, 3)
  LCc[:, :, 1] = -p.ν1 * g.KKsq.^(0.5*p.nν1)
  LCc[:, :, 2] = -p.ν1 * g.KKsq.^(0.5*p.nν1)

  PassiveAPVEquation(LCc, LCr, calcNL!)
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




"""
Defines a type for problems in which the barotropic APV is used as the 
prognostic variable and the baroclinic APV is assumed equal to zero.
"""
type PrognosticAPVVars <: TwoModeVars
  t::Float64
  solr::Array{Complex128, 2}
  solc::Array{Complex128, 3}

  # Auxiliary zeroth-mode vars
  Q::Array{Float64, 2}
  U::Array{Float64, 2}
  V::Array{Float64, 2}
  Ux::Array{Float64, 2}
  Uy::Array{Float64, 2}
  Vx::Array{Float64, 2}
  Vy::Array{Float64, 2}
  UQ::Array{Float64, 2}
  VQ::Array{Float64, 2}
  up::Array{Float64, 2}
  vp::Array{Float64, 2}
  Psi::Array{Float64, 2}

  # Auxiliary first-mode vars
  u::Array{Complex{Float64}, 2}
  v::Array{Complex{Float64}, 2}
  p::Array{Complex{Float64}, 2}

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
  Qh::Array{Complex{Float64}, 2}
  Uh::Array{Complex{Float64}, 2}
  Vh::Array{Complex{Float64}, 2}
  Uxh::Array{Complex{Float64}, 2}
  Uyh::Array{Complex{Float64}, 2}
  Vxh::Array{Complex{Float64}, 2}
  Vyh::Array{Complex{Float64}, 2}
  UQh::Array{Complex{Float64}, 2}
  VQh::Array{Complex{Float64}, 2}
  uph::Array{Complex{Float64}, 2}
  vph::Array{Complex{Float64}, 2}
  Psih::Array{Complex{Float64}, 2}

  # First-mode transforms
  uh::Array{Complex{Float64}, 2}
  vh::Array{Complex{Float64}, 2}
  ph::Array{Complex{Float64}, 2}

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

function PrognosticAPVVars(g::TwoDGrid)
  # Initialize with t=0
  t = 0.0
  solc = zeros(Complex{Float64}, g.nk, g.nl, 3)
  solr = zeros(Complex{Float64}, g.nkr, g.nl)

  # Auxiliary zeroth-mode vars
  Q      = zeros(Float64, g.nx, g.ny)
  U      = zeros(Float64, g.nx, g.ny)
  V      = zeros(Float64, g.nx, g.ny)
  Ux     = zeros(Float64, g.nx, g.ny)
  Uy     = zeros(Float64, g.nx, g.ny)
  Vx     = zeros(Float64, g.nx, g.ny)
  Vy     = zeros(Float64, g.nx, g.ny)
  UQ     = zeros(Float64, g.nx, g.ny)
  VQ     = zeros(Float64, g.nx, g.ny)
  up     = zeros(Float64, g.nx, g.ny)
  vp     = zeros(Float64, g.nx, g.ny)
  Psi    = zeros(Float64, g.nx, g.ny)
  
  # Auxiliary first-mode vars
  u      = zeros(Complex{Float64}, g.nx, g.ny)
  v      = zeros(Complex{Float64}, g.nx, g.ny)
  p      = zeros(Complex{Float64}, g.nx, g.ny)

  Uu     = zeros(Complex{Float64}, g.nx, g.ny)
  Uv     = zeros(Complex{Float64}, g.nx, g.ny)
  Up     = zeros(Complex{Float64}, g.nx, g.ny)
  Vu     = zeros(Complex{Float64}, g.nx, g.ny)
  Vv     = zeros(Complex{Float64}, g.nx, g.ny)
  Vp     = zeros(Complex{Float64}, g.nx, g.ny)

  uUxvUy = zeros(Complex{Float64}, g.nx, g.ny)
  uVxvVy = zeros(Complex{Float64}, g.nx, g.ny)

  # Transforms
  Qh     = zeros(Complex{Float64}, g.nkr, g.nl)
  Uh     = zeros(Complex{Float64}, g.nkr, g.nl)
  Vh     = zeros(Complex{Float64}, g.nkr, g.nl)
  Uxh    = zeros(Complex{Float64}, g.nkr, g.nl)
  Uyh    = zeros(Complex{Float64}, g.nkr, g.nl)
  Vxh    = zeros(Complex{Float64}, g.nkr, g.nl)
  Vyh    = zeros(Complex{Float64}, g.nkr, g.nl)
  UQh    = zeros(Complex{Float64}, g.nkr, g.nl)
  VQh    = zeros(Complex{Float64}, g.nkr, g.nl)
  uph    = zeros(Complex{Float64}, g.nkr, g.nl)
  vph    = zeros(Complex{Float64}, g.nkr, g.nl)
  Psih   = zeros(Complex{Float64}, g.nkr, g.nl)

  uh     = zeros(Complex{Float64}, g.nk, g.nl)
  vh     = zeros(Complex{Float64}, g.nk, g.nl)
  ph     = zeros(Complex{Float64}, g.nk, g.nl)
  
  Uuh    = zeros(Complex{Float64}, g.nk, g.nl)
  Uvh    = zeros(Complex{Float64}, g.nk, g.nl)
  Uph    = zeros(Complex{Float64}, g.nk, g.nl)
  Vuh    = zeros(Complex{Float64}, g.nk, g.nl)
  Vvh    = zeros(Complex{Float64}, g.nk, g.nl)
  Vph    = zeros(Complex{Float64}, g.nk, g.nl)

  uUxvUyh= zeros(Complex{Float64}, g.nk, g.nl)
  uVxvVyh= zeros(Complex{Float64}, g.nk, g.nl)

  return PrognosticAPVVars(t, solr, solc, 
    Q, U, V, Ux, Uy, Vx, Vy, UQ, VQ, up, vp, Psi, 
    u, v, p, Uu, Uv, Up, Vu, Vv, Vp, uUxvUy, uVxvVy,
    Qh, Uh, Vh, Uxh, Uyh, Vxh, Vyh, UQh, VQh, uph, vph, Psih, 
    uh, vh, ph, Uuh, Uvh, Uph, Vuh, Vvh, Vph, uUxvUyh, uVxvVyh,
    )
end





type PassiveAPVVars <: TwoModeVars

  _vars::Vars

  t::Float64
  solr::Array{Complex128, 3}
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

  Q::Array{Float64, 2}
  UL::Array{Float64, 2}
  VL::Array{Float64, 2}
  ULQ::Array{Float64, 2}
  VLQ::Array{Float64, 2}
  Uσ2::Array{Float64, 2}
  Vσ2::Array{Float64, 2}
  UσVσ::Array{Float64, 2}
  Juu::Array{Float64, 2}
  Jvv::Array{Float64, 2}
  Juv::Array{Float64, 2}
  up::Array{Float64, 2}
  vp::Array{Float64, 2}

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
  uUx::Array{Complex{Float64}, 2}
  uVx::Array{Complex{Float64}, 2}
  vUy::Array{Complex{Float64}, 2}
  vVy::Array{Complex{Float64}, 2}

  Uσ::Array{Complex{Float64}, 2}
  Vσ::Array{Complex{Float64}, 2}
  Uσx::Array{Complex{Float64}, 2}
  Vσx::Array{Complex{Float64}, 2}
  Uσy::Array{Complex{Float64}, 2}
  Vσy::Array{Complex{Float64}, 2}

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

  Qh::Array{Complex{Float64}, 2}
  Qwh::Array{Complex{Float64}, 2}
  ULh::Array{Complex{Float64}, 2}
  VLh::Array{Complex{Float64}, 2}
  ULQh::Array{Complex{Float64}, 2}
  VLQh::Array{Complex{Float64}, 2}
  Uσ2h::Array{Complex{Float64}, 2}
  Vσ2h::Array{Complex{Float64}, 2}
  UσVσh::Array{Complex{Float64}, 2}
  Juuh::Array{Complex{Float64}, 2}
  Jvvh::Array{Complex{Float64}, 2}
  Juvh::Array{Complex{Float64}, 2}
  uph::Array{Complex{Float64}, 2}
  vph::Array{Complex{Float64}, 2}

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
  uUxh::Array{Complex{Float64}, 2}
  uVxh::Array{Complex{Float64}, 2}
  vUyh::Array{Complex{Float64}, 2}
  vVyh::Array{Complex{Float64}, 2}

  uth::Array{Complex{Float64}, 2}
  vth::Array{Complex{Float64}, 2}
  Uσh::Array{Complex{Float64}, 2}
  Vσh::Array{Complex{Float64}, 2}
  Uσxh::Array{Complex{Float64}, 2}
  Vσxh::Array{Complex{Float64}, 2}
  Uσyh::Array{Complex{Float64}, 2}
  Vσyh::Array{Complex{Float64}, 2}
end

function PassiveAPVVars(g::TwoDGrid)
  # Initialize with t=0
  t = 0.0
  solc = zeros(Complex{Float64}, g.nk, g.nl, 3)
  solr = zeros(Complex{Float64}, g.nkr, g.nl, 2)

  v = Vars(g)

  Q      = zeros(Float64, g.nx, g.ny)
  UL     = zeros(Float64, g.nx, g.ny)
  VL     = zeros(Float64, g.nx, g.ny)
  ULQ    = zeros(Float64, g.nx, g.ny)
  VLQ    = zeros(Float64, g.nx, g.ny)
  Uσ2    = zeros(Float64, g.nx, g.ny)
  Vσ2    = zeros(Float64, g.nx, g.ny)
  UσVσ   = zeros(Float64, g.nx, g.ny)
  Juu    = zeros(Float64, g.nx, g.ny)
  Jvv    = zeros(Float64, g.nx, g.ny)
  Juv    = zeros(Float64, g.nx, g.ny)
  up     = zeros(Float64, g.nx, g.ny)
  vp     = zeros(Float64, g.nx, g.ny)

  Uσ     = zeros(Complex{Float64}, g.nx, g.ny)
  Vσ     = zeros(Complex{Float64}, g.nx, g.ny)
  Uσx    = zeros(Complex{Float64}, g.nx, g.ny)
  Vσx    = zeros(Complex{Float64}, g.nx, g.ny)
  Uσy    = zeros(Complex{Float64}, g.nx, g.ny)
  Vσy    = zeros(Complex{Float64}, g.nx, g.ny)

  Qh     = zeros(Complex{Float64}, g.nkr, g.nl)
  Qwh    = zeros(Complex{Float64}, g.nkr, g.nl)
  ULh    = zeros(Complex{Float64}, g.nkr, g.nl)
  VLh    = zeros(Complex{Float64}, g.nkr, g.nl)
  ULQh   = zeros(Complex{Float64}, g.nkr, g.nl)
  VLQh   = zeros(Complex{Float64}, g.nkr, g.nl)

  Uσ2h   = zeros(Complex{Float64}, g.nkr, g.nl)
  Vσ2h   = zeros(Complex{Float64}, g.nkr, g.nl)
  UσVσh  = zeros(Complex{Float64}, g.nkr, g.nl)
  Juuh   = zeros(Complex{Float64}, g.nkr, g.nl)
  Jvvh   = zeros(Complex{Float64}, g.nkr, g.nl)
  Juvh   = zeros(Complex{Float64}, g.nkr, g.nl)
  uph    = zeros(Complex{Float64}, g.nkr, g.ny)
  vph    = zeros(Complex{Float64}, g.nkr, g.ny)

  uth    = zeros(Complex{Float64}, g.nk, g.ny)
  vth    = zeros(Complex{Float64}, g.nk, g.ny)
  Uσh    = zeros(Complex{Float64}, g.nk, g.ny)
  Vσh    = zeros(Complex{Float64}, g.nk, g.ny)
  Uσxh   = zeros(Complex{Float64}, g.nk, g.ny)
  Vσxh   = zeros(Complex{Float64}, g.nk, g.ny)
  Uσyh   = zeros(Complex{Float64}, g.nk, g.ny)
  Vσyh   = zeros(Complex{Float64}, g.nk, g.ny)


  PassiveAPVVars(v, t, solr, solc, 
    v.Z, v.U, v.V, v.UZuz, v.VZvz, v.Ux, v.Uy, v.Vx, v.Vy, v.Psi,
    Q, UL, VL, ULQ, VLQ, Uσ2, Vσ2, UσVσ, Juu, Jvv, Juv, up, vp,
    v.u, v.v, v.w, v.p, v.zeta, v.Uu, v.Uv, v.Up, v.Vu, v.Vv, v.Vp, 
    v.uUxvUy, v.uVxvVy,
    Uσ, Vσ, Uσx, Vσx, Uσy, Vσy, 
    v.Zh, v.Uh, v.Vh, v.UZuzh, v.VZvzh, v.Uxh, v.Uyh, v.Vxh, v.Vyh, v.Psih,
    Qh, Qwh, ULh, VLh, ULQh, VLQh, Uσ2h, Vσ2h, UσVσh, Juuh, Jvvh, Juvh, uph, 
    vph,
    uh, v.vh, v.wh, v.ph, v.zetah, v.Uuh, v.Uvh, v.Uph, v.Vuh, v.Vvh, 
    v.Vph, v.uUxvUyh, v.uVxvVyh,
    uth, vth, Uσh, Vσh, Uσxh, Vσxh, Uσyh, Vσyh
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


  # ---------------------------------------------------------------------------   
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




""" 
Calculate the streamfunction of the barotropic flow.
"""
function calcPsih!(
  v::PrognosticAPVVars, p::PrognosticAPVParams, g::TwoDGrid,
  solc::Array{Complex{Float64}, 3}, solr::Array{Complex{Float64}, 2})
  
  @views A_mul_B!(v.u, g.ifftplan, solc[:, :, 1])
  @views A_mul_B!(v.v, g.ifftplan, solc[:, :, 2])
  @views A_mul_B!(v.p, g.ifftplan, solc[:, :, 3])

  @. v.up = real(v.u*conj(v.p) + conj(v.u)*v.p)
  @. v.vp = real(v.v*conj(v.p) + conj(v.v)*v.p)

  A_mul_B!(v.uph, g.rfftplan, v.up)
  A_mul_B!(v.vph, g.rfftplan, v.vp)

  @. v.Psih = -g.invKKrsq*(solr 
    + im*p.m^2.0/p.N^2.0*(g.Kr*v.vph - g.Lr*v.uph))

  nothing
end

function calcPsih!(v::PrognosticAPVVars, p::PrognosticAPVParams, g::TwoDGrid)
  calcPsih!(v, p, g, v.solc, v.solr)
end
  




"""
Calculate the nonlinear term for a problem where the prognostic barotropic
quantity is available potential vorticity under the assumption that the
baroclinic available potential vorticity is zero.
"""
function calcNL!(
  NLc::Array{Complex{Float64}, 3},  NLr::Array{Complex{Float64}, 2}, 
  solc::Array{Complex{Float64}, 3}, solr::Array{Complex{Float64}, 2}, 
  t::Float64, v::PrognosticAPVVars, p::PrognosticAPVParams, g::TwoDGrid)

  calcPsih!(v, p, g, solc, solr)

  @. v.Uh = -im*g.Lr*v.Psih
  @. v.Vh =  im*g.Kr*v.Psih

  @. v.Uxh = im*g.Kr*v.Uh
  @. v.Vxh = im*g.Kr*v.Vh

  @. v.Uyh = im*g.Lr*v.Uh
  @. v.Vyh = im*g.Lr*v.Vh

  v.Uh[1, 1] += p.Us*g.nx*g.ny
  v.Vh[1, 1] += p.Vs*g.nx*g.ny

  # Inverse transforms
  v.Qh .= solr
  A_mul_B!(v.Q, g.irfftplan, v.Qh)
  A_mul_B!(v.U, g.irfftplan, v.Uh)
  A_mul_B!(v.V, g.irfftplan, v.Vh)

  A_mul_B!(v.Ux, g.irfftplan, v.Uxh)
  A_mul_B!(v.Uy, g.irfftplan, v.Uyh)
  A_mul_B!(v.Vx, g.irfftplan, v.Vxh)
  A_mul_B!(v.Vy, g.irfftplan, v.Vyh)

  @. v.UQ = v.U * v.Q
  @. v.VQ = v.V * v.Q
  @. v.Uu = v.U * v.u
  @. v.Vu = v.V * v.u
  @. v.Uv = v.U * v.v
  @. v.Vv = v.V * v.v
  @. v.Up = v.U * v.p
  @. v.Vp = v.V * v.p

  @. v.uUxvUy = v.u*v.Ux + v.v*v.Uy
  @. v.uVxvVy = v.u*v.Vx + v.v*v.Vy

  # Forward transforms
  A_mul_B!(v.UQh, g.rfftplan, v.UQ)
  A_mul_B!(v.VQh, g.rfftplan, v.VQ)

  A_mul_B!(v.Uuh, g.fftplan, v.Uu)
  A_mul_B!(v.Uvh, g.fftplan, v.Uv)
  A_mul_B!(v.Vuh, g.fftplan, v.Vu)
  A_mul_B!(v.Vvh, g.fftplan, v.Vv)
  A_mul_B!(v.Uph, g.fftplan, v.Up)
  A_mul_B!(v.Vph, g.fftplan, v.Vp)

  A_mul_B!(v.uUxvUyh, g.fftplan, v.uUxvUy)
  A_mul_B!(v.uVxvVyh, g.fftplan, v.uVxvVy)


  # Zeroth-mode nonlinear term
  @. NLr = - im*g.Kr*v.UQh - im*g.Lr*v.VQh

  # First-mode nonlinear terms:
  # u
  @views @. NLc[:, :, 1] = ( p.f*solc[:, :, 2] 
    - im*(g.K*(solc[:, :, 3] + v.Uuh) + g.L*v.Vuh) - v.uUxvUyh
  )

  # v
  @views @. NLc[:, :, 2] = ( -p.f*solc[:, :, 1] 
      - im*(g.L*(solc[:, :, 3] + v.Vvh) + g.K*v.Uvh) - v.uVxvVyh
  )

  # p
  @views @. NLc[:, :, 3] = -im*( 
     g.K*(p.N^2.0/p.m^2.0*solc[:, :, 1] + v.Uph)
   + g.L*(p.N^2.0/p.m^2.0*solc[:, :, 2] + v.Vph)
  )

  dealias!(NLr, g)
  cubicdealias!(NLc, g)

  nothing
end

function calcNL_nowavecubics!(
  NLc::Array{Complex{Float64}, 3},  NLr::Array{Complex{Float64}, 2}, 
  solc::Array{Complex{Float64}, 3}, solr::Array{Complex{Float64}, 2}, 
  t::Float64, v::PrognosticAPVVars, p::PrognosticAPVParams, g::TwoDGrid)

  # APV advection
  calcPsih!(v, p, g, solc, solr)

  @. v.Uh = -im*g.Lr*v.Psih
  @. v.Vh =  im*g.Kr*v.Psih

  v.Uh[1, 1] += p.Us*g.nx*g.ny
  v.Vh[1, 1] += p.Vs*g.nx*g.ny

  v.Qh .= solr
  A_mul_B!(v.Q, g.irfftplan, v.Qh)
  A_mul_B!(v.U, g.irfftplan, v.Uh)
  A_mul_B!(v.V, g.irfftplan, v.Vh)

  @. v.UQ = v.U * v.Q
  @. v.VQ = v.V * v.Q

  A_mul_B!(v.UQh, g.rfftplan, v.UQ)
  A_mul_B!(v.VQh, g.rfftplan, v.VQ)

  @. NLr = - im*g.Kr*v.UQh - im*g.Lr*v.VQh


  # Wave terms, redefining Psi and neglecting waves
  @. v.Psih = -g.invKKrsq*solr

  A_mul_B!(v.U, g.irfftplan, v.Uh)
  A_mul_B!(v.V, g.irfftplan, v.Vh)

  @. v.Uh = -im*g.Lr*v.Psih
  @. v.Vh =  im*g.Kr*v.Psih

  v.Uh[1, 1] += p.Us*g.nx*g.ny
  v.Vh[1, 1] += p.Vs*g.nx*g.ny

  @. v.Uxh = im*g.Kr*v.Uh
  @. v.Vxh = im*g.Kr*v.Vh
  @. v.Uyh = im*g.Lr*v.Uh
  @. v.Vyh = im*g.Lr*v.Vh

  A_mul_B!(v.Ux, g.irfftplan, v.Uxh)
  A_mul_B!(v.Uy, g.irfftplan, v.Uyh)
  A_mul_B!(v.Vx, g.irfftplan, v.Vxh)
  A_mul_B!(v.Vy, g.irfftplan, v.Vyh)

  @. v.Uu = v.U * v.u
  @. v.Vu = v.V * v.u
  @. v.Uv = v.U * v.v
  @. v.Vv = v.V * v.v
  @. v.Up = v.U * v.p
  @. v.Vp = v.V * v.p

  @. v.uUxvUy = v.u*v.Ux + v.v*v.Uy
  @. v.uVxvVy = v.u*v.Vx + v.v*v.Vy

  # Forward transforms
  A_mul_B!(v.Uuh, g.fftplan, v.Uu)
  A_mul_B!(v.Uvh, g.fftplan, v.Uv)
  A_mul_B!(v.Vuh, g.fftplan, v.Vu)
  A_mul_B!(v.Vvh, g.fftplan, v.Vv)
  A_mul_B!(v.Uph, g.fftplan, v.Up)
  A_mul_B!(v.Vph, g.fftplan, v.Vp)

  A_mul_B!(v.uUxvUyh, g.fftplan, v.uUxvUy)
  A_mul_B!(v.uVxvVyh, g.fftplan, v.uVxvVy)


  # ---------------------------------------------------------------------------   
  # Zeroth-mode nonlinear term

  # First-mode nonlinear terms:
  # u
  @views @. NLc[:, :, 1] = ( p.f*solc[:, :, 2] 
    - im*(g.K*(solc[:, :, 3] + v.Uuh) + g.L*v.Vuh) - v.uUxvUyh
  )

  # v
  @views @. NLc[:, :, 2] = ( -p.f*solc[:, :, 1] 
      - im*(g.L*(solc[:, :, 3] + v.Vvh) + g.K*v.Uvh) - v.uVxvVyh
  )

  # p
  @views @. NLc[:, :, 3] = -im*( 
     g.K*(p.N^2.0/p.m^2.0*solc[:, :, 1] + v.Uph)
   + g.L*(p.N^2.0/p.m^2.0*solc[:, :, 2] + v.Vph)
  )

  dealias!(NLr, g)
  dealias!(NLc, g)

  nothing
end




function calcNL_nowaveinducedflow!(
  NLc::Array{Complex{Float64}, 3},  NLr::Array{Complex{Float64}, 2}, 
  solc::Array{Complex{Float64}, 3}, solr::Array{Complex{Float64}, 2}, 
  t::Float64, v::PrognosticAPVVars, p::PrognosticAPVParams, g::TwoDGrid)


  # APV-induced flow; no wave-induced flow
  @. v.Psih = -g.invKKrsq*solr

  v.Qh .= solr
  A_mul_B!(v.Q, g.irfftplan, v.Qh)
  A_mul_B!(v.U, g.irfftplan, v.Uh)
  A_mul_B!(v.V, g.irfftplan, v.Vh)

  @. v.Uh = -im*g.Lr*v.Psih
  @. v.Vh =  im*g.Kr*v.Psih

  v.Uh[1, 1] += p.Us*g.nx*g.ny
  v.Vh[1, 1] += p.Vs*g.nx*g.ny

  @. v.Uxh = im*g.Kr*v.Uh
  @. v.Vxh = im*g.Kr*v.Vh
  @. v.Uyh = im*g.Lr*v.Uh
  @. v.Vyh = im*g.Lr*v.Vh

  A_mul_B!(v.Ux, g.irfftplan, v.Uxh)
  A_mul_B!(v.Uy, g.irfftplan, v.Uyh)
  A_mul_B!(v.Vx, g.irfftplan, v.Vxh)
  A_mul_B!(v.Vy, g.irfftplan, v.Vyh)

  @. v.UQ = v.U * v.Q
  @. v.VQ = v.V * v.Q
  @. v.Uu = v.U * v.u
  @. v.Vu = v.V * v.u
  @. v.Uv = v.U * v.v
  @. v.Vv = v.V * v.v
  @. v.Up = v.U * v.p
  @. v.Vp = v.V * v.p

  @. v.uUxvUy = v.u*v.Ux + v.v*v.Uy
  @. v.uVxvVy = v.u*v.Vx + v.v*v.Vy

  # Forward transforms
  A_mul_B!(v.UQh, g.rfftplan, v.UQ)
  A_mul_B!(v.VQh, g.rfftplan, v.VQ)
  A_mul_B!(v.Uuh, g.fftplan, v.Uu)
  A_mul_B!(v.Uvh, g.fftplan, v.Uv)
  A_mul_B!(v.Vuh, g.fftplan, v.Vu)
  A_mul_B!(v.Vvh, g.fftplan, v.Vv)
  A_mul_B!(v.Uph, g.fftplan, v.Up)
  A_mul_B!(v.Vph, g.fftplan, v.Vp)

  A_mul_B!(v.uUxvUyh, g.fftplan, v.uUxvUy)
  A_mul_B!(v.uVxvVyh, g.fftplan, v.uVxvVy)


  # Zeroth-mode nonlinear term
  @. NLr = - im*g.Kr*v.UQh - im*g.Lr*v.VQh


  # First-mode nonlinear terms:
  # u
  @views @. NLc[:, :, 1] = ( p.f*solc[:, :, 2] 
    - im*(g.K*(solc[:, :, 3] + v.Uuh) + g.L*v.Vuh) - v.uUxvUyh
  )

  # v
  @views @. NLc[:, :, 2] = ( -p.f*solc[:, :, 1] 
      - im*(g.L*(solc[:, :, 3] + v.Vvh) + g.K*v.Uvh) - v.uVxvVyh
  )

  # p
  @views @. NLc[:, :, 3] = -im*( 
     g.K*(p.N^2.0/p.m^2.0*solc[:, :, 1] + v.Uph)
   + g.L*(p.N^2.0/p.m^2.0*solc[:, :, 2] + v.Vph)
  )

  dealias!(NLr, g)
  dealias!(NLc, g)

  nothing
end


function calcNL_linearwaves!(
  NLc::Array{Complex{Float64}, 3},  NLr::Array{Complex{Float64}, 2}, 
  solc::Array{Complex{Float64}, 3}, solr::Array{Complex{Float64}, 2}, 
  t::Float64, v::PrognosticAPVVars, p::PrognosticAPVParams, g::TwoDGrid)

  calcPsih!(v, p, g, solc, solr)

  @. v.Uh = -im*g.Lr*v.Psih
  @. v.Vh =  im*g.Kr*v.Psih

  v.Uh[1, 1] += p.Us*g.nx*g.ny
  v.Vh[1, 1] += p.Vs*g.nx*g.ny

  # Inverse transforms
  v.Qh .= solr
  A_mul_B!(v.Q, g.irfftplan, v.Qh)
  A_mul_B!(v.U, g.irfftplan, v.Uh)
  A_mul_B!(v.V, g.irfftplan, v.Vh)

  @. v.UQ = v.U * v.Q
  @. v.VQ = v.V * v.Q

  # Forward transforms
  A_mul_B!(v.UQh, g.rfftplan, v.UQ)
  A_mul_B!(v.VQh, g.rfftplan, v.VQ)
  
  # ---------------------------------------------------------------------------   
  # Zeroth-mode nonlinear term
  @. NLr = - im*g.Kr*v.UQh - im*g.Lr*v.VQh

  # First-mode nonlinear terms:
  # u
  @views @. NLc[:, :, 1] = p.f*solc[:, :, 2] - im*g.K*solc[:, :, 3]
  @views @. NLc[:, :, 2] = -p.f*solc[:, :, 1] - im*g.L*solc[:, :, 3]
  @views @. NLc[:, :, 3] = -im*( 
     p.N^2.0/p.m^2.0*(g.K*solc[:, :, 1]  + g.L*solc[:, :, 2]))

  dealias!(NLr, g)

  nothing
end




""" Calculate the nonlinear right side for the two-mode Boussinesq equations
and a passive tracer-advection equation, where the tracer is available 
potential vorticity. """
function calcNL!(
  NLc::Array{Complex{Float64}, 3},  NLr::Array{Complex{Float64}, 3}, 
  solc::Array{Complex{Float64}, 3}, solr::Array{Complex{Float64}, 3}, 
  t::Float64, v::PassiveAPVVars, p::PassiveAPVParams, g::TwoDGrid)
  
  @views calcNL!(NLc, NLr[:, :, 1], solc, solr[:, :, 1], t, v._vars, 
    p._params, g)

  calcULVL!(v, solr, colc, p, g) 

  @. v.ULQ = v.UL*v.Q
  @. v.ULQ = v.VL*v.Q

  A_mul_B!(v.ULQh, g.rfftplan, v.ULQ)
  A_mul_B!(v.VLQh, g.rfftplan, v.VLQ)

  @views @. NLr[:, :, 2] = -im*g.Kr*v.ULQh - im*g.Lr*v.VLQh

  dealias!(NLr, g)
  dealias!(NLc, g)

  nothing
end




function calcULVL!(v::PassiveAPVVars, solr, solc, p::PassiveAPVParams, 
  g::TwoDGrid)

  @views @. v.Qh = solr[:, :, 2]
  @views @. v.Uσh = exp(im*p.σ*v.t) * solc[:, :, 1]
  @views @. v.Vσh = exp(im*p.σ*v.t) * solc[:, :, 2]

  # Calculate Qwh
  @. v.Uσxh = im*g.K*v.Uσh
  @. v.Uσyh = im*g.L*v.Uσh
  @. v.Vσxh = im*g.K*v.Vσh
  @. v.Vσyh = im*g.L*v.Vσh

  A_mul_B!(v.Uσ,  g.ifftplan, v.Uσh)
  A_mul_B!(v.Uσx, g.ifftplan, v.Uσxh)
  A_mul_B!(v.Uσy, g.ifftplan, v.Uσyh)
  A_mul_B!(v.Vσ,  g.ifftplan, v.Vσh)
  A_mul_B!(v.Vσx, g.ifftplan, v.Vσxh)
  A_mul_B!(v.Vσy, g.ifftplan, v.Vσyh)

  @. v.Juu = real(im*conj(v.Uσx)*v.Uσy - im*conj(v.Uσy)*v.Uσx)
  @. v.Jvv = real(im*conj(v.Vσx)*v.Vσy - im*conj(v.Vσy)*v.Vσx)

  @. v.Juv = real(
    conj(v.Vσx)*v.Uσy - conj(v.Vσy)*v.Uσx + 
    v.Vσx*conj(v.Uσy) - v.Vσy*conj(v.Uσx)
  )

  # Non-Jacobian terms
  @. v.Uσ2 = abs2(v.Uσ)
  @. v.Vσ2 = abs2(v.Vσ)
  @. v.UσVσ = real(v.Uσ*conj(v.Vσ) + conj(v.Uσ)*v.Vσ)

  A_mul_B!(v.Juuh,  g.rfftplan, v.Juu)
  A_mul_B!(v.Jvvh,  g.rfftplan, v.Jvv)
  A_mul_B!(v.Juvh,  g.rfftplan, v.Juv)
  A_mul_B!(v.Uσ2h,  g.rfftplan, v.Uσ2)
  A_mul_B!(v.Vσ2h,  g.rfftplan, v.Vσ2)
  A_mul_B!(v.UσVσh, g.rfftplan, v.UσVσ)

  @. v.Qwh = -2.0/p.σ*(v.Juuh + v.Jvvh) - p.f/p.σ^2.0*(
       v.Juvh - g.Kr2*v.Uσ2h - g.Lr2*v.Vσ2h - g.KLr*v.UσVσh)
  
  @. v.ULh =  im*g.Lr*g.invKKrsq*(v.Qh - v.Qwh)
  @. v.VLh = -im*g.Kr*g.invKKrsq*(v.Qh - v.Qwh)

  v.ULh[1, 1] += p.Us*g.nx*g.ny
  v.VLh[1, 1] += p.Vs*g.nx*g.ny

  # Calculate nonlinear term.
  A_mul_B!(v.Q, g.irfftplan, v.ULh) 
  A_mul_B!(v.UL, g.irfftplan, v.ULh) 
  A_mul_B!(v.VL, g.irfftplan, v.VLh) 

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

function updatevars!(v::PassiveAPVVars, p::PassiveAPVParams, g::TwoDGrid)
  @views @. v.Zh = v.solr[:, :, 1]
  @views @. v.Qh = v.solr[:, :, 2]
  v.Q = irfft(v.Qh, g.nx)
  updatevars!(v, p, g, v.Zh)
  calcULVL!(v, v.solr, v.solc, p, g)
  nothing
end

function updatevars!(v::PrognosticAPVVars, p::TwoModeParams, g::TwoDGrid) 

  calcPsih!(v, p, g, v.solc, v.solr)
  @. v.Uh = -im*g.Lr*v.Psih
  @. v.Vh =  im*g.Kr*v.Psih
 
  # Update vars destroys Psih, Uh, Vh, Qh
  v.Qh .= v.solr
  A_mul_B!(v.Psi, g.irfftplan, v.Psih)
  A_mul_B!(v.U, g.irfftplan, v.Uh)
  A_mul_B!(v.V, g.irfftplan, v.Vh)
  A_mul_B!(v.Q, g.irfftplan, v.Qh)

  nothing
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

function set_Z!(v::PassiveAPVVars, p::PassiveAPVParams, g::TwoDGrid, Z)
  A_mul_B!(v.Zh, g.rfftplan, Z)
  @views @. v.solr[:, :, 1] = v.Zh

  Q = mode0apv(v, p, g)
  Qh = rfft(Q)

  @. v.Q = Q
  @views @. v.solr[:, :, 2] = Qh

  updatevars!(v, p, g)
  nothing
end




"""
Set the zeroth mode available potential vorticity.
"""
function set_Q!(v::PrognosticAPVVars, p::TwoModeParams, g::TwoDGrid, Q)
  A_mul_B!(v.solr, g.rfftplan, Q)
  updatevars!(v, p, g)
  nothing
end

function set_Q!(prob::AbstractProblem, Q)
  set_Q!(prob.vars, prob.params, prob.grid, Q)
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

function set_uvp!(vs::PassiveAPVVars, pr::PassiveAPVParams, g::TwoDGrid, 
  u, v, p)

  @views A_mul_B!(vs.solc[:, :, 1], g.fftplan, u)
  @views A_mul_B!(vs.solc[:, :, 2], g.fftplan, v)
  @views A_mul_B!(vs.solc[:, :, 3], g.fftplan, p)

  Q = mode0apv(vs, pr, g)
  Qh = rfft(Q)

  @. v.Q = Q
  @views @. vs.solr[:, :, 2] = Qh

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

function mode0energy(v::PassiveAPVVars, p::TwoModeParams, g::TwoDGrid)
  0.5*FourierFlows.parsevalsum(real.(g.invKKrsq).*abs2.(v.solr[:, :, 1]), g)
end

function mode0energy(v::PrognosticAPVVars, p::PrognosticAPVParams, g::TwoDGrid)
  calcPsih!(v, p, g, v.solc, v.solr)
  0.5*FourierFlows.parsevalsum(real.(g.KKrsq).*abs2.(v.Psih), g)
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

function mode0apv(v::PassiveAPVVars, p::TwoModeParams, g::TwoDGrid)
  v.Z = irfft(v.solr[:, :, 1], g.nx)
  @views A_mul_B!(v.u, g.ifftplan, v.solc[:, :, 1])
  @views A_mul_B!(v.v, g.ifftplan, v.solc[:, :, 2])
  @views A_mul_B!(v.p, g.ifftplan, v.solc[:, :, 3])
  mode0apv(v.Z, v.u, v.v, v.p, p, g)
end

function mode0apv(prob::AbstractProblem)
  mode0apv(prob.vars, prob.params, prob.grid)
end

function apv(prob::AbstractProblem)
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


function mode1apv(v::Vars, p::TwoModeParams, g::TwoDGrid)
  v.Z = irfft(v.solr, g.nx)
  mode1apv(v.Z, v, p, g)
end

function mode1apv(v::PassiveAPVVars, p::PassiveAPVParams, g::TwoDGrid)
  @views v.Z = irfft(v.solr[:, :, 1], g.nx)
  mode1apv(v.Z, v, p, g)
end

function mode1apv(prob::AbstractProblem)
  mode1apv(prob.vars, prob.params, prob.grid)
end


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
mode0speed(prob) = sqrt.(prob.vars.U.^2.0 + prob.vars.V.^2.0)

mode1buoyancy(v, p) = real.(im*p.m*v.p - im*p.m*conj.(v.p))
mode1buoyancy(prob::AbstractProblem) = mode1buoyancy(prob.vars, prob.params)

function mode1w(v::PrognosticAPVVars, p::PrognosticAPVParams, g::TwoDGrid)
  @views w = -1.0./p.m.*ifft(g.K.*v.solc[:, :, 1] .+ g.L.*v.solc[:, :, 2])
  real.(w + conj.(w))
end




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



  
""" 
Returns Q^chi, the quadratic mode-1 contribution to barotropic 
available potential vorticity. 
"""
function Qchi(u, v, p, pr::TwoModeParams, g::TwoDGrid)
  pr.m.^2.0./pr.N.^2.0 .* irfft(
       im.*g.Kr .* rfft(real.(v.*conj.(p) .+ conj.(v).*p))
    .- im.*g.Lr .* rfft(real.(u.*conj.(p) .+ conj.(u).*p)), g.nx)
end

function Qchi(v::TwoModeVars, p::TwoModeParams, g::TwoDGrid)
  @views A_mul_B!(v.u, g.ifftplan, v.solc[:, :, 1])
  @views A_mul_B!(v.v, g.ifftplan, v.solc[:, :, 2])
  @views A_mul_B!(v.p, g.ifftplan, v.solc[:, :, 3])

  Qchi(v.u, v.v, v.p, p, g)
end

function Qchi(v::PrognosticAPVVars, p::PrognosticAPVParams, g::TwoDGrid)
  @views A_mul_B!(v.u, g.ifftplan, v.solc[:, :, 1])
  @views A_mul_B!(v.v, g.ifftplan, v.solc[:, :, 2])
  @views A_mul_B!(v.p, g.ifftplan, v.solc[:, :, 3])

  @. v.up = real(v.u*conj(v.p) + conj(v.u)*v.p)
  @. v.vp = real(v.v*conj(v.p) + conj(v.v)*v.p)

  A_mul_B!(v.uph, g.rfftplan, v.up)
  A_mul_B!(v.vph, g.rfftplan, v.vp)

  irfft(p.m^2.0./p.N^2.0.*im.*(g.Kr.*v.vph .- g.Lr.*v.uph), g.nx)
end

function Qchi(prob::AbstractProblem)
  Qchi(prob.vars, prob.params, prob.grid)
end




""" 
Returns chi, the "mode-1-induced" barotropic streamfunction.
"""
function chi(Qchi, g::TwoDGrid)
  -irfft(g.invKKrsq.*rfft(Qchi), g.nx)
end

function chi(v::TwoModeVars, p::TwoModeParams, g::TwoDGrid)
  chi(Qchi(v, p, g), g)
end

function chi(prob::AbstractProblem)
  chi(prob.vars, prob.params, prob.grid)
end




""" 
Returns the velocity field associated with the "mode-1-induced" barotropic 
streamfunction. 
"""
function UVchi(Qchi::AbstractArray, g::TwoDGrid)
  Uchi =  irfft(im.*g.Lr.*g.invKKrsq.*rfft(Qchi), g.nx)
  Vchi = -irfft(im.*g.Kr.*g.invKKrsq.*rfft(Qchi), g.nx)
  Uchi, Vchi
end

function UVchi(v::TwoModeVars, p::TwoModeParams, g::TwoDGrid)
  UVchi(Qchi(v, p, g), g)
end

function UVchi(prob::AbstractProblem)
  UVchi(prob.vars, prob.params, prob.grid)
end




""" 
Returns the barotropic vorticity, Z.
"""
function get_Z(v::PrognosticAPVVars, p::PrognosticAPVParams, g::TwoDGrid)
  calcPsih!(v, p, g)
  -irfft(g.KKrsq.*v.Psih, g.nx)
end

function get_Z(prob::AbstractProblem)
  get_Z(prob.vars, prob.params, prob.grid)
end
 



# Wave-induced flow and potential vorticity ----------------------------------- 
""" 
Calculate the wave-induced streamfunction and velocity fields.
"""
function wave_induced_uv(qw, g::TwoDGrid)

  qwh = rfft(qw)

  Psiwh = g.invKKrsq.*qwh
  uwh   = -im*g.Lr.*Psiwh
  vwh   =  im*g.Kr.*Psiwh

  Psiw = irfft(Psiwh, g.nx)
  uw   = irfft(uwh, g.nx)
  vw   = irfft(vwh, g.nx)

  uw, vw
end

function wave_induced_uv(sig::Real, v::TwoModeVars, p::TwoModeParams, 
  g::TwoDGrid)
  wave_induced_uv(calc_qw(sig, v, p, g), g)
end

function wave_induced_uv(sig, prob::AbstractProblem)
  wave_induced_uv(sig, prob.vars, prob.params, prob.grid)
end

""" 
Returns the wave-induced streamfunction.
"""
function wave_induced_psi(qw, g::TwoDGrid)
  irfft(g.invKKrsq.*rfft(qw), g.nx)
end

function wave_induced_psi(sig::Real, v::TwoModeVars, p::TwoModeParams, 
  g::TwoDGrid)
  wave_induced_psi(calc_qw(sig, v, p, g), g)
end

function wave_induced_psi(sig, prob::AbstractProblem)
  wave_induced_psi(sig, prob.vars, prob.params, prob.grid)
end

""" 
Calculate the wave contribution to PV, qw.
"""
function calc_qw(sig::Real, v::TwoModeVars, p::TwoModeParams, g::TwoDGrid)

  usig, vsig = calc_usigvsig(sig, v, p, g)

  # Non-Jacobian terms
  usig2xx = ifft(-g.K.^2.0.*fft(abs2.(usig)))
  vsig2yy = ifft(-g.L.^2.0.*fft(abs2.(vsig)))
  usigvsigxy = ifft(-g.K.*g.L.*fft(usig.*conj.(vsig) + conj.(usig).*vsig))

  # Assemble contributions
  qw = -real.( 
       2.0*im/sig * (
        FourierFlows.jacobian(conj.(usig), usig, g) 
      + FourierFlows.jacobian(conj.(vsig), vsig, g))
    + p.f/sig^2.0 * (
        FourierFlows.jacobian(conj.(vsig), usig, g) 
      + FourierFlows.jacobian(vsig, conj.(usig), g))
    + p.f/sig^2.0 * (usig2xx + vsig2yy + usigvsigxy)
  )

  qw
end

""" 
Calculate the wave contribution to PV, qw.
"""
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
        FourierFlows.jacobian(conj.(usig), usig, g) 
      + FourierFlows.jacobian(conj.(vsig), vsig, g))
    + p.f/sig^2.0 * (
        FourierFlows.jacobian(conj.(vsig), usig, g) 
      + FourierFlows.jacobian(vsig, conj.(usig), g))
    + p.f/sig^2.0 * (usig2xx + vsig2yy + usigvsigxy)
  )

  qw
end

function calc_qw(sig::Real, prob::AbstractProblem)
  calc_qw(sig, prob.vars, prob.params, prob.grid)
end

function calc_usigvsig(sig, Zh, v::TwoModeVars, p::TwoModeParams, g::TwoDGrid)
  @views @. v.uh = v.solc[:, :, 1]
  @views @. v.vh = v.solc[:, :, 2]
  A_mul_B!(v.u,  g.ifftplan, v.uh)
  A_mul_B!(v.v,  g.ifftplan, v.vh)
  usig = exp(im*sig*v.t) * v.u
  vsig = exp(im*sig*v.t) * v.v
  usig, vsig
end

function calc_usigvsig(sig, v::Vars, p::TwoModeParams, g::TwoDGrid)
  v.Zh .= v.solr
  calc_usigvsig(sig, v.Zh, v, p, g)
end

function calc_usigvsig(sig, v::PassiveAPVVars, p::TwoModeParams, g::TwoDGrid)
  @views @. v.Zh = v.solr[:, :, 1]
  calc_usigvsig(sig, v.Zh, v, p, g)
end

""" 
Returns the wave-induced speed.
"""
function wave_induced_speed(sig, vs::AbstractVars, pr::AbstractParams, 
  g::AbstractGrid)
  uw, vw = wave_induced_uv(sig, vs, pr, g)
  sqrt.(uw.^2.0 + vw.^2.0)
end

function wave_induced_speed(sig, prob::AbstractProblem)
  wave_induced_speed(sig, prob.vars, prob.params, prob.grid)
end

""" 
Returns the wave-induced x-velocity. 
"""
function wave_induced_u(sig, vs::AbstractVars, pr::AbstractParams, 
  g::AbstractGrid)
  wave_induced_uv(sig, vs, pr, g)[1]
end

function wave_induced_u(sig, prob::AbstractProblem)
  wave_induced_u(prob.vars, prob.params, prob.grid)
end

""" 
Returns the wave-induced y-velocity. 
"""
function wave_induced_v(vs, pr, g)
  wave_induced_uv(sig, vs, pr, g)[2]
end

function wave_induced_v(prob::AbstractProblem)
  wave_induced_v(prob.vars, prob.params, prob.grid)
end

""" 
Returns the APV-induced streamfunction.
"""
function apv_induced_psi(Q, g::TwoDGrid)
  PsiQh = -g.invKKrsq.*rfft(Q)
  irfft(PsiQh, g.nx)
end

function apv_induced_psi(v, p, g)
  apv_induced_psi(mode0apv(v, p, g), g)
end

function apv_induced_psi(prob::AbstractProblem)
  apv_induced_psi(prob.vars, prob.params, prob.grid)
end

""" 
Returns the speed of the flow induced by the available potential 
vorticity field. 
"""
function apv_induced_speed(vs, pr, g)
  Q = mode0apv(vs, pr, g)
  PsiQh = -g.invKKrsq.*rfft(Q)
  uQ = irfft(-im*g.Lr.*PsiQh, g.nx)
  vQ = irfft( im*g.Kr.*PsiQh, g.nx)
  sqrt.(uQ.^2.0+vQ.^2.0)
end

function apv_induced_speed(prob::AbstractProblem)
  apv_induced_speed(prob.vars, prob.params, prob.grid)
end

""" 
Return the total Lagrangian-mean flow. 
"""
function lagrangian_mean_uv(sig, vs::AbstractVars, pr::AbstractParams, 
  g::AbstractGrid)
  PsiLh = lagrangian_mean_psih(sig, vs, pr, g)
  uL = irfft(-im*g.Lr.*PsiLh, g.nx)
  vL = irfft( im*g.Kr.*PsiLh, g.nx)
  uL, vL
end

function lagrangian_mean_uv(sig, prob::AbstractProblem)
  lagrangian_mean_uv(sig, prob.vars, prob.params, prob.grid)
end

""" Return the Lagrangian-mean streamfunction. """
function lagrangian_mean_psih(σ, vs::AbstractVars, pr::AbstractParams, 
  g::AbstractGrid)
  Q  = mode0apv(vs, pr, g)
  Qw = calc_qw(σ, vs, pr, g)
  -g.invKKrsq.*rfft(Q-Qw)
end

function lagrangian_mean_psih(σ, prob::AbstractProblem)
  lagrangian_mean_psih(σ, prob.vars, prob.params, prob.grid)
end

function lagrangian_mean_psi(σ, prob::AbstractProblem)
  irfft(lagrangian_mean_psih(σ, prob), prob.grid.nx)
end



# End module
end
