__precompile__()


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# T W O M O D E B O U S S I N E S Q >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
module TwoModeBoussinesq

using FourierFlows

export Grid,
       Params,
       Vars,
       Equation

export set_q!, updatevars!


# P A R A M S ----------------------------------------------------------------- 
type Params <: AbstractParams
  nu0::Float64                    # Mode-0 viscosity
  nu0n::Int                       # Mode-0 hyperviscous order
  nu1::Float64                    # Mode-1 viscosity
  nu1n::Int                       # Mode-1 hyperviscous order
  m::Float64                      # Mode-one wavenumber
  f::Float64                      # Planetary vorticity
  N::Float64                      # Buoyancy frequency
end

function Params(nu::Real, nun::Int)
  nu0 = nu
  nu1 = nu
  nu0n = nun
  nu1n = nun
  m = 4.0/100.0
  f = 1.0
  N = 10.0
  Params(nu0, nu0n, nu1, nu1n, m, f, N)
end



# E Q U A T I O N S ----------------------------------------------------------- 
type Equation <: AbstractEquation
  LCc::Array{Complex{Float64}, 3}  # Element-wise coeff of the eqn's linear part
  LCr::Array{Complex{Float64}, 2}  # Element-wise coeff of the eqn's linear part
  calcNL!::Function               # Function to calculate eqn's nonlinear part
end

function Equation(p::Params, g::TwoDGrid)
  LCr = -p.nu0 * g.KKrsq.^(0.5*p.nu0n)

  LCc = zeros(g.nx, g.ny, 3)
  LCc[:, :, 1] = -p.nu1 * g.KKsq.^(0.5*p.nu1n)
  LCc[:, :, 2] = -p.nu1 * g.KKsq.^(0.5*p.nu1n)
  #LCc[:, :, 3] = -p.nu1 * g.KKsq.^(0.5*p.nu1n)

  # Function calcNL! is defined below.
  Equation(LCc, LCr, calcNL!)
end




# V A R S --------------------------------------------------------------------- 
type Vars <: AbstractVars

  # q = zeta
  # u0, v0, psi0 = U, V, psi
  # u1, v1, w1, p1 = u, v, w, p

  t::Float64
  solr::Array{Complex128, 2}
  solc::Array{Complex128, 3}

  # Auxiliary zeroth-mode vars
  q::Array{Float64, 2}
  U::Array{Float64, 2}
  V::Array{Float64, 2}
  Uq::Array{Float64, 2}
  Vq::Array{Float64, 2}
  Ux::Array{Float64, 2}
  Uy::Array{Float64, 2}
  Vx::Array{Float64, 2}
  Vy::Array{Float64, 2}
  uw::Array{Float64, 2}
  vw::Array{Float64, 2}
  psi::Array{Float64, 2}

  # Auxiliary first-mode vars
  u::Array{Complex{Float64}, 2}
  v::Array{Complex{Float64}, 2}
  w::Array{Complex{Float64}, 2}
  p::Array{Complex{Float64}, 2}

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
  qh::Array{Complex{Float64}, 2}
  Uh::Array{Complex{Float64}, 2}
  Vh::Array{Complex{Float64}, 2}
  Uqh::Array{Complex{Float64}, 2}
  Vqh::Array{Complex{Float64}, 2}
  Uxh::Array{Complex{Float64}, 2}
  Uyh::Array{Complex{Float64}, 2}
  Vxh::Array{Complex{Float64}, 2}
  Vyh::Array{Complex{Float64}, 2}
  uwh::Array{Complex{Float64}, 2}
  vwh::Array{Complex{Float64}, 2}
  psih::Array{Complex{Float64}, 2}

  # First-mode transforms
  uh::Array{Complex{Float64}, 2}
  vh::Array{Complex{Float64}, 2}
  wh::Array{Complex{Float64}, 2}
  ph::Array{Complex{Float64}, 2}

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
  solr  = zeros(Complex{Float64}, g.nkr, g.nl)
  solc  = zeros(Complex{Float64}, g.nk, g.nl, 3)

  # Auxiliary zeroth-mode vars
  q    = zeros(Float64, g.nx, g.ny)
  U    = zeros(Float64, g.nx, g.ny)
  V    = zeros(Float64, g.nx, g.ny)
  Uq   = zeros(Float64, g.nx, g.ny)
  Vq   = zeros(Float64, g.nx, g.ny)
  Ux   = zeros(Float64, g.nx, g.ny)
  Uy   = zeros(Float64, g.nx, g.ny)
  Vx   = zeros(Float64, g.nx, g.ny)
  Vy   = zeros(Float64, g.nx, g.ny)
  uw   = zeros(Float64, g.nx, g.ny)
  vw   = zeros(Float64, g.nx, g.ny)
  psi  = zeros(Float64, g.nx, g.ny)
  
  # Auxiliary first-mode vars
  u    = zeros(Complex{Float64}, g.nx, g.ny)
  v    = zeros(Complex{Float64}, g.nx, g.ny)
  w    = zeros(Complex{Float64}, g.nx, g.ny)
  p    = zeros(Complex{Float64}, g.nx, g.ny)

  Uu   = zeros(Complex{Float64}, g.nx, g.ny)
  Uv   = zeros(Complex{Float64}, g.nx, g.ny)
  Up   = zeros(Complex{Float64}, g.nx, g.ny)
  Vu   = zeros(Complex{Float64}, g.nx, g.ny)
  Vv   = zeros(Complex{Float64}, g.nx, g.ny)
  Vp   = zeros(Complex{Float64}, g.nx, g.ny)

  uUx  = zeros(Complex{Float64}, g.nx, g.ny)
  uVx  = zeros(Complex{Float64}, g.nx, g.ny)
  vUy  = zeros(Complex{Float64}, g.nx, g.ny)
  vVy  = zeros(Complex{Float64}, g.nx, g.ny)

  # Transforms
  qh   = zeros(Complex{Float64}, g.nkr, g.nl)
  Uh   = zeros(Complex{Float64}, g.nkr, g.nl)
  Vh   = zeros(Complex{Float64}, g.nkr, g.nl)
  Uqh  = zeros(Complex{Float64}, g.nkr, g.nl)
  Vqh  = zeros(Complex{Float64}, g.nkr, g.nl)
  Uxh  = zeros(Complex{Float64}, g.nkr, g.nl)
  Uyh  = zeros(Complex{Float64}, g.nkr, g.nl)
  Vxh  = zeros(Complex{Float64}, g.nkr, g.nl)
  Vyh  = zeros(Complex{Float64}, g.nkr, g.nl)
  uwh  = zeros(Complex{Float64}, g.nkr, g.ny)
  vwh  = zeros(Complex{Float64}, g.nkr, g.ny)
  psih = zeros(Complex{Float64}, g.nkr, g.nl)

  uh   = zeros(Complex{Float64}, g.nk, g.nl)
  vh   = zeros(Complex{Float64}, g.nk, g.nl)
  wh   = zeros(Complex{Float64}, g.nk, g.nl)
  ph   = zeros(Complex{Float64}, g.nk, g.nl)
  
  Uuh  = zeros(Complex{Float64}, g.nk, g.nl)
  Uvh  = zeros(Complex{Float64}, g.nk, g.nl)
  Uph  = zeros(Complex{Float64}, g.nk, g.nl)
  Vuh  = zeros(Complex{Float64}, g.nk, g.nl)
  Vvh  = zeros(Complex{Float64}, g.nk, g.nl)
  Vph  = zeros(Complex{Float64}, g.nk, g.nl)

  uUxh = zeros(Complex{Float64}, g.nk, g.nl)
  uVxh = zeros(Complex{Float64}, g.nk, g.nl)
  vUyh = zeros(Complex{Float64}, g.nk, g.nl)
  vVyh = zeros(Complex{Float64}, g.nk, g.nl)

  return Vars(t, solr, solc, 
    q, U, V, Uq, Vq, Ux, Uy, Vx, Vy, uw, vw, psi, 
    u, v, w, p, Uu, Uv, Up, Vu, Vv, Vp, uUx, uVx, vUy, vVy,
    qh, Uh, Vh, Uqh, Vqh, Uxh, Uyh, Vxh, Vyh, uwh, vwh, psih, 
    uh, vh, wh, ph, Uuh, Uvh, Uph, Vuh, Vvh, Vph, uUxh, uVxh, vUyh, vVyh,
    )
end




# S O L V E R S ---------------------------------------------------------------

function calcNL!(
  NLc::Array{Complex{Float64}, 3},  NLr::Array{Complex{Float64}, 2}, 
  solc::Array{Complex{Float64}, 3}, solr::Array{Complex{Float64}, 2}, 
  t::Float64, v::Vars, p::Params, g::TwoDGrid)
  

  # This copy is necessary because calling A_mul_B(v.q, g.irfftplan, sol) 
  # a few lines below destroys sol when using Julia's FFTW.
  v.qh .= solr
  A_mul_B!(v.q, g.irfftplan, solr)


  # ---------------------------------------------------------------------------   
  # Zeroth-mode calcs and inverse transforms
  @. v.psih = -g.invKKrsq*v.qh

  @. v.Uh = -im*g.Lr*v.psih
  @. v.Vh =  im*g.Kr*v.psih

  @. v.Uxh = im*g.Kr*v.Uh
  @. v.Vxh = im*g.Kr*v.Vh

  @. v.Uyh = im*g.Lr*v.Uh
  @. v.Vyh = im*g.Lr*v.Vh
 
  A_mul_B!(v.U, g.irfftplan, v.Uh)
  A_mul_B!(v.V, g.irfftplan, v.Vh)

  A_mul_B!(v.Ux, g.irfftplan, v.Uxh)
  A_mul_B!(v.Uy, g.irfftplan, v.Uyh)
  A_mul_B!(v.Vx, g.irfftplan, v.Vxh)
  A_mul_B!(v.Vy, g.irfftplan, v.Vyh)


  # First mode calcs and inverse transforms
  @views A_mul_B!(v.u, g.ifftplan, solc[:, :, 1])
  @views A_mul_B!(v.v, g.ifftplan, solc[:, :, 2])
  @views A_mul_B!(v.p, g.ifftplan, solc[:, :, 3])

  @views @. v.wh = -1.0/p.m*g.K*solc[:, :, 1] + g.L*solc[:, :, 2]
  A_mul_B!(v.w, g.ifftplan, v.wh)


  # ---------------------------------------------------------------------------   
  # Zeroth-mode multiplies and forward transforms
  @. v.Uq = v.U*v.q
  @. v.Vq = v.V*v.q

  @. v.uw = real(im*p.m * v.u*conj(v.w) - im*p.m * conj(v.u)*v.w )
  @. v.vw = real(im*p.m * v.v*conj(v.w) - im*p.m * conj(v.v)*v.w )
               
  A_mul_B!(v.Uqh, g.rfftplan, v.Uq)
  A_mul_B!(v.Vqh, g.rfftplan, v.Vq)

  A_mul_B!(v.uwh, g.rfftplan, v.uw)
  A_mul_B!(v.vwh, g.rfftplan, v.vw)

  
  # ---------------------------------------------------------------------------   
  # First mode multiplies and forward transforms
  @. v.Uu = v.U*v.u
  @. v.Vu = v.V*v.u
  @. v.Uv = v.U*v.v
  @. v.Vv = v.V*v.v
  @. v.Up = v.U*v.p
  @. v.Vp = v.V*v.p

  @. v.uUx = v.u*v.Ux
  @. v.uVx = v.u*v.Vx
  @. v.vUy = v.v*v.Uy
  @. v.vVy = v.v*v.Vy

  A_mul_B!(v.Uuh, g.fftplan, v.Uu)
  A_mul_B!(v.Uvh, g.fftplan, v.Uv)
  A_mul_B!(v.Vuh, g.fftplan, v.Vu)
  A_mul_B!(v.Vvh, g.fftplan, v.Vv)
  A_mul_B!(v.Uph, g.fftplan, v.Up)
  A_mul_B!(v.Vph, g.fftplan, v.Vp)

  A_mul_B!(v.uUxh, g.ifftplan, v.uUx)
  A_mul_B!(v.uVxh, g.ifftplan, v.uVx)
  A_mul_B!(v.vUyh, g.ifftplan, v.vUy)
  A_mul_B!(v.vVyh, g.ifftplan, v.vVy)


  # ---------------------------------------------------------------------------   
  # Zeroth-mode nonlinear term
  @. NLr = ( - im*g.Kr*v.Uqh - im*g.Lr*v.Vqh
  #           + im*g.Lr*v.uwh + im*g.Kr*v.vwh
  )

  # First-mode nonlinear terms
  @views @. NLc[:, :, 1] = ( p.f*v.vh - im*g.K*solc[:, :, 3]
    - im*g.K*v.Uuh - im*g.L*v.Vuh - v.uUxh - v.vUyh
  )

  @views @. NLc[:, :, 2] = ( -p.f*v.uh - im*g.L*solc[:, :, 3]
    - im*g.K*v.Uvh - im*g.L*v.Vvh - v.uVxh - v.vVyh
  )

  @views @. NLc[:, :, 3] = ( im*p.N^2.0/p.m*v.wh
    - im*g.K*v.Uph - im*g.L*v.Vph
  )

  #NLc[:, :, 1] .=   p.f  .* v.vh .- im .* g.K .* solc[:, :, 3]
  #NLc[:, :, 2] .= (-p.f) .* v.uh .- im .* g.L .* solc[:, :, 3]
  #NLc[:, :, 3] .= im .* p.N.^2.0 ./ p.m .* v.wh
  
  nothing
end




# H E L P E R   F U N C T I O N S --------------------------------------------- 
function updatevars!(v::Vars, p::Params, g::TwoDGrid)

  v.qh .= v.solr

  # We don't use A_mul_B here because irfft destroys its input.
  # A_mul_B!(v.q, g.irfftplan, v.qh)
  v.q = irfft(v.qh, g.nx)

  @. v.psih =         -g.invKKrsq*v.qh
  @. v.Uh   =  im*g.Lr*g.invKKrsq*v.qh
  @. v.Vh   = -im*g.Kr*g.invKKrsq*v.qh
 
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




function set_q!(v::Vars, p::Params, g::TwoDGrid, q)
  # Set zeroth mode vorticity
  A_mul_B!(v.solr, g.rfftplan, q)
  updatevars!(v, p, g)
  nothing
end

function set_uvp!(vs::Vars, pr::Params, g::TwoDGrid, u, v, p)
  uh = fft(u)
  vh = fft(v)
  ph = fft(p)

  vs.solc[:, :, 1] .= uh
  vs.solc[:, :, 2] .= vh
  vs.solc[:, :, 3] .= ph

  updatevars!(vs, pr, g)
  nothing
end




end
# E N D   T W O M O D E B O U S S I N E S Q >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
