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
  kap::Float64  # Potential vorticity hyperdiffusivity
  nkap::Int     # Order of Potential vorticity hyperdiffusity
  nu::Float64   # Wave hyperviscosity
  nnu::Int      # Order of wave hyperviscosity
  eta::Float64  # Wave dispersivity
  f::Float64    # Planetary vorticity
end

function Params(kap::Real, nu::Real, eta::Real, f::Real)
  nkap = nnu = 4
  Params(kap, nkap, nu, nnu, eta, f)
end




# E Q U A T I O N S ----------------------------------------------------------- 
type Equation <: AbstractEquation
  LCc::Array{Complex{Float64}, 2}  # Coeff of the complex eqn's linear part
  LCr::Array{Complex{Float64}, 2}  # Coeff of the real eqn's linear part
  calcNL!::Function                # Function to calculate eqn's nonlinear part
end

function Equation(p::Params, g::TwoDGrid)
  LCr = -p.kap * g.KKrsq.^(0.5*p.nkap)
  LCc = -im*0.5*p.eta*g.KKsq - p.nu*g.KKsq.^(0.5*p.nnu)

  # Function calcNL! is defined below.
  Equation(LCc, LCr, calcNL!)
end




# V A R S --------------------------------------------------------------------- 
type Vars <: AbstractVars

  t::Float64
  solr::Array{Complex128, 2} # solr = q
  solc::Array{Complex128, 2} # solc = phi

  # Auxiliary vorticity vars
  q::Array{Float64, 2}
  U::Array{Float64, 2}
  V::Array{Float64, 2}
  zeta::Array{Float64, 2}
  psi::Array{Float64, 2}
  Uq::Array{Float64, 2}
  Vq::Array{Float64, 2}
  modphi::Array{Float64, 2}
  jacphi::Array{Float64, 2}

  # Auxiliary wave vars
  phi::Array{Complex{Float64}, 2}
  Uphi::Array{Complex{Float64}, 2}
  Vphi::Array{Complex{Float64}, 2}
  zetaphi::Array{Complex{Float64}, 2}
  
  # Vorticity transforms
  qh::Array{Complex{Float64}, 2}
  Uh::Array{Complex{Float64}, 2}
  Vh::Array{Complex{Float64}, 2}
  zetah::Array{Complex{Float64}, 2}
  psih::Array{Complex{Float64}, 2}
  Uqh::Array{Complex{Float64}, 2}
  Vqh::Array{Complex{Float64}, 2}
  modphih::Array{Complex{Float64}, 2}
  jacphih::Array{Complex{Float64}, 2}

  # Wave transforms
  phih::Array{Complex{Float64}, 2}
  Uphih::Array{Complex{Float64}, 2}
  Vphih::Array{Complex{Float64}, 2}
  zetaphih::Array{Complex{Float64}, 2}

end

function Vars(g::TwoDGrid)

  # Initialize with t=0
  t = 0.0
  solr  = zeros(Complex{Float64}, g.nkr, g.nl)
  solc  = zeros(Complex{Float64}, g.nk, g.nl)

  # Auxiliary vorticity vars
  q      = zeros(Float64, g.nx, g.ny)
  U      = zeros(Float64, g.nx, g.ny)
  V      = zeros(Float64, g.nx, g.ny)
  zeta   = zeros(Float64, g.nx, g.ny)
  psi    = zeros(Float64, g.nx, g.ny)
  Uq     = zeros(Float64, g.nx, g.ny)
  Vq     = zeros(Float64, g.nx, g.ny)
  modphi = zeros(Float64, g.nx, g.ny)
  jacphi = zeros(Float64, g.nx, g.ny)
  
  # Auxiliary wave vars
  phi     = zeros(Complex{Float64}, g.nx, g.ny)
  Uphi    = zeros(Complex{Float64}, g.nx, g.ny)
  Vphi    = zeros(Complex{Float64}, g.nx, g.ny)
  zetaphi = zeros(Complex{Float64}, g.nx, g.ny)

  # Transforms
  qh       = zeros(Complex{Float64}, g.nkr, g.nl)
  Uh       = zeros(Complex{Float64}, g.nkr, g.nl)
  Vh       = zeros(Complex{Float64}, g.nkr, g.nl)
  zetah    = zeros(Complex{Float64}, g.nkr, g.nl)
  psih     = zeros(Complex{Float64}, g.nkr, g.nl)
  Uqh      = zeros(Complex{Float64}, g.nkr, g.nl)
  Vqh      = zeros(Complex{Float64}, g.nkr, g.nl)
  modphih  = zeros(Complex{Float64}, g.nkr, g.nl)
  jacphih  = zeros(Complex{Float64}, g.nkr, g.nl)

  phih     = zeros(Complex{Float64}, g.nk, g.nl)
  Uphih    = zeros(Complex{Float64}, g.nk, g.nl)
  Vphih    = zeros(Complex{Float64}, g.nk, g.nl)
  zetaphih = zeros(Complex{Float64}, g.nk, g.nl)
  
  return Vars(t, solr, solc, 
    q, U, V, zeta, psi, Uq, Vq, modphi, jacphi,
    phi, Uphi, Vphi, zetaphi,
    qh, Uh, Vh, zetah, psih, Uqh, Vqh, modphih, jacphih,
    phih, Uphih, Vphih, zetaphih,
    )
end




# S O L V E R S ---------------------------------------------------------------

function calcNL!(
  NLc::Array{Complex{Float64}, 3},  NLr::Array{Complex{Float64}, 2}, 
  solc::Array{Complex{Float64}, 3}, solr::Array{Complex{Float64}, 2}, 
  t::Float64, v::Vars, p::Params, g::TwoDGrid)
  

  # Calcs associated with wave contrib to PV
  @. v.phixh = im*g.K*solc
  @. v.phiyh = im*g.L*solc

  A_mul_B!(v.phix, g.ifftplan, v.phixh)
  A_mul_B!(v.phiy, g.ifftplan, v.phiyh)

  @. v.modphi = abs2(v.phi)
  @. v.jacphi = real(im*conj(v.phix)*v.phiy - im*conj(v.phiy)*v.phix)

  A_mul_B!(v.modphih, g.rfftplan, v.modphi)
  A_mul_B!(v.jacphih, g.rfftplan, v.jacphi)


  # Mean flow calcs and inverse transforms
  @. v.qh = solr        # Necessary because irfft destroys input
  @. v.zetah = (v.qh 
    + 0.25/p.f*g.KKrsq*v.modphih - 0.5/p.f*v.jacphih
  )

  @. v.psih = -g.invKKrsq*v.zetah
  @. v.Uh = -im*g.Lr*v.psih
  @. v.Vh =  im*g.Kr*v.psih

  A_mul_B!(v.U, g.irfftplan, v.Uh)
  A_mul_B!(v.V, g.irfftplan, v.Vh)
  A_mul_B!(v.zeta, g.irfftplan, v.zetah)

  A_mul_B!(v.q, g.irfftplan, solr)
  A_mul_B!(v.phi, g.ifftplan, solc)


  # Multiplies and forward transforms
  @. v.Uq = v.U*v.q
  @. v.Vq = v.V*v.q

  @. v.Uphi    = v.U*v.phi
  @. v.Vphi    = v.V*v.phi
  @. v.zetaphi = v.zeta*v.phi

  A_mul_B!(v.Uqh, g.rfftplan, v.Uq)
  A_mul_B!(v.Vqh, g.rfftplan, v.Vq)

  A_mul_B!(v.Uphih, g.fftplan, v.Uphi)
  A_mul_B!(v.Vphih, g.fftplan, v.Vphi)
  A_mul_B!(v.zetaphih, g.fftplan, v.zetaphi)


  # Nonlinear terms
  @. NLr = -im*g.Kr*v.Uqh - im*g.Lr*v.Vqh
  @. NLc = -im*g.K*v.Uphi - im*g.L*v.Vphi - 0.5*im*v.zetaphi

  nothing
end




# H E L P E R   F U N C T I O N S --------------------------------------------- 
function updatevars!(v::Vars, p::Params, g::TwoDGrid)

  @. v.qh = v.solr
  @. v.phih = v.solc 

  # Wave calcs and inverse transforms
  @. v.phixh = im*g.K*v.phih
  @. v.phiyh = im*g.L*v.phih

  A_mul_B!(v.phi, g.ifftplan, v.phih)
  A_mul_B!(v.phix, g.ifftplan, v.phixh)
  A_mul_B!(v.phiy, g.ifftplan, v.phiyh)

  # Wave parts of PV. J(phi', phi) = conj(phi'x)*phiy) - conj(phi'y)*phix)
  @. v.modphi = v.phi*conj(v.phi)
  @. v.jacphi = real(im*conj(v.phix)*v.phiy - im*conj(v.phiy)*v.phix)

  A_mul_B!(v.modphih, g.rfftplan, v.modphi)
  A_mul_B!(v.jacphih, g.rfftplan, v.jacphi)

# Vorticity calcs and inverse transforms
  @. v.zetah = (v.qh 
    + 0.25/p.f*g.KKrsq*v.modphih - 0.5/p.f*v.jacphih
  )

  @. v.psih = -g.invKKrsq*v.zetah
  @. v.Uh = -im*g.Lr*v.psih
  @. v.Vh =  im*g.Kr*v.psih

  # We don't use A_mul_B here because irfft destroys its input.
  # A_mul_B!(v.q, g.irfftplan, v.qh)
  v.q = irfft(v.qh, g.nx)
  v.U = irfft(v.Uh, g.nx) 
  v.V = irfft(v.Vh, g.nx) 
  v.zeta = irfft(v.zetah, g.nx) 

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
