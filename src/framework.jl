include("domain.jl")

module Framework
using Domain

export Grid, Vars, Params 


# ----------------------------------------------------------------------------- 
# Params ---------------------------------------------------------------------- 
# ----------------------------------------------------------------------------- 
struct Params

  # Physical parameters
  f0::Float64                           # Inertial period
  sig::Float64                          # Wave frequency
  kap::Float64                          # Modewise wavenumber
  alpha::Float64                        # Frequency parameter
  nuq::Float64                          # Vorticity viscosity
  nua::Float64                          # Wavefield viscosity
  nuqn::Int                             # Vorticity hyperviscous order
  nuan::Int                             # Wave hyperviscous order

  # Linear left hand sides of the qh and ah equation
  LCq::Array{Complex128, 2}
  LCa::Array{Complex128, 2}

  # Computational convenience attributes
  E::Array{Complex128, 2}
  invE::Array{Complex128, 2}

  invEk::Array{Complex128, 2}
  invEl::Array{Complex128, 2}

  invEf0k::Array{Complex128, 2}
  invEf0l::Array{Complex128, 2}

  invE2sigk::Array{Complex128, 2}
  invE2sigl::Array{Complex128, 2}

end

function Params(f0::Float64, sig::Float64, kap::Float64,
  nuq::Float64, nua::Float64, nuqn::Int, nuan::Int, g::Grid)

  alpha = sig^2/f0^2 - 1.0

  E = -0.5*alpha*( g.KKsq + kap^2.0*(4.0 + 3.0*alpha) )
  invE = 1.0 ./ E

  invEk = g.K ./ E
  invEl = g.L ./ E

  invEf0k = f0*g.K ./ E
  invEf0l = f0*g.L ./ E

  invE2sigk = 2.0*sig/f0^2 * g.K ./ E
  invE2sigl = 2.0*sig/f0^2 * g.L ./ E

  # Linear coefficients:
  # Dissipation of mean vorticity
  LCq = -nuq * g.KKrsq.^(0.5*nuqn)

  # Hyperdissipation and dispersion of wave field
  diss = -nua * g.KKsq.^(0.5*nuan)
  disp = (im*alpha*sig) * invE .* (g.KKsq - alpha*kap^2.0)
  LCa = diss + disp

  Params(f0, sig, kap, alpha, nuq, nua, nuqn, nuan, LCq, LCa, 
    E, invE, invEk, invEl, invEf0k, invEf0l, invE2sigk, invE2sigl)
end




# ----------------------------------------------------------------------------- 
# Vars ------------------------------------------------------------------------ 
# ----------------------------------------------------------------------------- 

# Phyiscal variables, including intermediates for calc of nonlinear terms
mutable struct Vars
  t::Float64

  # Solution
  qh::Array{Complex128, 2}
  Ah::Array{Complex128, 2}

  # Vorticity auxiliary vars
  q::Array{Float64, 2}
  U::Array{Float64, 2}
  V::Array{Float64, 2}
  Uq::Array{Float64, 2}
  Vq::Array{Float64, 2}

  psih::Array{Complex128, 2}
  Uh::Array{Complex128, 2}
  Vh::Array{Complex128, 2}
  Uqh::Array{Complex128, 2}
  Vqh::Array{Complex128, 2}

  # Wave auxiliary vars
  A::Array{Complex128, 2}
  Ax::Array{Complex128, 2}
  Ay::Array{Complex128, 2}
  Axx::Array{Complex128, 2}
  Ayy::Array{Complex128, 2}
  Axy::Array{Complex128, 2}
  EA::Array{Complex128, 2}
  UEA::Array{Complex128, 2}
  VEA::Array{Complex128, 2}

  qrefxA::Array{Complex128, 2}
  qrefyA::Array{Complex128, 2}
  UVjacxA::Array{Complex128, 2}
  UVjacyA::Array{Complex128, 2}

  Axh::Array{Complex128, 2}
  Ayh::Array{Complex128, 2}
  Axxh::Array{Complex128, 2}
  Ayyh::Array{Complex128, 2}
  Axyh::Array{Complex128, 2}
  EAh::Array{Complex128, 2}
  UEAh::Array{Complex128, 2}
  VEAh::Array{Complex128, 2}

  qrefxAh::Array{Complex128, 2}
  qrefyAh::Array{Complex128, 2}
  UVjacxAh::Array{Complex128, 2}
  UVjacyAh::Array{Complex128, 2}

  uh::Array{Complex128, 2}
  vh::Array{Complex128, 2}
  u::Array{Float64, 2}
  v::Array{Float64, 2}
end

function Vars(p::Params, g::Grid)
  t = 0.0

  # Solution
  qh = zeros(Complex128, g.nkr, g.nl)
  Ah = zeros(Complex128, g.nk, g.nl)

  # Vorticity auxiliary vars
  q  = zeros(Float64, g.nx, g.ny)
  U  = zeros(Float64, g.nx, g.ny)
  V  = zeros(Float64, g.nx, g.ny)
  Uq = zeros(Float64, g.nx, g.ny)
  Vq = zeros(Float64, g.nx, g.ny)

  psih = zeros(Complex128, g.nkr, g.nl)
  Uh   = zeros(Complex128, g.nkr, g.nl)
  Vh   = zeros(Complex128, g.nkr, g.nl)
  Uqh  = zeros(Complex128, g.nkr, g.nl)
  Vqh  = zeros(Complex128, g.nkr, g.nl)
  
  # Wave auxiliary vars
  A  = zeros(Complex128, g.nx, g.ny)
  Ax = zeros(Complex128, g.nx, g.ny)
  Ay = zeros(Complex128, g.nx, g.ny)

  Axx = zeros(Complex128, g.nx, g.ny)
  Ayy = zeros(Complex128, g.nx, g.ny)
  Axy = zeros(Complex128, g.nx, g.ny)
  EA  = zeros(Complex128, g.nx, g.ny)
  UEA = zeros(Complex128, g.nx, g.ny)
  VEA = zeros(Complex128, g.nx, g.ny)

  qrefxA  = zeros(Complex128, g.nx, g.ny)
  qrefyA  = zeros(Complex128, g.nx, g.ny)
  UVjacxA = zeros(Complex128, g.nx, g.ny)
  UVjacyA = zeros(Complex128, g.nx, g.ny)

  Axh  = zeros(Complex128, g.nk, g.nl)
  Ayh  = zeros(Complex128, g.nk, g.nl)
  Axxh = zeros(Complex128, g.nk, g.nl)
  Ayyh = zeros(Complex128, g.nk, g.nl)
  Axyh = zeros(Complex128, g.nk, g.nl)
  EAh  = zeros(Complex128, g.nk, g.nl)
  UEAh = zeros(Complex128, g.nk, g.nl)
  VEAh = zeros(Complex128, g.nk, g.nl)

  qrefxAh  = zeros(Complex128, g.nk, g.nl)
  qrefyAh  = zeros(Complex128, g.nk, g.nl)
  UVjacxAh = zeros(Complex128, g.nk, g.nl)
  UVjacyAh = zeros(Complex128, g.nk, g.nl)

  uh = zeros(Complex128, g.nk, g.nl)
  vh = zeros(Complex128, g.nk, g.nl)
  u  = zeros(Float64, g.nx, g.ny)
  v  = zeros(Float64, g.nx, g.ny)

  # Random initial condition
  qh = exp.( 2.0*pi*im*rand(g.nkr, g.nl) )

  # Plane wave
  k0 = 2.0*pi/g.Lx * round(g.Lx*sqrt(p.alpha)*p.kap/(2.0*pi))
  Ah = fft(exp.( im*k0*g.X ))

  return Vars(t, qh, Ah, q, U, V, Uq, Vq, psih, Uh, Vh, Uqh, Vqh, 
    A, Ax, Ay, Axx, Ayy, Axy, EA, UEA, VEA, qrefxA, qrefyA, UVjacxA, UVjacyA,
    Axh, Ayh, Axxh, Ayyh, Axyh, EAh, UEAh, VEAh, qrefxAh, qrefyAh, 
    UVjacxAh, UVjacyAh, uh, vh, u, v)
end




end
