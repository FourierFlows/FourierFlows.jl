__precompile__()
println("framework.jl uses precompile")

include("domain.jl")

module Framework
using Domain

export Grid, Vars, Params


# -----------------------------------------------------------------------------
# Params ----------------------------------------------------------------------
# -----------------------------------------------------------------------------
immutable Params

  # Physical parameters
  f0::Float64                           # Inertial period
  nuq::Float64                          # Vorticity viscosity
  nuqn::Int                             # Vorticity hyperviscous order

  # Linear left hand sides of the qh and ah equation
  LCq::Array{Complex128, 2}

end

function Params(f0::Float64, nuq::Float64, nuqn::Int, g::Grid)

  # Linear coefficients:
  # Dissipation of mean vorticity
  LCq = -nuq * g.KKrsq.^(0.5*nuqn)

  Params(f0, nuq, nuqn, LCq)
end




# -----------------------------------------------------------------------------
# Vars ------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# Phyiscal variables, including intermediates for calc of nonlinear terms
type Vars
  t::Float64
  # Solution
  qh::Array{Complex128, 2}

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

end

function Vars(p::Params, g::Grid)
  t = 0.0

  # Solution
  qh = zeros(Complex128, g.nkr, g.nl)

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

  # Random initial condition
  qh = exp.( 2.0*pi*im*rand(g.nkr, g.nl) )

  ampl = 0.119321207356767;
  ampl = 1.131562576275490e-04;
  qh = ampl*rfft( 200.0*exp.(-((g.X-1).^2-0.4*g.X.*g.Y)./.3^2-(g.Y-1).^2./.5^2)  - 100.0* exp.(-((g.X+1).^2-0.4*g.X.*g.Y)./.3^2-(g.Y+1).^2./.5^2) );


  # qh = rfft(sin.(3*g.X+4*g.Y) + sin.(2*g.X+5*g.Y));
  # qh = rfft(sin(3*g.X+4*g.Y))
  qh[1,1]=0;

  return Vars(t, qh, q, U, V, Uq, Vq, psih, Uh, Vh, Uqh, Vqh)
end




end
