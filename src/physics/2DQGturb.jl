__precompile__()
println("framework.jl uses precompile")

include("../domain.jl")
include("../timesteppers.jl")




# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# B E G I N    M O D U L E    F R A M E W O R K --------------------------------


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
  LC::Array{Complex128, 2}

end

function Params(f0::Float64, nuq::Float64, nuqn::Int, g::Grid)

  # Linear coefficients:
  # Dissipation of mean vorticity
  LC = -nuq * g.KKrsq.^(0.5*nuqn)

  Params(f0, nuq, nuqn, LC)
end




# -----------------------------------------------------------------------------
# Vars ------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# Phyiscal variables, including intermediates for calc of nonlinear terms
type Vars
  t::Float64
  # Solution
  qh::Array{Complex128, 2}

  # Auxiliary vars
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
  # qh = exp.( 2.0*pi*im*rand(g.nkr, g.nl) )

  # Initial condition in the form of two ellipsoid vortices used for comparisson with Matlab code
  ampl = 0.119321207356767;
  ampl = 1.131562576275490e-04;
  qh = ampl*rfft( 200.0*exp.(-((g.X-1).^2-0.4*g.X.*g.Y)./.3^2-(g.Y-1).^2./.5^2)  - 100.0* exp.(-((g.X+1).^2-0.4*g.X.*g.Y)./.3^2-(g.Y+1).^2./.5^2) );


  # qh = rfft(sin.(3*g.X+4*g.Y) + sin.(2*g.X+5*g.Y));
  # qh = rfft(sin(3*g.X+4*g.Y))
  qh[1,1]=0;

  return Vars(t, qh, q, U, V, Uq, Vq, psih, Uh, Vh, Uqh, Vqh)
end



end

# E N D    M O D U L E    F R A M E W O R K ------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------





# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# B E G I N    M O D U L E    S O L V E R --------------------------------------

module Solver
using Domain, Framework, TimeSteppers

export Grid, Vars, Params, ETDRK4TimeStepper, ForwardEulerTimeStepper
export updatevars!, calc_nl!, calc_nl, stepforward!

# -----------------------------------------------------------------------------
# Solver ----------------------------------------------------------------------
# -----------------------------------------------------------------------------

function calc_nl!(
  NLqh::Array{Complex{Float64},2},
  qh::Array{Complex128, 2},
  t::Float64, v::Vars, p::Params, g::Grid)

  # ON NAVID'S LAPTOP A_mul_B! messes up!!

  # A_mul_B!( v.q, g.irfftplan, v.qh )
  v.q = irfft(qh, g.nx)

  v.Uh = +im*real.(g.Lr).*qh.*real.(g.invKKrsq)
  v.Vh = -im*real.(g.Kr).*qh.*real.(g.invKKrsq)

  # A_mul_B!( v.U, g.irfftplan, v.Uh )
  # A_mul_B!( v.V, g.irfftplan, v.Vh )
  v.U = irfft(v.Uh, g.nx)
  v.V = irfft(v.Vh, g.nx)

  # v.Uq = v.U .* v.q
  # v.Vq = v.V .* v.q

  # A_mul_B!( v.Uqh, g.rfftplan, v.Uq )
  # A_mul_B!( v.Vqh, g.rfftplan, v.Vq )
  v.Uqh = rfft(v.U .* v.q)
  v.Vqh = rfft(v.V .* v.q)

  NLqh = -im.*g.Kr.*v.Uqh -im.*g.Lr.*v.Vqh
  # dealias!(NLqh, g)
end


function calc_nl(
  NLqh::Array{Complex{Float64},2},
  qh::Array{Complex128, 2},
  t::Float64, v::Vars, p::Params, g::Grid)

  # ON NAVID'S LAPTOP A_mul_B! messes up!!

  # A_mul_B!( v.q, g.irfftplan, v.qh )
  v.q = irfft(qh, g.nx)

  v.Uh .= +im*real.(g.Lr).*qh.*real.(g.invKKrsq)
  v.Vh .= -im*real.(g.Kr).*qh.*real.(g.invKKrsq)

  # A_mul_B!( v.U, g.irfftplan, v.Uh )
  # A_mul_B!( v.V, g.irfftplan, v.Vh )
  v.U = irfft(v.Uh, g.nx)
  v.V = irfft(v.Vh, g.nx)

  v.Uq .= v.U .* v.q
  v.Vq .= v.V .* v.q

  # A_mul_B!( v.Uqh, g.rfftplan, v.Uq )
  # A_mul_B!( v.Vqh, g.rfftplan, v.Vq )
  v.Uqh = rfft(v.Uq)
  v.Vqh = rfft(v.Vq)

  NLqh .= -im.*g.Kr.*v.Uqh -im.*g.Lr.*v.Vqh
  # dealias!(NLqh, g)

  return NLqh
end



function stepforward!(nsteps::Int,
  qts::ForwardEulerTimeStepper,
  v::Vars, p::Params, g::Grid)

  for step = 1:nsteps
    # println(v.qh[10,20])
    # calc_nl!(qts.NL, v.qh, v.t, v, p, g)
    # println("why is qts.NL not updated?")
    # println(qts.NL[10,20])

    qts.NL = calc_nl(qts.NL, v.qh, v.t, v, p, g)

    v.qh = v.qh + qts.dt * (qts.NL + qts.LC.*v.qh)
  end
end


function stepforward!(nsteps::Int,
  qts::ETDRK4TimeStepper,
  v::Vars, p::Params, g::Grid)

  for step = 1:nsteps
    # add the ETDRK4 time-stepping scheme
  end

end

# -----------------------------------------------------------------------------
# Helper functions ------------------------------------------------------------
# -----------------------------------------------------------------------------

function updatevars!(v::Vars, p::Params, g::Grid)

  v.q = irfft(v.qh, g.nx)
  v.psih = -v.qh .* g.invKKrsq
  v.U = -irfft(im*g.Lr.*v.psih, g.nx)
  v.V =  irfft(im*g.Kr.*v.psih, g.nx)

end

end


# E N D    M O D U L E    S O L V E R ------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
