__precompile__()

include("../fourierflows.jl")

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# B E G I N    M O D U L E    T W O D T U R B ----------------------------------

module TwoDTurb

using FourierFlowTypes, Domains, TimeSteppers

export Grid, Vars, Params, Equation
export ForwardEulerTimeStepper, ETDRK4TimeStepper
export RK4TimeStepper, AB3TimeStepper

export set_q!, updatevars!, calc_NL!, calc_NL, stepforward!

Grid = TwoDGrid

# Params type: this type defines physical parameters of the equation.
type Params <: AbstractParams
  nu::Float64                     # Vorticity viscosity
  nun::Int                        # Vorticity hyperviscous order
end

# Equation type: this type defines the problems' linear and nonlinear parts
type Equation <: AbstractEquation
  LC::Array{Complex{Float64}, 2}  # Element-wise coeff of the eqn's linear part
  calcNL!::Function               # Function to calculate eqn's nonlinear part
end

function Equation(p::Params, g::Grid)
  # Function calcNL! is defined below.
  LC = -p.nu * g.KKrsq.^(0.5*p.nun)
  Equation(LC, calcNL!)
end

# Vars type:
type Vars <: AbstractVars

  t::Float64
  sol::Array{Complex128, 2}

  # Auxiliary vars
  q::Array{Float64, 2}
  U::Array{Float64, 2}
  V::Array{Float64, 2}
  Uq::Array{Float64, 2}
  Vq::Array{Float64, 2}
  psi::Array{Float64, 2}

  # Solution
  qh::Array{Complex128, 2}
  Uh::Array{Complex128, 2}
  Vh::Array{Complex128, 2}
  Uqh::Array{Complex128, 2}
  Vqh::Array{Complex128, 2}
  psih::Array{Complex128, 2}

end

function Vars(g::Grid)
  # Initialize with t=0
  t = 0.0
  sol  = zeros(Complex128, g.nkr, g.nl)

  # Vorticity auxiliary vars
  q    = zeros(Float64, g.nx, g.ny)
  U    = zeros(Float64, g.nx, g.ny)
  V    = zeros(Float64, g.nx, g.ny)
  Uq   = zeros(Float64, g.nx, g.ny)
  Vq   = zeros(Float64, g.nx, g.ny)
  psi  = zeros(Float64, g.nx, g.ny)

  qh   = zeros(Complex128, g.nkr, g.nl)
  Uh   = zeros(Complex128, g.nkr, g.nl)
  Vh   = zeros(Complex128, g.nkr, g.nl)
  Uqh  = zeros(Complex128, g.nkr, g.nl)
  Vqh  = zeros(Complex128, g.nkr, g.nl)
  psih = zeros(Complex128, g.nkr, g.nl)

  # Random initial condition
  sol = exp.( 2.0*pi*im*rand(g.nkr, g.nl) )

  return Vars(t, sol, q, U, V, Uq, Vq, psi, qh, Uh, Vh, Uqh, Vqh, psih)
end

function build_problem(nx::Int, Lx::Float64, nu::Float64, nun::Int)
  g  = Grid(nx, Lx)
  p  = Params(nu, nun, g)
  v  = Vars(g)
  eq = Equation(p, g)
  return eq, v, p, g
end


# -----------------------------------------------------------------------------
# Solver ----------------------------------------------------------------------
# -----------------------------------------------------------------------------
function calcNL!(NL::Array{Complex{Float64}, 2}, sol::Array{Complex{Float64}, 2},
  t::Float64, v::Vars, p::Params, g::Grid)

  # ON NAVID'S LAPTOP A_mul_B! messes up!!
  #v.q = irfft(sol[:, :, 1], g.nx)
  #v.U = irfft(v.Uh, g.nx)
  #v.V = irfft(v.Vh, g.nx)
  #v.Uqh = rfft(v.Uq)
  #v.Vqh = rfft(v.Vq)

  # Solution key:
  # sol = qh

  A_mul_B!( v.q, g.irfftplan, sol )

  v.Uh .=    im .* g.Lr .* g.invKKrsq .* sol
  v.Vh .= (-im) .* g.Kr .* g.invKKrsq .* sol
 
  A_mul_B!( v.U, g.irfftplan, v.Uh )
  A_mul_B!( v.V, g.irfftplan, v.Vh )

  v.Uq .= v.U.*v.q
  v.Vq .= v.V.*v.q

  A_mul_B!( v.Uqh, g.rfftplan, v.Uq )
  A_mul_B!( v.Vqh, g.rfftplan, v.Vq )

  NL .= (-im) .* g.Kr.*v.Uqh .- im .* g.Lr.*v.Vqh

end


function calcNL(NL::Array{Complex{Float64}, 2}, sol::Array{Complex{Float64}, 2},
  t::Float64, v::Vars, p::Params, g::Grid)

  # ON NAVID'S LAPTOP A_mul_B! messes up!!
  #v.q = irfft(sol[:, :, 1], g.nx)
  #v.U = irfft(v.Uh, g.nx)
  #v.V = irfft(v.Vh, g.nx)
  #v.Uqh = rfft(v.Uq)
  #v.Vqh = rfft(v.Vq)

  # Solution key:
  # sol = qh

  A_mul_B!( v.q, g.irfftplan, sol )

  v.Uh .=    im .* g.Lr .* g.invKKrsq .* sol
  v.Vh .= (-im) .* g.Kr .* g.invKKrsq .* sol
 
  A_mul_B!( v.U, g.irfftplan, v.Uh )
  A_mul_B!( v.V, g.irfftplan, v.Vh )

  v.Uq .= v.U.*v.q
  v.Vq .= v.V.*v.q

  A_mul_B!( v.Uqh, g.rfftplan, v.Uq )
  A_mul_B!( v.Vqh, g.rfftplan, v.Vq )

  NL .= (-im) .* g.Kr.*v.Uqh .- im .* g.Lr.*v.Vqh
end




function calcNL!(NL::Array{Complex{Float64}, 3}, sol::Array{Complex128, 3},
  t::Float64, v::Vars, p::Params, g::Grid)

  # ON NAVID'S LAPTOP A_mul_B! messes up!!
  #v.q = irfft(sol[:, :, 1], g.nx)
  #v.U = irfft(v.Uh, g.nx)
  #v.V = irfft(v.Vh, g.nx)
  #v.Uqh = rfft(v.Uq)
  #v.Vqh = rfft(v.Vq)

  # Solution key:
  # sol[:, ;, 1] : qh

  A_mul_B!( v.q, g.irfftplan, sol[:, :, 1] )

  v.Uh .=    im .* g.Lr .* g.invKKrsq .* sol[:, :, 1]
  v.Vh .= (-im) .* g.Kr .* g.invKKrsq .* sol[:, :, 1]
 
  A_mul_B!( v.U, g.irfftplan, v.Uh )
  A_mul_B!( v.V, g.irfftplan, v.Vh )

  v.Uq .= v.U .* v.q
  v.Vq .= v.V .* v.q

  A_mul_B!( v.Uqh, g.rfftplan, v.Uq )
  A_mul_B!( v.Vqh, g.rfftplan, v.Vq )

  NL[:, :, 1] .= (-im) .* g.Kr.*v.Uqh  .-  im .* g.Lr.*v.Vqh

end

# -----------------------------------------------------------------------------
# Helper functions ------------------------------------------------------------
# -----------------------------------------------------------------------------
function updatevars!(v::Vars, g::Grid)

  v.qh .= v.sol

  A_mul_B!( v.q, g.irfftplan, v.qh )

  v.psih .= .- v.qh .* g.invKKrsq

  v.Uh .=    im .* g.Lr .* g.invKKrsq .* v.qh
  v.Vh .= (-im) .* g.Kr .* g.invKKrsq .* v.qh
 
  A_mul_B!( v.U, g.irfftplan, v.Uh )
  A_mul_B!( v.V, g.irfftplan, v.Vh )

  v.Uq .= v.U .* v.q
  v.Vq .= v.V .* v.q

  A_mul_B!( v.Uqh, g.rfftplan, v.Uq )
  A_mul_B!( v.Vqh, g.rfftplan, v.Vq )

end


function updatevars!(v::Vars, p::Params, g::Grid)

  v.qh .= v.sol

  A_mul_B!( v.q, g.irfftplan, v.qh )

  v.psih .= .- v.qh .* g.invKKrsq

  v.Uh .=    im .* g.Lr .* g.invKKrsq .* v.qh
  v.Vh .= (-im) .* g.Kr .* g.invKKrsq .* v.qh
 
  A_mul_B!( v.U, g.irfftplan, v.Uh )
  A_mul_B!( v.V, g.irfftplan, v.Vh )

  v.Uq .= v.U .* v.q
  v.Vq .= v.V .* v.q

  A_mul_B!( v.Uqh, g.rfftplan, v.Uq )
  A_mul_B!( v.Vqh, g.rfftplan, v.Vq )

end

# This function exists only to test the speed of fused vs hand-coded loops.
function updatevars!(v::Vars, p::Params, g::Grid, withloops::Bool)

  for j = 1:g.nl, i = 1:g.nkr 
    v.qh[i, j] = v.sol[i, j]
  end

  A_mul_B!( v.q, g.irfftplan, v.qh )

  for j = 1:g.nl, i = 1:g.nkr 
    v.psih[i, j] = -v.qh[i, j] * g.invKKrsq[i, j]
  end

  for j = 1:g.nl, i = 1:g.nkr 
    v.Uh[i, j] =  im*g.Lr[i, j]*g.invKKrsq[i, j]*v.qh[i, j]
    v.Vh[i, j] = -im*g.Lr[i, j]*g.invKKrsq[i, j]*v.qh[i, j]
  end
 
  A_mul_B!( v.U, g.irfftplan, v.Uh )
  A_mul_B!( v.V, g.irfftplan, v.Vh )

  for j = 1:g.ny, i = 1:g.nx
    v.Uq[i, j] = v.U[i, j]*v.q[i, j]
    v.Vq[i, j] = v.V[i, j]*v.q[i, j]
  end

  A_mul_B!( v.Uqh, g.rfftplan, v.Uq )
  A_mul_B!( v.Vqh, g.rfftplan, v.Vq )

end


function set_q!(v::Vars, g::Grid, q::Array{Float64, 2})
  # Set vorticity
  A_mul_B!( v.sol, g.rfftplan, q )
  updatevars!(v, g)
end


# Include solver-related functions from solvers.jl.
#include("../solvers.jl")

end


# E N D    M O D U L E    T W O D T U R B --------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
