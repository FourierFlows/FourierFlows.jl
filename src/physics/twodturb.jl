__precompile__()


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# T W O D T U R B >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
module TwoDTurb

using FourierFlows

export Grid,
       Params,
       Vars,
       Equation

export set_q!, updatevars!

# 2D grids for Two-D turbulence.
Grid = TwoDGrid




# P A R A M S ----------------------------------------------------------------- 
type Params <: AbstractParams
  nu::Float64                     # Vorticity viscosity
  nun::Int                        # Vorticity hyperviscous order
end




# E Q U A T I O N S ----------------------------------------------------------- 
type Equation <: AbstractEquation
  LC::Array{Complex{Float64}, 2}  # Element-wise coeff of the eqn's linear part
  calcNL!::Function               # Function to calculate eqn's nonlinear part
end

function Equation(p::Params, g::TwoDGrid)
  # Function calcNL! is defined below.
  LC = -p.nu * g.KKrsq.^(0.5*p.nun)
  Equation(LC, calcNL!)
end




# V A R S --------------------------------------------------------------------- 
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

function Vars(g::TwoDGrid)
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




# S O L V E R S ---------------------------------------------------------------

function calcNL!(NL::Array{Complex{Float64}, 2}, sol::Array{Complex{Float64}, 2},
  t::Float64, v::Vars, p::Params, g::TwoDGrid)

  # This copy is necessary because calling A_mul_B(v.q, g.irfftplan, sol) 
  # a few lines below destroys sol when using Julia's FFTW.
  v.qh .= sol

  A_mul_B!(v.q, g.irfftplan, sol)

  v.Uh .=    im .* g.Lr .* g.invKKrsq .* v.qh
  v.Vh .= (-im) .* g.Kr .* g.invKKrsq .* v.qh
 
  A_mul_B!(v.U, g.irfftplan, v.Uh)
  A_mul_B!(v.V, g.irfftplan, v.Vh)

  v.Uq .= v.U.*v.q
  v.Vq .= v.V.*v.q

  A_mul_B!(v.Uqh, g.rfftplan, v.Uq)
  A_mul_B!(v.Vqh, g.rfftplan, v.Vq)

  NL .= (-im) .* g.Kr.*v.Uqh .- im .* g.Lr.*v.Vqh

end




# H E L P E R   F U N C T I O N S --------------------------------------------- 
function updatevars!(v::Vars, g::TwoDGrid)

  v.qh .= v.sol

  # We don't use A_mul_B here because irfft destroys its input.
  # A_mul_B!(v.q, g.irfftplan, v.qh)
  v.q = irfft(v.qh, g.nx)

  v.psih .= .- v.qh .* g.invKKrsq

  v.Uh .=    im .* g.Lr .* g.invKKrsq .* v.qh
  v.Vh .= (-im) .* g.Kr .* g.invKKrsq .* v.qh
 
  # We don't use A_mul_B here because irfft destroys its input.
  #A_mul_B!(v.U, g.irfftplan, v.Uh)
  #A_mul_B!(v.V, g.irfftplan, v.Vh)
  v.U = irfft(v.Uh, g.nx)
  v.V = irfft(v.Vh, g.nx)

  v.Uq .= v.U .* v.q
  v.Vq .= v.V .* v.q

  A_mul_B!(v.Uqh, g.rfftplan, v.Uq)
  A_mul_B!(v.Vqh, g.rfftplan, v.Vq)

end




# This function exists only to test the speed of fused vs hand-coded loops.
function updatevars!(v::Vars, p::Params, g::TwoDGrid, withloops::Bool)

  for j = 1:g.nl, i = 1:g.nkr 
    v.qh[i, j] = v.sol[i, j]
  end

  # We don't use A_mul_B here because irfft destroys its input.
  #A_mul_B!(v.q, g.irfftplan, v.qh)
  v.q = irfft(v.qh)

  for j = 1:g.nl, i = 1:g.nkr 
    v.psih[i, j] = -v.qh[i, j] * g.invKKrsq[i, j]
  end

  for j = 1:g.nl, i = 1:g.nkr 
    v.Uh[i, j] =  im*g.Lr[i, j]*g.invKKrsq[i, j]*v.qh[i, j]
    v.Vh[i, j] = -im*g.Lr[i, j]*g.invKKrsq[i, j]*v.qh[i, j]
  end
 
  # We don't use A_mul_B here because irfft destroys its input.
  #A_mul_B!(v.U, g.irfftplan, v.Uh)
  #A_mul_B!(v.V, g.irfftplan, v.Vh)
  v.U = irfft(v.Uh)
  v.V = irfft(v.Vh)

  for j = 1:g.ny, i = 1:g.nx
    v.Uq[i, j] = v.U[i, j]*v.q[i, j]
    v.Vq[i, j] = v.V[i, j]*v.q[i, j]
  end

  A_mul_B!(v.Uqh, g.rfftplan, v.Uq)
  A_mul_B!(v.Vqh, g.rfftplan, v.Vq)

end




function set_q!(v::Vars, g::TwoDGrid, q::Array{Float64, 2})
  # Set vorticity
  A_mul_B!(v.sol, g.rfftplan, q)
  updatevars!(v, g)
end




# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# S E T U P S >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
module Setups

using FourierFlows.TwoDTurb

export simplenondim

function simplenondim(nx::Int; nu=1e-6, nun=4)
  # Construct a barotropic QG problem that should reproduce results obtained
  # from a bare 2D turbulence simulation when beta=0.

  Lx     = 2.0*pi                   # Domain size (meters)
   
  g  = Grid(nx, Lx)
  p  = Params(nu, nun)
  v  = Vars(g)
  eq = Equation(p, g)

  set_q!(v, g, rand(nx, nx))        # Random initial condition

  return g, p, v, eq
end

end
# E N D   S E T U P S >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>




end
# E N D   T W O D T U R B >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
