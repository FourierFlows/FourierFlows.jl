__precompile__()


# Types and functions associated with time-stepping Fourier doubly-periodic
# partial differential equations / flows
module TimeSteppers

export ForwardEulerTimeStepper, ETDRK4TimeStepper, RK4TimeStepper
export AB3TimeStepper

# ----------------------------------------------------------------------------- 
# Forward Euler
# ----------------------------------------------------------------------------- 
# The simplest time-stepping method in the books. Explicit and 1st-order
# accurate.

type ForwardEulerTimeStepper{dim}
  step::Int
  dt::Float64

  # Linear part of the problem
  LC::Array{Complex{Float64}, dim}

  # Pre-allocated nonlinear term
  NL::Array{Complex{Float64}, dim}
end

function ForwardEulerTimeStepper(dt::Float64, LC::Array{Complex{Float64}})
  NL = zeros(LC)
  ForwardEulerTimeStepper{ndims(LC)}(0, dt, LC, NL)
end


# ----------------------------------------------------------------------------- 
# ETDRK4
# ----------------------------------------------------------------------------- 
# The Rolls-Royce of time-stepping. Exact treatment of linear part of 
# the equation, explicit and 4th-order accurate integration of nonlinear
# parts of equation.

type ETDRK4TimeStepper
  step::Int
  dt::Float64

  # Linear part of the problem
  LC::Array{Complex128, 2}

  # ETDRK4 coefficents
  zeta::Array{Complex128, 2}
  alph::Array{Complex128, 2}
  beta::Array{Complex128, 2}
  gamm::Array{Complex128, 2}

  # Precomputed arrays: exp(LC*dt) and exp(LC*dt/2)
  expLCdt::Array{Complex128, 2}
  expLCdt2::Array{Complex128, 2}

  # Intermediate times, solutions, and nonlinear evaluations
  ti::Float64

  sol1::Array{Complex128, 2}
  sol2::Array{Complex128, 2}

  NL1::Array{Complex128, 2}
  NL2::Array{Complex128, 2}
  NL3::Array{Complex128, 2}
  NL4::Array{Complex128, 2}
end



function ETDRK4TimeStepper(dt::Float64, LC::Array{Complex128, 2})

  # dt is the time step size
  # LC is the linear coefficient of the equation, such that
  #   A_t + LC*A = NL

  expLCdt  = exp.(dt*LC)
  expLCdt2 = exp.(0.5*dt*LC)

  alph, beta, gamm, zeta = get_etd_coeffs(dt, LC)

  ti = 0.0

  # Array size needed locally
  ni, nj = size(LC)

  sol1 = zeros(LC)
  sol2 = zeros(LC)

  NL1 = zeros(LC)
  NL2 = zeros(LC)
  NL3 = zeros(LC)
  NL4 = zeros(LC)

  ETDRK4TimeStepper(0, dt, LC, zeta, alph, beta, gamm, expLCdt, expLCdt2,
    ti, sol1, sol2, NL1, NL2, NL3, NL4)
end


function get_etd_coeffs(dt::Float64, LC::Array{Complex128, 2})

  # Calculate ETDRK4 coefficients by integrating over a small circle
  # in complex space.

  # Circle parameters
  ncirc = 32
  rcirc = 1.0

  # Make circle
  circ  = Array{Complex128}(1, 1, ncirc)
  circ[1, 1, :]  = rcirc * exp.( 2.0*pi*im*( (1:ncirc)-0.5 ) / ncirc )

  # Construct intermediate vars
  zc = broadcast(+, dt*LC, circ)

  # Four coefficients: zeta, alph, beta, gamm
  zeta = dt*squeeze(mean( (exp.(0.5*zc) - 1.0)./zc, 3), 3)

  alph = dt*squeeze(mean(
    (-4.0 - zc + exp.(zc).*(4.0 - 3.0*zc + zc.^2.0))./zc.^3.0, 3), 3)

  beta = dt*squeeze(mean(
    (2.0 + zc + exp.(zc).*(-2.0 + zc))./zc.^3.0, 3), 3)

  gamm = dt*squeeze(mean(
    (-4.0 - 3.0*zc - zc.^2.0 + exp.(zc).*(4.0 - zc))./zc.^3.0, 3), 3)

  return zeta, alph, beta, gamm
end


# ----------------------------------------------------------------------------- 
# RK4
# ----------------------------------------------------------------------------- 
# RK4 is the classical explicit 4th-order Runge-Kutta time-stepping
# method. It uses a series of substeps/estimators to achieve 4th-order
# accuracy over each individual time-step, at the cost of requiring
# relatively more evaluations of the nonlinear right hand side.
# It is described, among other places, in Bewley's Numerical
# Renaissance.

type RK4TimeStepper
  step::Int
  dt::Float64

  # Linear part of the problem
  LC::Array{Complex128, 2}

  # Intermediate times, solutions, and nonlinear evaluations
  ti::Float64

  sol1::Array{Complex128, 2}

  RHS1::Array{Complex128, 2}
  RHS2::Array{Complex128, 2}
  RHS3::Array{Complex128, 2}
  RHS4::Array{Complex128, 2}
end

function RK4TimeStepper(dt::Float64, LC::Array{Complex128, 2})
  ti = 0.0
  sol1 = zeros(LC)
  RHS1  = zeros(LC)
  RHS2  = zeros(LC)
  RHS3  = zeros(LC)
  RHS4  = zeros(LC)
  RK4TimeStepper(0, dt, LC, ti, sol1, RHS1, RHS2, RHS3, RHS4)
end






# ----------------------------------------------------------------------------- 
# AB3
# ----------------------------------------------------------------------------- 
# 3rd order Adams-Bashforth time stepping is an explicit scheme that uses
# solutions from two previous time-steps to achieve 3rd order accuracy.

type AB3TimeStepper
  step::Int
  dt::Float64

  # Linear part of the problem
  LC::Array{Complex128, 2}

  # Intermediate times, solutions, and nonlinear evaluations
  RHS::Array{Complex128, 2}
  RHSm1::Array{Complex128, 2}
  RHSm2::Array{Complex128, 2}
end

function AB3TimeStepper(dt::Float64, LC::Array{Complex128, 2})
  RHS   = zeros(LC)
  RHSm1 = zeros(LC)
  RHSm2 = zeros(LC)
  AB3TimeStepper(0, dt, LC, RHS, RHSm1, RHSm2)
end



end
