module TimeSteppers

export ForwardEulerTimeStepper, ETDRK4TimeStepper

mutable struct ForwardEulerTimeStepper
  dt::Float64
  LC::Array{Complex128, 2}
  NL::Array{Complex128, 2}
end

function ForwardEulerTimeStepper(dt::Float64, LC::Array{Complex128, 2})
  NL = zeros(LC)
  ForwardEulerTimeStepper(dt, LC, NL)
end


mutable struct ETDRK4TimeStepper
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

  # Intermediate times, arrays and solutions
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

  ETDRK4TimeStepper(dt, LC, zeta, alph, beta, gamm, expLCdt, expLCdt2,
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
    (2.0 + zc + exp.(zc).*(2.0 + zc))./zc.^3.0, 3), 3)

  gamm = dt*squeeze(mean(
    (-4.0 - 3.0*zc - zc.^2.0 + exp.(zc).*(4.0 - zc))./zc.^3.0, 3), 3)

  return zeta, alph, beta, gamm
end



end
