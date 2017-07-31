__precompile__()




export ForwardEulerTimeStepper,
       AB3TimeStepper,
       RK4TimeStepper,
       ETDRK4TimeStepper

export stepforward!,
       step_nsteps!




# Looping stepforward functions -----------------------------------------------
function stepforward!(v::AbstractVars, ts::AbstractTimeStepper,
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid, nsteps::Int)

  for i = 1:nsteps
    stepforward!(v, ts, eq, p, g)
  end
end

function stepforward!(v::AbstractVars, ts::AbstractTimeStepper,
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid; nsteps=1)
  for i = 1:nsteps
    stepforward!(v, ts, eq, p, g)
  end
end

# Let's depreciate this form!
function stepforward!(v::AbstractVars, nsteps::Int, ts::AbstractTimeStepper,
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid)

  for i = 1:nsteps
    stepforward!(v, ts, eq, p, g)
  end
end




# Forward Euler ---------------------------------------------------------------
# The simplest time-stepping method in the books. Explicit and 1st-order
# accurate.

type ForwardEulerTimeStepper{dim} <: AbstractTimeStepper
  step::Int
  dt::Float64
  NL::Array{Complex{Float64}, dim}    # Nonlinear term
end

function ForwardEulerTimeStepper(dt::Float64, v::AbstractVars)
  NL = zeros(v.sol)
  ForwardEulerTimeStepper{ndims(NL)}(0, dt, NL)
end

function ForwardEulerTimeStepper(dt::Float64, LC::Array{Complex{Float64}, 2})
  NL = zeros(LC)
  ForwardEulerTimeStepper{ndims(LC)}(0, dt, NL)
end




function stepforward!(v::AbstractVars, ts::ForwardEulerTimeStepper,
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid)

  eq.calcNL!(ts.NL, v.sol, v.t, v, p, g)

  v.sol .+= ts.dt .* (ts.NL .+ eq.LC.*v.sol)
  v.t    += ts.dt

  ts.step += 1
end








# ETDRK4 ----------------------------------------------------------------------
# The Rolls-Royce of time-stepping. Exact treatment of linear part of
# the equation, explicit and 4th-order accurate integration of nonlinear
# parts of equation.

type ETDRK4TimeStepper <: AbstractTimeStepper
  step::Int
  dt::Float64

  LC::Array{Complex{Float64}, 2}          # Linear coefficient

  # ETDRK4 coefficents
  zeta::Array{Complex{Float64}, 2}
  alph::Array{Complex{Float64}, 2}
  beta::Array{Complex{Float64}, 2}
  gamm::Array{Complex{Float64}, 2}
  expLCdt::Array{Complex{Float64}, 2}     # Precomputed exp(LC*dt)
  expLCdt2::Array{Complex{Float64}, 2}    # Precomputed exp(LC*dt/2)

  # Intermediate times, solutions, and nonlinear evaluations
  ti::Float64
  sol1::Array{Complex{Float64}, 2}
  sol2::Array{Complex{Float64}, 2}
  NL1::Array{Complex{Float64}, 2}
  NL2::Array{Complex{Float64}, 2}
  NL3::Array{Complex{Float64}, 2}
  NL4::Array{Complex{Float64}, 2}
end

function ETDRK4TimeStepper(dt::Float64, LC::Array{Complex{Float64}, 2})
  expLCdt  = exp.(dt.*LC)
  expLCdt2 = exp.(0.5.*dt.*LC)

  zeta, alph, beta, gamm = get_etd_coeffs(dt, LC)

  ti = 0.0

  sol1 = zeros(LC)
  sol2 = zeros(LC)
  NL1  = zeros(LC)
  NL2  = zeros(LC)
  NL3  = zeros(LC)
  NL4  = zeros(LC)

  ETDRK4TimeStepper(0, dt, LC, zeta, alph, beta, gamm, expLCdt, expLCdt2,
    ti, sol1, sol2, NL1, NL2, NL3, NL4)
end


function get_etd_coeffs(dt::Float64, LC::Array{Complex{Float64}, 2})

  # Calculate ETDRK4 coefficients by integrating over a small circle
  # in complex space.

  # Circle parameters
  ncirc = 32
  rcirc = 1.0

  # Make circle
  circ  = Array{Complex{Float64}}(1, 1, ncirc)
  circ[1, 1, :]  = rcirc * exp.( 2.0*pi*im*(0.5:1.0:(ncirc-0.5))/ncirc )

  # Construct intermediate vars
  zc = broadcast(+, dt.*LC, circ)

  # Four coefficients: zeta, alph, beta, gamm
  zeta = dt.*squeeze(mean( (exp.(0.5.*zc)-1.0)./zc, 3), 3)

  alph = dt.*squeeze(mean(
    ( (-4.0) .- zc .+ exp.(zc).*(4.0.-3.0.*zc.+zc.^2.0))./zc.^3.0,      3), 3)
  beta = dt.*squeeze(mean(
    (   2.0  .+ zc .+ exp.(zc).*((-2.0).+zc))./zc.^3.0,                 3), 3)
  gamm = dt.*squeeze(mean(
    ( (-4.0) .- 3.0.*zc .- zc.^2.0 .+ exp.(zc).*(  4.0 .-zc))./zc.^3.0, 3), 3)

  return zeta, alph, beta, gamm
end




function stepforward!(v::AbstractVars, ts::ETDRK4TimeStepper,
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid)

  # Substep 1
  eq.calcNL!(ts.NL1, v.sol, v.t, v, p, g)
  ts.sol1 .= ts.expLCdt2.*v.sol .+ ts.zeta.*ts.NL1

  # Substep 2
  ts.ti = v.t + 0.5*ts.dt
  eq.calcNL!(ts.NL2, ts.sol1, ts.ti, v, p, g)
  ts.sol2 .= ts.expLCdt2.*v.sol .+ ts.zeta.*ts.NL2

  # Substep 3
  eq.calcNL!(ts.NL3, ts.sol2, ts.ti, v, p, g)
  ts.sol2 .= ts.expLCdt2.*ts.sol1 .+ ts.zeta.*(2.0.*ts.NL3 .- ts.NL1)

  # Substep 4
  ts.ti = v.t + ts.dt
  eq.calcNL!(ts.NL4, ts.sol2, ts.ti, v, p, g)

  # Update
  v.sol .= (ts.expLCdt.*v.sol .+      ts.alph .* ts.NL1
                              .+ 2.0.*ts.beta .* (ts.NL2 .+ ts.NL3)
                              .+      ts.gamm .* ts.NL4 )
  v.t   += ts.dt
  ts.step += 1

end








# RK4 -------------------------------------------------------------------------
# RK4 is the classical explicit 4th-order Runge-Kutta time-stepping
# method. It uses a series of substeps/estimators to achieve 4th-order
# accuracy over each individual time-step, at the cost of requiring
# relatively more evaluations of the nonlinear right hand side.
# It is described, among other places, in Bewley's Numerical
# Renaissance.

type RK4TimeStepper <: AbstractTimeStepper
  step::Int
  dt::Float64

  # Intermediate times, solutions, and nonlinear evaluations
  ti::Float64
  sol1::Array{Complex{Float64}, 2}
  RHS1::Array{Complex{Float64}, 2}
  RHS2::Array{Complex{Float64}, 2}
  RHS3::Array{Complex{Float64}, 2}
  RHS4::Array{Complex{Float64}, 2}
end

function RK4TimeStepper(dt::Float64, v::AbstractVars)
  ti = 0.0
  sol1 = zeros(v.sol)
  RHS1 = zeros(v.sol)
  RHS2 = zeros(v.sol)
  RHS3 = zeros(v.sol)
  RHS4 = zeros(v.sol)
  RK4TimeStepper(0, dt, ti, sol1, RHS1, RHS2, RHS3, RHS4)
end

function RK4TimeStepper(dt::Float64, LC::Array{Complex{Float64}, 2})
  ti = 0.0
  sol1 = zeros(LC)
  RHS1 = zeros(LC)
  RHS2 = zeros(LC)
  RHS3 = zeros(LC)
  RHS4 = zeros(LC)
  RK4TimeStepper(0, dt, ti, sol1, RHS1, RHS2, RHS3, RHS4)
end




function stepforward!(v::AbstractVars, ts::RK4TimeStepper,
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid)

  eq.calcNL!(ts.RHS1, v.sol, v.t, v, p, g)
  ts.RHS1 .+= eq.LC.*v.sol

  # Substep 1
  ts.ti = v.t + 0.5*ts.dt
  ts.sol1 .= v.sol .+ (0.5*ts.dt).*ts.RHS1
  eq.calcNL!(ts.RHS2, ts.sol1, v.t, v, p, g)
  ts.RHS2 .+= eq.LC.*ts.sol1

  # Substep 2
  ts.sol1 .= v.sol .+ (0.5*ts.dt).*ts.RHS2
  eq.calcNL!(ts.RHS3, ts.sol1, v.t, v, p, g)
  ts.RHS3 .+= eq.LC.*ts.sol1

  # Substep 3
  ts.ti = v.t + ts.dt
  ts.sol1 .= v.sol .+ ts.dt.*ts.RHS3
  eq.calcNL!(ts.RHS4, ts.sol1, v.t, v, p, g)
  ts.RHS4 .+= eq.LC.*ts.sol1

  # Substep 4 and final step
  v.sol .+= ts.dt.*(
       (1.0/6.0).*ts.RHS1 .+ (1.0/3.0).*ts.RHS2
    .+ (1.0/3.0).*ts.RHS3 .+ (1.0/6.0).*ts.RHS4 )

  v.t += ts.dt
  ts.step += 1
end








# AB3 -------------------------------------------------------------------------
# 3rd order Adams-Bashforth time stepping is an explicit scheme that uses
# solutions from two previous time-steps to achieve 3rd order accuracy.

type AB3TimeStepper <: AbstractTimeStepper
  step::Int
  dt::Float64

  # Intermediate times, solutions, and nonlinear evaluations
  RHS::Array{Complex{Float64}, 2}
  RHSm1::Array{Complex{Float64}, 2}
  RHSm2::Array{Complex{Float64}, 2}
end

function AB3TimeStepper(dt::Float64, v::AbstractVars)
  RHS   = zeros(v.sol)
  RHSm1 = zeros(v.sol)
  RHSm2 = zeros(v.sol)
  AB3TimeStepper(0, dt, RHS, RHSm1, RHSm2)
end

function AB3TimeStepper(dt::Float64, LC::Array{Complex{Float64}, 2})
  RHS   = zeros(LC)
  RHSm1 = zeros(LC)
  RHSm2 = zeros(LC)
  AB3TimeStepper(0, dt, RHS, RHSm1, RHSm2)
end




function stepforward!(v::AbstractVars, ts::AB3TimeStepper,
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid)

  eq.calcNL!(ts.RHS, v.sol, v.t, v, p, g)
  ts.RHS  .+= eq.LC.*v.sol

  v.sol .+= ts.dt .* (
    (23.0/12.0).*ts.RHS .- (16.0/12.0).*ts.RHSm1 .+ (5.0/12.0).*ts.RHSm2 )

  v.t += ts.dt
  ts.step += 1

  ts.RHSm2 .= ts.RHSm1
  ts.RHSm1 .= ts.RHS

end




function stepforward!(v::AbstractVars, nsteps::Int,
  ts::AB3TimeStepper, eq::AbstractEquation,
  p::AbstractParams, g::AbstractGrid;
  initstepper=false)

  # AB3 initialization
  initstepper ? initsteps=0 : initsteps=ts.step

  while initsteps < 2
    # Take forward Euler steps to initialize AB3
    eq.calcNL!(ts.RHS, v.sol, v.t, v, p, g)
    ts.RHS .+= eq.LC.*v.sol

    # Update
    ts.RHSm2 .= ts.RHSm1
    ts.RHSm1 .= ts.RHS

    v.sol     .+= ts.dt.*ts.RHS
    v.t        += ts.dt
    ts.step    += 1
    initsteps  += 1

    # Account for steps taken
    nsteps -= 1
  end

  # Loop
  for step = 1:nsteps
    stepforward!(v, ts, eq, p, g)
  end

end
