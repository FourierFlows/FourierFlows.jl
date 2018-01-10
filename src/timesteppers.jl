__precompile__()


export ForwardEulerTimeStepper, FilteredForwardEulerTimeStepper,
       AB3TimeStepper,
       RK4TimeStepper,
       ETDRK4TimeStepper, FilteredETDRK4TimeStepper

export stepforward!




# Looping stepforward function ------------------------------------------------
function stepforward!(v::AbstractVars, ts::AbstractTimeStepper,
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid; nsteps=1)
  for step = 1:nsteps
    stepforward!(v, ts, eq, p, g)
  end
  nothing
end

function stepforward!(prob::Problem; nsteps=1)
  for step = 1:nsteps
    stepforward!(prob.vars, prob.ts, prob.eqn, prob.params, prob.grid)
    prob.t = prob.vars.t
    prob.step = prob.ts.step
  end
  nothing
end

function stepforward!(prob::Problem, diags::AbstractArray; nsteps=1)
  # Initialize diagnostics for speed
  for diag in diags
    newnum = ceil(Int, (diag.count+nsteps)/diag.freq)
    if newnum > diag.num
      warn("Resizing diags before stepping forward...")
      resize!(diag, newnum)
    end
  end

  for step = 1:nsteps
    stepforward!(prob.vars, prob.ts, prob.eqn, prob.params, prob.grid)
    prob.t = prob.vars.t
    prob.step = prob.ts.step
    for diag in diags
      if (prob.step+1) % diag.freq == 0.0
        increment!(diag)
      end
    end

  end
  nothing
end

function stepforward!(prob::Problem, diag::AbstractDiagnostic; nsteps=1)
  stepforward!(prob, [diag]; nsteps=nsteps)
end



# Utilities -------------------------------------------------------------------
"""
Calculate ETDRK4 coefficients by integrating over a small circle
in complex space.
"""
function get_etd_coeffs(dt::Float64, LC::Array{Complex{Float64}, 2};
  ncirc=32, rcirc=1.0)

  # Make circle
  circ  = Array{Complex{Float64}}(1, 1, ncirc)
  circ[1, 1, :]  = rcirc * exp.( 2.0*pi*im*(0.5:1.0:(ncirc-0.5))/ncirc )

  # Construct intermediate vars
  zc = broadcast(+, dt.*LC, circ)

  # Four coefficients: ζ, α, β, Γ
  ζ = dt.*squeeze(mean( (exp.(0.5.*zc)-1.0)./zc, 3), 3)

  α = dt.*squeeze(mean(
    ( (-4.0) .- zc .+ exp.(zc).*(4.0.-3.0.*zc.+zc.^2.0))./zc.^3.0,      3), 3)
  β = dt.*squeeze(mean(
    (   2.0  .+ zc .+ exp.(zc).*((-2.0).+zc))./zc.^3.0,                 3), 3)
  Γ = dt.*squeeze(mean(
    ( (-4.0) .- 3.0.*zc .- zc.^2.0 .+ exp.(zc).*(  4.0 .-zc))./zc.^3.0, 3), 3)

  return ζ, α, β, Γ
end


"""
Calculate ETDRK4 coefficients by integrating over a small circle
in complex space.
"""
function get_etd_coeffs(dt::Float64, LC::Array{Complex{Float64}, 3};
  ncirc=32, rcirc=1.0)

  # Make circle
  circ  = Array{Complex{Float64}}(1, 1, 1, ncirc)
  circ[1, 1, 1, :]  = rcirc * exp.( 2.0*pi*im*(0.5:1.0:(ncirc-0.5))/ncirc )

  # Construct intermediate vars
  zc = broadcast(+, dt.*LC, circ)

  # Four coefficients: ζ, α, β, Γ
  ζ = dt.*squeeze(mean( (exp.(0.5.*zc)-1.0)./zc, 4), 4)

  α = dt.*squeeze(mean(
    ( (-4.0) .- zc .+ exp.(zc).*(4.0.-3.0.*zc.+zc.^2.0))./zc.^3.0,      4), 4)
  β = dt.*squeeze(mean(
    (   2.0  .+ zc .+ exp.(zc).*((-2.0).+zc))./zc.^3.0,                 4), 4)
  Γ = dt.*squeeze(mean(
    ( (-4.0) .- 3.0.*zc .- zc.^2.0 .+ exp.(zc).*(  4.0 .-zc))./zc.^3.0, 4), 4)

  return ζ, α, β, Γ
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

function ForwardEulerTimeStepper(dt::Float64, LC::AbstractArray)
  NL = zeros(LC)
  ForwardEulerTimeStepper{ndims(LC)}(0, dt, NL)
end




function stepforward!(v::AbstractVars, ts::ForwardEulerTimeStepper,
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid)

  eq.calcNL!(ts.NL, v.sol, v.t, v, p, g)

  @. v.sol += ts.dt*(ts.NL + eq.LC*v.sol)
  v.t += ts.dt

  ts.step += 1
end








# Filtered Forward Euler ------------------------------------------------------
# The simplest time-stepping method in the books. Explicit and 1st-order
# accurate.

type FilteredForwardEulerTimeStepper{dim} <: AbstractTimeStepper
  step::Int
  dt::Float64
  NL::Array{Complex{Float64},dim}   # Nonlinear term
  filter::Array{Float64,dim}        # Filter for solution
end

function FilteredForwardEulerTimeStepper(dt::Float64, g::AbstractGrid,
  sol::AbstractArray; filterorder=4.0, innerfilterK=0.65, outerfilterK=0.95)
  NL = zeros(sol)

  if size(sol)[1] == g.nkr
    realvars = true
  else
    realvars = false
  end

  filter = makefilter(g; order=filterorder, innerK=innerfilterK,
    outerK=outerfilterK, realvars=realvars)

  # Broadcast to correct size
  filter = ones(real.(sol)) .* filter

  FilteredForwardEulerTimeStepper{ndims(NL)}(0, dt, NL, filter)
end

function FilteredForwardEulerTimeStepper(dt::Float64, g::AbstractGrid,
  v::AbstractVars; filterorder=4.0, innerfilterK=0.65, outerfilterK=0.95)
  FilteredForwardEulerTimeStepper(dt, g, v.sol; filterorder=filterorder,
    innerfilterK=innerfilterK, outerfilterK=outerfilterK)
end



function stepforward!(v::AbstractVars, ts::FilteredForwardEulerTimeStepper,
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid)

  eq.calcNL!(ts.NL, v.sol, v.t, v, p, g)

  @. v.sol = ts.filter * ( v.sol + ts.dt*(ts.NL + eq.LC.*v.sol) )
  v.t += ts.dt
  ts.step += 1
end








# ETDRK4 ----------------------------------------------------------------------
# The Rolls-Royce of time-stepping. Exact treatment of linear part of
# the equation, explicit and 4th-order accurate integration of nonlinear
# parts of equation.

type ETDRK4TimeStepper{dim} <: AbstractTimeStepper
  step::Int
  dt::Float64

  LC::Array{Complex{Float64}, dim}          # Linear coefficient

  # ETDRK4 coefficents
  ζ::Array{Complex{Float64}, dim}
  α::Array{Complex{Float64}, dim}
  β::Array{Complex{Float64}, dim}
  Γ::Array{Complex{Float64}, dim}
  expLCdt::Array{Complex{Float64}, dim}     # Precomputed exp(LC*dt)
  expLCdt2::Array{Complex{Float64}, dim}    # Precomputed exp(LC*dt/2)

  # Intermediate times, solutions, and nonlinear evaluations
  ti::Float64
  sol1::Array{Complex{Float64}, dim}
  sol2::Array{Complex{Float64}, dim}
  NL1::Array{Complex{Float64}, dim}
  NL2::Array{Complex{Float64}, dim}
  NL3::Array{Complex{Float64}, dim}
  NL4::Array{Complex{Float64}, dim}
end

function ETDRK4TimeStepper(dt::Float64, LC::AbstractArray)
  expLCdt  = exp.(dt*LC)
  expLCdt2 = exp.(0.5*dt*LC)

  ζ, α, β, Γ = get_etd_coeffs(dt, LC)

  ti = 0.0

  sol1 = zeros(LC)
  sol2 = zeros(LC)
  NL1  = zeros(LC)
  NL2  = zeros(LC)
  NL3  = zeros(LC)
  NL4  = zeros(LC)

  ETDRK4TimeStepper{ndims(LC)}(0, dt, LC, ζ, α, β, Γ, expLCdt,
    expLCdt2, ti, sol1, sol2, NL1, NL2, NL3, NL4)
end




# Filtered ETDRK4 --------------------------------------------------------------
# The Rolls-Royce of time-stepping. Exact treatment of linear part of
# the equation, explicit and 4th-order accurate integration of nonlinear
# parts of equation.

type FilteredETDRK4TimeStepper{dim} <: AbstractTimeStepper
  step::Int
  dt::Float64

  LC::Array{Complex{Float64}, dim}          # Linear coefficient

  # ETDRK4 coefficents
  ζ::Array{Complex{Float64}, dim}
  α::Array{Complex{Float64}, dim}
  β::Array{Complex{Float64}, dim}
  Γ::Array{Complex{Float64}, dim}
  expLCdt::Array{Complex{Float64}, dim}     # Precomputed exp(LC*dt)
  expLCdt2::Array{Complex{Float64}, dim}    # Precomputed exp(LC*dt/2)

  # Intermediate times, solutions, and nonlinear evaluations
  ti::Float64
  sol1::Array{Complex{Float64}, dim}
  sol2::Array{Complex{Float64}, dim}
  NL1::Array{Complex{Float64}, dim}
  NL2::Array{Complex{Float64}, dim}
  NL3::Array{Complex{Float64}, dim}
  NL4::Array{Complex{Float64}, dim}
  filter::Array{Complex{Float64}, dim}    # Filter for solution
end


function FilteredETDRK4TimeStepper(dt::Float64, LC::AbstractArray,
    g::AbstractGrid; filterorder=4.0, innerfilterK=0.65, outerfilterK=0.95)
  expLCdt  = exp.(dt*LC)
  expLCdt2 = exp.(0.5*dt*LC)

  ζ, α, β, Γ = get_etd_coeffs(dt, LC)

  ti = 0.0

  sol1 = zeros(LC)
  sol2 = zeros(LC)
  NL1  = zeros(LC)
  NL2  = zeros(LC)
  NL3  = zeros(LC)
  NL4  = zeros(LC)

  NL = zeros(LC)

  if size(LC)[1] == g.nkr
    realvars = true
  else
    realvars = false
  end

  filter = makefilter(g; order=filterorder, innerK=innerfilterK,
    outerK=outerfilterK, realvars=realvars)

  # Broadcast to correct size
  filter = ones(real.(LC)) .* filter

  FilteredETDRK4TimeStepper{ndims(LC)}(0, dt, LC, ζ, α, β, Γ,
        expLCdt, expLCdt2, ti, sol1, sol2, NL1, NL2, NL3, NL4, filter)
end







type DualETDRK4TimeStepper{dimc, dimr} <: AbstractTimeStepper
  step::Int
  dt::Float64
  c::ETDRK4TimeStepper{dimc}
  r::ETDRK4TimeStepper{dimr}
end

function ETDRK4TimeStepper(dt::Float64, LCc::AbstractArray, LCr::AbstractArray)
  c = ETDRK4TimeStepper(dt::Float64, LCc)
  r = ETDRK4TimeStepper(dt::Float64, LCr)
  DualETDRK4TimeStepper{ndims(LCc), ndims(LCr)}(0, dt, c, r)
end






function stepforward!(v::AbstractVars, ts::ETDRK4TimeStepper,
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid)

  # Substep 1
  eq.calcNL!(ts.NL1, v.sol, v.t, v, p, g)
  @. ts.sol1 = ts.expLCdt2*v.sol + ts.ζ*ts.NL1

  # Substep 2
  ts.ti = v.t + 0.5*ts.dt
  eq.calcNL!(ts.NL2, ts.sol1, ts.ti, v, p, g)
  @. ts.sol2 = ts.expLCdt2*v.sol + ts.ζ*ts.NL2

  # Substep 3
  eq.calcNL!(ts.NL3, ts.sol2, ts.ti, v, p, g)
  @. ts.sol2 = ts.expLCdt2*ts.sol1 + ts.ζ*(2.0*ts.NL3 - ts.NL1)

  # Substep 4
  ts.ti = v.t + ts.dt
  eq.calcNL!(ts.NL4, ts.sol2, ts.ti, v, p, g)

  # Update
  @. v.sol = (ts.expLCdt.*v.sol +     ts.α * ts.NL1
                                + 2.0*ts.β * (ts.NL2 + ts.NL3)
                                +     ts.Γ * ts.NL4 )

  v.t += ts.dt
  ts.step += 1

end





function stepforward!(v::AbstractVars, ts::FilteredETDRK4TimeStepper,
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid)

  # Substep 1
  eq.calcNL!(ts.NL1, v.sol, v.t, v, p, g)
  @. ts.sol1 = ts.expLCdt2*v.sol + ts.ζ*ts.NL1

  # Substep 2
  ts.ti = v.t + 0.5*ts.dt
  eq.calcNL!(ts.NL2, ts.sol1, ts.ti, v, p, g)
  @. ts.sol2 = ts.expLCdt2*v.sol + ts.ζ*ts.NL2

  # Substep 3
  eq.calcNL!(ts.NL3, ts.sol2, ts.ti, v, p, g)
  @. ts.sol2 = ts.expLCdt2*ts.sol1 + ts.ζ*(2.0*ts.NL3 - ts.NL1)

  # Substep 4
  ts.ti = v.t + ts.dt
  eq.calcNL!(ts.NL4, ts.sol2, ts.ti, v, p, g)

  # Update
  @. v.sol = ts.filter*(ts.expLCdt*v.sol +     ts.α * ts.NL1
                                         + 2.0*ts.β * (ts.NL2 + ts.NL3)
                                         +     ts.Γ * ts.NL4 )
  v.t += ts.dt
  ts.step += 1

end





function stepforward!(v::AbstractVars, ts::DualETDRK4TimeStepper,
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid)

  # ---------------------------------------------------------------------------
  # Substep 1
  eq.calcNL!(ts.c.NL1, ts.r.NL1, v.solc, v.solr, v.t, v, p, g)
  @. ts.c.sol1 = ts.c.expLCdt2*v.solc + ts.c.ζ*ts.c.NL1
  @. ts.r.sol1 = ts.r.expLCdt2*v.solr + ts.r.ζ*ts.r.NL1

  # Substep 2
  ts.c.ti = v.t + 0.5*ts.c.dt
  eq.calcNL!(ts.c.NL2, ts.r.NL2, ts.c.sol1, ts.r.sol1, ts.c.ti, v, p, g)
  @. ts.c.sol2 = ts.c.expLCdt2*v.solc + ts.c.ζ*ts.c.NL2
  @. ts.r.sol2 = ts.r.expLCdt2*v.solr + ts.r.ζ*ts.r.NL2

  # Substep 3
  eq.calcNL!(ts.c.NL3, ts.r.NL3, ts.c.sol2, ts.r.sol2, ts.c.ti, v, p, g)
  @. ts.c.sol2 = (ts.c.expLCdt2*ts.c.sol1
    + ts.c.ζ*(2.0*ts.c.NL3 - ts.c.NL1))
  @. ts.r.sol2 = (ts.r.expLCdt2*ts.r.sol1
    + ts.r.ζ*(2.0*ts.r.NL3 - ts.r.NL1))

  # Substep 4
  ts.c.ti = v.t + ts.c.dt
  eq.calcNL!(ts.c.NL4, ts.r.NL4, ts.c.sol2, ts.r.sol2, ts.c.ti, v, p, g)

  # Update
  @. v.solc = (ts.c.expLCdt.*v.solc +     ts.c.α * ts.c.NL1
                                    + 2.0*ts.c.β * (ts.c.NL2 + ts.c.NL3)
                                    +     ts.c.Γ * ts.c.NL4 )

  @. v.solr = (ts.r.expLCdt.*v.solr +     ts.r.α * ts.r.NL1
                                    + 2.0*ts.r.β * (ts.r.NL2 + ts.r.NL3)
                                    +     ts.r.Γ * ts.r.NL4 )

  v.t += ts.dt
  ts.step += 1
  ts.c.step += 1
  ts.r.step += 1

end















# RK4 -------------------------------------------------------------------------
# RK4 is the classical explicit 4th-order Runge-Kutta time-stepping
# method. It uses a series of substeps/estimators to achieve 4th-order
# accuracy over each individual time-step, at the cost of requiring
# relatively more evaluations of the nonlinear right hand side.
# It is described, among other places, in Bewley's Numerical
# Renaissance.

type RK4TimeStepper{T,dim} <: AbstractTimeStepper
  step::Int
  dt::Float64

  # Intermediate times, solutions, and nonlinear evaluations
  ti::Float64
  sol1::Array{T, dim}
  RHS1::Array{T, dim}
  RHS2::Array{T, dim}
  RHS3::Array{T, dim}
  RHS4::Array{T, dim}
end

function RK4TimeStepper(dt::Float64, v::AbstractVars)
  ti = 0.0
  sol1 = zeros(v.sol)
  RHS1 = zeros(v.sol)
  RHS2 = zeros(v.sol)
  RHS3 = zeros(v.sol)
  RHS4 = zeros(v.sol)
  RK4TimeStepper{eltype(v.sol),ndims(v.sol)}(
    0, dt, ti, sol1, RHS1, RHS2, RHS3, RHS4)
end

function RK4TimeStepper(dt::Float64, LC::AbstractArray)
  ti = 0.0
  sol1 = zeros(LC)
  RHS1 = zeros(LC)
  RHS2 = zeros(LC)
  RHS3 = zeros(LC)
  RHS4 = zeros(LC)
  RK4TimeStepper{eltype(LC),ndims(LC)}(
    0, dt, ti, sol1, RHS1, RHS2, RHS3, RHS4)
end




function stepforward!(v::AbstractVars, ts::RK4TimeStepper,
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid)

  eq.calcNL!(ts.RHS1, v.sol, v.t, v, p, g)
  ts.RHS1 .+= eq.LC.*v.sol

  # Substep 1
  ts.ti = v.t + 0.5*ts.dt
  @. ts.sol1 = v.sol + (0.5*ts.dt)*ts.RHS1
  eq.calcNL!(ts.RHS2, ts.sol1, v.t, v, p, g)
  @. ts.RHS2 += eq.LC*ts.sol1

  # Substep 2
  @. ts.sol1 = v.sol + (0.5*ts.dt)*ts.RHS2
  eq.calcNL!(ts.RHS3, ts.sol1, v.t, v, p, g)
  @. ts.RHS3 += eq.LC*ts.sol1

  # Substep 3
  ts.ti = v.t + ts.dt
  @. ts.sol1 = v.sol + ts.dt*ts.RHS3
  eq.calcNL!(ts.RHS4, ts.sol1, v.t, v, p, g)
  @. ts.RHS4 += eq.LC*ts.sol1

  # Substep 4 and final step
  @. v.sol += ts.dt*(
       (1.0/6.0)*ts.RHS1 + (1.0/3.0)*ts.RHS2
     + (1.0/3.0)*ts.RHS3 + (1.0/6.0)*ts.RHS4 )

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
  @. ts.RHS += eq.LC.*v.sol

  @. v.sol += ts.dt * (
    (23.0/12.0)*ts.RHS - (16.0/12.0)*ts.RHSm1 + (5.0/12.0)*ts.RHSm2 )

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
    @. ts.RHS += eq.LC.*v.sol

    # Update
    ts.RHSm2 .= ts.RHSm1
    ts.RHSm1 .= ts.RHS

    @. v.sol   += ts.dt*ts.RHS
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
