export ForwardEulerTimeStepper, FilteredForwardEulerTimeStepper,
       AB3TimeStepper,
       RK4TimeStepper,
       ETDRK4TimeStepper, FilteredETDRK4TimeStepper

export stepforward!




# Looping stepforward function ------------------------------------------------
"""
    stepforward!(prob)

Step forward the Problem prob for one timestep.
"""
function stepforward!(prob::Problem)
    stepforward!(prob.state, prob.ts, prob.eqn, prob.vars, prob.params, 
                 prob.grid)
    prob.t = prob.state.t
    prob.step = prob.state.step
end


"""
    stepforward!(prob, nsteps)

Step forward the problem 'prob' for 'nsteps'.
"""
function stepforward!(prob::Problem, nsteps)
  for step = 1:nsteps
    stepforward!(prob)
  end
  nothing
end


"""
    stepforward!(prob, diags, nsteps)

Step forward the problem prob for nsteps while calculating the 
diagnostics in diags.
"""
function stepforward!(prob::Problem, diags::AbstractArray, nsteps)

  # Initialize diagnostics for speed
  for diag in diags
    newnum = ceil(Int, (diag.count+nsteps)/diag.freq)
    if newnum > diag.num
      warn("Resizing diags before stepping forward...")
      resize!(diag, newnum)
    end
  end

  for step = 1:nsteps
    stepforward!(prob)
    for diag in diags
      if (prob.step+1) % diag.freq == 0.0
        increment!(diag)
      end
    end

  end
  nothing
end



# Utilities -------------------------------------------------------------------
"""
Calculate ETDRK4 coefficients by integrating over a small circle
in complex space.
"""
function getetdcoeffs(dt::Float64, LC::Array{Complex{Float64},2};
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
function getetdcoeffs(dt::Float64, LC::Array{Complex{Float64},3};
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

struct ForwardEulerTimeStepper{dim} <: AbstractTimeStepper
  dt::Float64
  N::Array{Complex{Float64},dim}    # Explicit linear and nonlinear terms
end

function ForwardEulerTimeStepper(dt::Float64, sol::AbstractArray)
  N = zeros(sol)
  ForwardEulerTimeStepper{ndims(LC)}(dt, N)
end




function stepforward!(s::State, ts::ForwardEulerTimeStepper,
                      eq::AbstractEquation, v::AbstractVars, p::AbstractParams, 
                      g::AbstractGrid)

  eq.calcN!(ts.N, s.sol, s.t, s, v, p, g)

  @. s.sol += ts.dt*(ts.N + eq.LC*s.sol)
  s.t += ts.dt
  s.step += 1

  nothing
end








# Filtered Forward Euler ------------------------------------------------------
# The simplest time-stepping method in the books. Explicit and 1st-order
# accurate.

struct FilteredForwardEulerTimeStepper{dim} <: AbstractTimeStepper
  dt::Float64
  N::Array{Complex{Float64},dim}    # Explicit linear and nonlinear terms
  filter::Array{Float64,dim}        # Filter for solution
end

function FilteredForwardEulerTimeStepper(dt, LC, g; filterkwargs...)
  N = zeros(LC)
  filter = makefilter(g, typeof(dt), size(LC); filterkwargs...)
  FilteredForwardEulerTimeStepper{ndims(N)}(dt, N, filter)
end

function stepforward!(s::State, ts::FilteredForwardEulerTimeStepper,
                      eq::AbstractEquation, v::AbstractVars, p::AbstractParams, 
                      g::AbstractGrid)
  eq.calcN!(ts.N, s.sol, s.t, s, v, p, g)
  @. s.sol = ts.filter*(s.sol + ts.dt*(ts.N + eq.LC.*s.sol) )
  s.t += ts.dt
  s.step += 1
  nothing
end


# ETDRK4 ----------------------------------------------------------------------
# The Rolls-Royce of time-stepping. Exact treatment of the implicit linear part
# of the equation, explicit and 4th-order accurate integration of nonlinear
# parts of equation.

struct ETDRK4TimeStepper{dim} <: AbstractTimeStepper
  dt::Float64
  LC::Array{Complex{Float64},dim}          # Linear coefficient
  # ETDRK4 coefficents
  ζ::Array{Complex{Float64},dim}
  α::Array{Complex{Float64},dim}
  β::Array{Complex{Float64},dim}
  Γ::Array{Complex{Float64},dim}
  expLCdt::Array{Complex{Float64},dim}     # Precomputed exp(LC*dt)
  expLCdt2::Array{Complex{Float64},dim}    # Precomputed exp(LC*dt/2)

  # Intermediate times, solutions, and nonlinear evaluations
  sol₁::Array{Complex{Float64},dim}
  sol₂::Array{Complex{Float64},dim}
  N₁::Array{Complex{Float64},dim}
  N₂::Array{Complex{Float64},dim}
  N₃::Array{Complex{Float64},dim}
  N₄::Array{Complex{Float64},dim}
end

function ETDRK4TimeStepper(dt, LC)
  expLCdt  = exp.(dt*LC)
  expLCdt2 = exp.(0.5*dt*LC)
  ζ, α, β, Γ = getetdcoeffs(dt, LC)
  @createarrays eltype(LC) size(LC) sol₁ sol₂ N₁ N₂ N₃ N₄
  ETDRK4TimeStepper{ndims(LC)}(dt, LC, ζ, α, β, Γ, expLCdt,
    expLCdt2, sol₁, sol₂, N₁, N₂, N₃, N₄)
end




# Filtered ETDRK4 --------------------------------------------------------------
# The Rolls-Royce of time-stepping. Exact treatment of linear part of
# the equation, explicit and 4th-order accurate integration of nonlinear
# parts of equation.

struct FilteredETDRK4TimeStepper{dim} <: AbstractTimeStepper
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
  sol₁::Array{Complex{Float64}, dim}
  sol₂::Array{Complex{Float64}, dim}
  N₁::Array{Complex{Float64}, dim}
  N₂::Array{Complex{Float64}, dim}
  N₃::Array{Complex{Float64}, dim}
  N₄::Array{Complex{Float64}, dim}
  filter::Array{Complex{Float64}, dim}    # Filter for solution
end


function FilteredETDRK4TimeStepper(dt, LC, g; filterkwargs...)
  expLCdt  = exp.(dt*LC)
  expLCdt2 = exp.(0.5*dt*LC)
  ζ, α, β, Γ = getetdcoeffs(dt, LC)

  @createarrays eltype(LC) size(LC) sol₁ sol₂ N₁ N₂ N₃ N₄

  filter = makefilter(g, typeof(dt), size(LC); filterkwargs...)

  FilteredETDRK4TimeStepper{ndims(LC)}(dt, LC, ζ, α, β, Γ,
        expLCdt, expLCdt2, sol₁, sol₂, N₁, N₂, N₃, N₄, filter)
end



function stepforward!(s::State, ts::ETDRK4TimeStepper,
                      eq::AbstractEquation, v::AbstractVars, p::AbstractParams, 
                      g::AbstractGrid)
  # Substep 1
  eq.calcN!(ts.N₁, s.sol, s.t, s, v, p, g)
  @. ts.sol₁ = ts.expLCdt2*s.sol + ts.ζ*ts.N₁
  # Substep 2
  t2 = s.t + 0.5*ts.dt
  eq.calcN!(ts.N₂, ts.sol₁, t2, s, v, p, g)
  @. ts.sol₂ = ts.expLCdt2*s.sol + ts.ζ*ts.N₂
  # Substep 3
  eq.calcN!(ts.N₃, ts.sol₂, t2, s, v, p, g)
  @. ts.sol₂ = ts.expLCdt2*ts.sol₁ + ts.ζ*(2.0*ts.N₃ - ts.N₁)
  # Substep 4
  t3 = s.t + ts.dt
  eq.calcN!(ts.N₄, ts.sol₂, t3, s, v, p, g)

  # Update
  @. s.sol = (ts.expLCdt.*s.sol +     ts.α * ts.N₁
                                + 2.0*ts.β * (ts.N₂ + ts.N₃)
                                +     ts.Γ * ts.N₄ )
  s.t += ts.dt
  s.step += 1

  nothing
end

function stepforward!(s::State, ts::FilteredETDRK4TimeStepper,
                      eq::AbstractEquation, v::AbstractVars, p::AbstractParams, 
                      g::AbstractGrid)
  # Substep 1
  eq.calcN!(ts.N₁, s.sol, s.t, s, v, p, g)
  @. ts.sol₁ = ts.expLCdt2*s.sol + ts.ζ*ts.N₁
  # Substep 2
  t2 = s.t + 0.5*ts.dt
  eq.calcN!(ts.N₂, ts.sol₁, t2, s, v, p, g)
  @. ts.sol₂ = ts.expLCdt2*s.sol + ts.ζ*ts.N₂
  # Substep 3
  eq.calcN!(ts.N₃, ts.sol₂, t2, s, v, p, g)
  @. ts.sol₂ = ts.expLCdt2*ts.sol₁ + ts.ζ*(2.0*ts.N₃ - ts.N₁)
  # Substep 4
  t3 = s.t + ts.dt
  eq.calcN!(ts.N₄, ts.sol₂, t3, s, v, p, g)

  # Update
  @. s.sol = ts.filter*(ts.expLCdt*s.sol +     ts.α * ts.N₁
                                         + 2.0*ts.β * (ts.N₂ + ts.N₃)
                                         +     ts.Γ * ts.N₄ )
  s.t += ts.dt
  s.step += 1

  nothing
end



struct DualETDRK4TimeStepper{dimc, dimr} <: AbstractTimeStepper
  dt::Float64
  c::ETDRK4TimeStepper{dimc}
  r::ETDRK4TimeStepper{dimr}
end

function ETDRK4TimeStepper(dt, LCc, LCr)
  c = ETDRK4TimeStepper(dt, LCc)
  r = ETDRK4TimeStepper(dt, LCr)
  DualETDRK4TimeStepper{ndims(LCc), ndims(LCr)}(dt, c, r)
end

function stepforward!(s::DualState, ts::DualETDRK4TimeStepper,
                      eq::AbstractEquation, v::AbstractVars, p::AbstractParams, 
                      g::AbstractGrid)
  # Substep 1
  eq.calcN!(ts.c.N₁, ts.r.N₁, s.solc, s.solr, s.t, s, v, p, g)
  @. ts.c.sol₁ = ts.c.expLCdt2*s.solc + ts.c.ζ*ts.c.N₁
  @. ts.r.sol₁ = ts.r.expLCdt2*s.solr + ts.r.ζ*ts.r.N₁
  # Substep 2
  t2 = s.t + 0.5*ts.c.dt
  eq.calcN!(ts.c.N₂, ts.r.N₂, ts.c.sol₁, ts.r.sol₁, t2, s, v, p, g)
  @. ts.c.sol₂ = ts.c.expLCdt2*s.solc + ts.c.ζ*ts.c.N₂
  @. ts.r.sol₂ = ts.r.expLCdt2*s.solr + ts.r.ζ*ts.r.N₂
  # Substep 3
  eq.calcN!(ts.c.N₃, ts.r.N₃, ts.c.sol₂, ts.r.sol₂, t2, s, v, p, g)
  @. ts.c.sol₂ = (ts.c.expLCdt2*ts.c.sol₁
    + ts.c.ζ*(2.0*ts.c.N₃ - ts.c.N₁))
  @. ts.r.sol₂ = (ts.r.expLCdt2*ts.r.sol₁
    + ts.r.ζ*(2.0*ts.r.N₃ - ts.r.N₁))
  # Substep 4
  t3 = s.t + ts.c.dt
  eq.calcN!(ts.c.N₄, ts.r.N₄, ts.c.sol₂, ts.r.sol₂, t3, s, v, p, g)

  # Update
  @. s.solc = (ts.c.expLCdt.*s.solc +     ts.c.α * ts.c.N₁
                                    + 2.0*ts.c.β * (ts.c.N₂ + ts.c.N₃)
                                    +     ts.c.Γ * ts.c.N₄ )
  @. s.solr = (ts.r.expLCdt.*s.solr +     ts.r.α * ts.r.N₁
                                    + 2.0*ts.r.β * (ts.r.N₂ + ts.r.N₃)
                                    +     ts.r.Γ * ts.r.N₄ )
  s.t += ts.dt
  s.step += 1

  nothing
end


# RK4 -------------------------------------------------------------------------
# RK4 is the classical explicit 4th-order Runge-Kutta time-stepping
# method. It uses a series of substeps/estimators to achieve 4th-order
# accuracy over each individual time-step, at the cost of requiring
# relatively more evaluations of the nonlinear right hand side.
# It is described, among other places, in Bewley's Numerical
# Renaissance.

struct RK4TimeStepper{T,dim} <: AbstractTimeStepper
  dt::Float64
  sol₁::Array{T,dim}
  RHS₁::Array{T,dim}
  RHS₂::Array{T,dim}
  RHS₃::Array{T,dim}
  RHS₄::Array{T,dim}
end

function RK4TimeStepper(dt, LC)
  @createarrays eltype(LC) size(LC) sol₁ RHS₁ RHS₂ RHS₃ RHS₄
  RK4TimeStepper{eltype(LC),ndims(LC)}(dt, sol₁, RHS₁, RHS₂, RHS₃, RHS₄)
end

function stepforward!(s::State, ts::RK4TimeStepper,
                      eq::AbstractEquation, v::AbstractVars, p::AbstractParams, 
                      g::AbstractGrid)
  eq.calcN!(ts.RHS₁, s.sol, s.t, s, v, p, g)
  @. ts.RHS₁ += eq.LC*s.sol
  # Substep 1
  t2 = s.t + 0.5*ts.dt
  @. ts.sol₁ = s.sol + 0.5*ts.dt*ts.RHS₁
  eq.calcN!(ts.RHS₂, ts.sol₁, t2, s, v, p, g)
  @. ts.RHS₂ += eq.LC*ts.sol₁
  # Substep 2
  @. ts.sol₁ = s.sol + 0.5*ts.dt*ts.RHS₂
  eq.calcN!(ts.RHS₃, ts.sol₁, t2, s, v, p, g)
  @. ts.RHS₃ += eq.LC*ts.sol₁
  # Substep 3
  t3 = s.t + ts.dt
  @. ts.sol₁ = s.sol + ts.dt*ts.RHS₃
  eq.calcN!(ts.RHS₄, ts.sol₁, t3, s, v, p, g)
  @. ts.RHS₄ += eq.LC*ts.sol₁

  # Substep 4 and final step
  @. s.sol += ts.dt*(   1.0/6.0*ts.RHS₁ + 1.0/3.0*ts.RHS₂
                      + 1.0/3.0*ts.RHS₃ + 1.0/6.0*ts.RHS₄ )
  s.t += ts.dt
  s.step += 1

  nothing
end


struct FilteredRK4TimeStepper{T,dim} <: AbstractTimeStepper
  dt::Float64
  sol₁::Array{T,dim}
  RHS₁::Array{T,dim}
  RHS₂::Array{T,dim}
  RHS₃::Array{T,dim}
  RHS₄::Array{T,dim}
  filter::Array{T,dim}    # Filter for solution
end

function FilteredRK4TimeStepper(dt, LC, g; filterkwargs...)
  @createarrays eltype(LC) size(LC) sol₁ RHS₁ RHS₂ RHS₃ RHS₄
  filter = makefilter(g, typeof(dt), size(LC); filterkwargs...)
  FilteredRK4TimeStepper{eltype(LC),ndims(LC)}(dt, sol₁, RHS₁, RHS₂, RHS₃, 
    RHS₄, filter)
end

function stepforward!(s::State, ts::FilteredRK4TimeStepper,
                      eq::AbstractEquation, v::AbstractVars, p::AbstractParams, 
                      g::AbstractGrid)
  eq.calcN!(ts.RHS₁, s.sol, s.t, s, v, p, g)
  @. ts.RHS₁ += eq.LC*s.sol
  # Substep 1
  t2 = s.t + 0.5*ts.dt
  @. ts.sol₁ = s.sol + 0.5*ts.dt*ts.RHS₁
  eq.calcN!(ts.RHS₂, ts.sol₁, t2, s, v, p, g)
  @. ts.RHS₂ += eq.LC*ts.sol₁
  # Substep 2
  @. ts.sol₁ = s.sol + 0.5*ts.dt*ts.RHS₂
  eq.calcN!(ts.RHS₃, ts.sol₁, t2, s, v, p, g)
  @. ts.RHS₃ += eq.LC*ts.sol₁
  # Substep 3
  t3 = s.t + ts.dt
  @. ts.sol₁ = s.sol + ts.dt*ts.RHS₃
  eq.calcN!(ts.RHS₄, ts.sol₁, t3, s, v, p, g)
  @. ts.RHS₄ += eq.LC*ts.sol₁

  # Substep 4 and final step
  @. s.sol = ts.filter*(s.sol + ts.dt*(   1.0/6.0*ts.RHS₁ + 1.0/3.0*ts.RHS₂
                                        + 1.0/3.0*ts.RHS₃ + 1.0/6.0*ts.RHS₄ ))
  s.t += ts.dt
  s.step += 1

  nothing
end



# AB3 -------------------------------------------------------------------------
# 3rd order Adams-Bashforth time stepping is an explicit scheme that uses
# solutions from two previous time-steps to achieve 3rd order accuracy.

struct AB3TimeStepper <: AbstractTimeStepper
  dt::Float64
  RHS::Array{Complex{Float64}, 2}
  RHS₋₁::Array{Complex{Float64}, 2}
  RHS₋₂::Array{Complex{Float64}, 2}
end

function AB3TimeStepper(dt::Float64, sol::Array{Complex{Float64}, 2})
  RHS   = zeros(sol)
  RHS₋₁ = zeros(sol)
  RHS₋₂ = zeros(sol)
  AB3TimeStepper(dt, RHS, RHS₋₁, RHS₋₂)
end

function stepforward!(s::State, ts::AB3TimeStepper,
                      eq::AbstractEquation, v::AbstractVars, p::AbstractParams, 
                      g::AbstractGrid)
  if s.step < 2                 # forward Euler steps to initialize AB3
    eq.calcN!(ts.RHS, s.sol, s.t, s, v, p, g)
    @. ts.RHS += eq.LC.*s.sol   # Add linear term to RHS

    @. s.sol += ts.dt*ts.RHS    # Update
    s.t += ts.dt                # ...
    s.step += 1                 # ...

    ts.RHS₋₂ .= ts.RHS₋₁        # Store
    ts.RHS₋₁ .= ts.RHS          # ... previous values of RHS.
  end

  # Otherwise, stepforward with 3rd order Adams Bashforth:
  eq.calcN!(ts.RHS, s.sol, s.t, s, v, p, g)
  @. ts.RHS += eq.LC*s.sol      # Add linear term to RHS

  @. s.sol += ts.dt*(23.0/12.0*ts.RHS - 16.0/12.0*ts.RHS₋₁ + 5.0/12.0*ts.RHS₋₂)
  s.t += ts.dt
  s.step += 1

  ts.RHS₋₂ .= ts.RHS₋₁          # Store
  ts.RHS₋₁ .= ts.RHS            # ... previous values of RHS

  nothing
end
