"""
    stepforward!(prob)

Step forward the Problem `prob` for one timestep.
"""
function stepforward!(prob::Problem)
    stepforward!(prob.state, prob.ts, prob.eqn, prob.vars, prob.params, prob.grid)
    prob.t = prob.state.t
    prob.step = prob.state.step
end

"""
    stepforward!(prob, nsteps)

Step forward `prob` for `nsteps`.
"""
stepforward!(prob::Problem, nsteps) = for step=1:nsteps; stepforward!(prob); end

"""
    stepforward!(prob, diags, nsteps)

Step forward `prob` for `nsteps`, incrementing diagnostics in the array `diags` along the way.
"""
function stepforward!(prob::Problem, diags::AbstractArray, nsteps)
  for step = 1:nsteps
    stepforward!(prob)
    for diag in diags
      if (prob.step+1) % diag.freq == 0
        increment!(diag)
      end
    end
  end
  nothing
end

"""
    getetdcoeffs(dt, LC; ncirc=32, rcirc=1)

Calculate ETDRK4 coefficients associated with the (diagonal) linear coefficient
LC by integrating over a small circle in complex space.
Note: arbitrary-precision arithmetic might provide a more robust method for
calculating these coefficients.
"""
function getetdcoeffs(dt, LC; ncirc=32, rcirc=1)

  shape = Tuple(cat(ncirc, ones(Int, ndims(LC)), dims=1))

  circ = zeros(cxeltype(LC), shape)
  circ .= rcirc * exp.(2π*im/ncirc*(0.5:1:(ncirc-0.5)))
  circ = permutedims(circ, ndims(circ):-1:1)

  zc = dt*LC .+ circ
  M = ndims(LC)+1

  # Four coefficients: ζ, α, β, Γ
  ζc = @.          ( exp(zc/2)-1 ) / zc
  αc = @. ( -4 - zc + exp(zc)*(4 - 3zc + zc^2) ) / zc^3
  βc = @.    ( 2  + zc + exp(zc)*(-2 + zc) ) / zc^3
  Γc = @. ( -4 - 3zc - zc^2 + exp(zc)*(4 - zc) ) / zc^3

  if eltype(LC) <: Real
    ζ = dt*real.(dropdims(mean(ζc, dims=M), dims=M))
    α = dt*real.(dropdims(mean(αc, dims=M), dims=M))
    β = dt*real.(dropdims(mean(βc, dims=M), dims=M))
    Γ = dt*real.(dropdims(mean(Γc, dims=M), dims=M))
  else
    ζ = dt*dropdims(mean(ζc, dims=M), dims=M)
    α = dt*dropdims(mean(αc, dims=M), dims=M)
    β = dt*dropdims(mean(βc, dims=M), dims=M)
    Γ = dt*dropdims(mean(Γc, dims=M), dims=M)
  end

  ζ, α, β, Γ
end


# ------------
# Timesteppers
# ------------

abstract type AbstractForwardEulerTimeStepper <: AbstractTimeStepper end
abstract type AbstractFilteredForwardEulerTimeStepper <: AbstractTimeStepper end

abstract type AbstractETDRK4TimeStepper <: AbstractTimeStepper end
abstract type AbstractFilteredETDRK4TimeStepper <: AbstractTimeStepper end

abstract type AbstractRK4TimeStepper <: AbstractTimeStepper end
abstract type AbstractFilteredRK4TimeStepper <: AbstractTimeStepper end

abstract type AbstractDualRK4TimeStepper <: AbstractTimeStepper end
abstract type AbstractDualFilteredRK4TimeStepper <: AbstractTimeStepper end

abstract type AbstractDualETDRK4TimeStepper <: AbstractTimeStepper end
abstract type AbstractDualFilteredETDRK4TimeStepper <: AbstractTimeStepper end

abstract type AbstractAB3TimeStepper <: AbstractTimeStepper end
abstract type AbstractFilteredAB3TimeStepper <: AbstractTimeStepper end


# -------------
# Forward Euler
# -------------

"""
    ForwardEulerTimeStepper(dt, LC)

Initialize a forward Euler timestepper. The forward Euler method is
the simplest time-stepping method in the books and is explicit and 1st-order accurate.
"""
struct ForwardEulerTimeStepper{T,dim} <: AbstractForwardEulerTimeStepper
  dt::Float64
  N::Array{T,dim}    # Explicit linear and nonlinear terms
  ForwardEulerTimeStepper(dt, N::Array{T,dim}) where {T,dim} = (
    new{cxeltype(N),ndims(N)}(dt, zeros(cxeltype(N), size(N))))
end

function stepforward!(s, ts::AbstractForwardEulerTimeStepper, eq, v, p, g)
  eq.calcN!(ts.N, s.sol, s.t, s, v, p, g)
  @. s.sol += ts.dt*(ts.N + eq.LC*s.sol)
  s.t += ts.dt
  s.step += 1
  nothing
end

# ----------------------
# Filtered Forward Euler
# ----------------------

"""
    FilteredForwardEulerTimeStepper(dt, LC, g; filterkwargs...)

Initialize a forward Euler timestepper with spectral filtering. The forward Euler method is
the simplest time-stepping method in the books and is explicit and 1st-order accurate.
"""
struct FilteredForwardEulerTimeStepper{dim} <: AbstractFilteredForwardEulerTimeStepper
  dt::Float64
  N::Array{Complex{Float64},dim}    # Explicit linear and nonlinear terms
  filter::Array{Float64,dim}        # Filter for solution
end

function FilteredForwardEulerTimeStepper(dt, LC, g; filterkwargs...)
  @createarrays eltype(LC) size(LC) N
  filter = makefilter(g, typeof(dt), size(LC); filterkwargs...)
  FilteredForwardEulerTimeStepper{ndims(N)}(dt, N, filter)
end

function stepforward!(s, ts::AbstractFilteredForwardEulerTimeStepper, eq, v, p, g)
  eq.calcN!(ts.N, s.sol, s.t, s, v, p, g)
  @. s.sol = ts.filter*(s.sol + ts.dt*(ts.N + eq.LC*s.sol) )
  s.t += ts.dt
  s.step += 1
  nothing
end


# ------
# ETDRK4
# ------

"""
    ETDRK4TimeStepper(dt, LC)

The Rolls-Royce of time-stepping. Exact treatment of the implicit linear part
of the equation, explicit and 4th-order accurate integration of nonlinear
parts of equation.
"""
struct ETDRK4TimeStepper{T,dim} <: AbstractETDRK4TimeStepper
  dt::Float64
  LC::Array{T,dim}          # Linear coefficient
  # ETDRK4 coefficents
  ζ::Array{T,dim}
  α::Array{T,dim}
  β::Array{T,dim}
  Γ::Array{T,dim}
  expLCdt::Array{Complex{Float64},dim}     # Precomputed exp(LC*dt)
  expLCdt2::Array{Complex{Float64},dim}    # Precomputed exp(LC*dt/2)
  # Intermediate times, solutions, and nonlinear evaluations
  sol₁::Array{T,dim}
  sol₂::Array{T,dim}
  N₁::Array{T,dim}
  N₂::Array{T,dim}
  N₃::Array{T,dim}
  N₄::Array{T,dim}
end

function ETDRK4TimeStepper(dt, LC; cxsol=true)
  T = cxsol ? cxeltype(LC) : eltype(LC)
  expLCdt  = exp.(dt*LC)
  expLCdt2 = exp.(0.5*dt*LC)
  ζ, α, β, Γ = getetdcoeffs(dt, LC)
  @createarrays T size(LC) sol₁ sol₂ N₁ N₂ N₃ N₄
  ETDRK4TimeStepper{T,ndims(LC)}(dt, LC, ζ, α, β, Γ, expLCdt, expLCdt2, sol₁, sol₂, N₁, N₂, N₃, N₄)
end

function stepforward!(s, ts::AbstractETDRK4TimeStepper, eq, v, p, g)
  # Substep 1
  eq.calcN!(ts.N₁, s.sol, s.t, s, v, p, g)
  @. ts.sol₁ = ts.expLCdt2*s.sol + ts.ζ*ts.N₁
  # Substep 2
  t2 = s.t + 0.5*ts.dt
  eq.calcN!(ts.N₂, ts.sol₁, t2, s, v, p, g)
  @. ts.sol₂ = ts.expLCdt2*s.sol + ts.ζ*ts.N₂
  # Substep 3
  eq.calcN!(ts.N₃, ts.sol₂, t2, s, v, p, g)
  @. ts.sol₂ = ts.expLCdt2*ts.sol₁ + ts.ζ*(2*ts.N₃ - ts.N₁)
  # Substep 4
  t3 = s.t + ts.dt
  eq.calcN!(ts.N₄, ts.sol₂, t3, s, v, p, g)

  # Update
  @. s.sol = (ts.expLCdt.*s.sol +   ts.α * ts.N₁
                                + 2*ts.β * (ts.N₂ + ts.N₃)
                                +   ts.Γ * ts.N₄ )
  s.t += ts.dt
  s.step += 1

  nothing
end

struct DualETDRK4TimeStepper{Tc,Tr,dimc,dimr} <: AbstractDualETDRK4TimeStepper
  dt::Float64
  c::ETDRK4TimeStepper{Tc,dimc}
  r::ETDRK4TimeStepper{Tr,dimr}
end

function DualETDRK4TimeStepper(dt, LCc, LCr)
  c = ETDRK4TimeStepper(dt, LCc)
  r = ETDRK4TimeStepper(dt, LCr)
  DualETDRK4TimeStepper{cxeltype(LCc),cxeltype(LCr),ndims(LCc),ndims(LCr)}(dt, c, r)
end
ETDRK4TimeStepper(dt, LCc, LCr) = DualETDRK4TimeStepper(dt, LCc, LCr)

function stepforward!(s::DualState, ts::AbstractDualETDRK4TimeStepper, eq, v, p, g)
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
  @. ts.c.sol₂ = (ts.c.expLCdt2*ts.c.sol₁ + ts.c.ζ*(2*ts.c.N₃ - ts.c.N₁))
  @. ts.r.sol₂ = (ts.r.expLCdt2*ts.r.sol₁ + ts.r.ζ*(2*ts.r.N₃ - ts.r.N₁))
  # Substep 4
  t3 = s.t + ts.c.dt
  eq.calcN!(ts.c.N₄, ts.r.N₄, ts.c.sol₂, ts.r.sol₂, t3, s, v, p, g)

  # Update
  @. s.solc = (ts.c.expLCdt.*s.solc +   ts.c.α * ts.c.N₁
                                    + 2*ts.c.β * (ts.c.N₂ + ts.c.N₃)
                                    +   ts.c.Γ * ts.c.N₄ )
  @. s.solr = (ts.r.expLCdt.*s.solr +   ts.r.α * ts.r.N₁
                                    + 2*ts.r.β * (ts.r.N₂ + ts.r.N₃)
                                    +   ts.r.Γ * ts.r.N₄ )
  s.t += ts.dt
  s.step += 1

  nothing
end


# ---------------
# Filtered ETDRK4
# ---------------

struct FilteredETDRK4TimeStepper{T,dim} <: AbstractFilteredETDRK4TimeStepper
  dt::Float64
  LC::Array{T,dim}          # Linear coefficient
  # ETDRK4 coefficents
  ζ::Array{T,dim}
  α::Array{T,dim}
  β::Array{T,dim}
  Γ::Array{T,dim}
  expLCdt::Array{T,dim}     # Precomputed exp(LC*dt)
  expLCdt2::Array{T,dim}    # Precomputed exp(LC*dt/2)
  # Intermediate times, solutions, and nonlinear evaluations
  sol₁::Array{T,dim}
  sol₂::Array{T,dim}
  N₁::Array{T,dim}
  N₂::Array{T,dim}
  N₃::Array{T,dim}
  N₄::Array{T,dim}
  filter::Array{T,dim}    # Filter for solution
end


function FilteredETDRK4TimeStepper(dt, LC, g; cxsol=true, filterkwargs...)
  T = cxsol ? cxeltype(LC) : eltype(LC)
  expLCdt  = exp.(dt*LC)
  expLCdt2 = exp.(0.5*dt*LC)
  ζ, α, β, Γ = getetdcoeffs(dt, LC)
  @createarrays T size(LC) sol₁ sol₂ N₁ N₂ N₃ N₄
  filter = makefilter(g, typeof(dt), size(LC); filterkwargs...)
  FilteredETDRK4TimeStepper{T,ndims(LC)}(dt, LC, ζ, α, β, Γ, expLCdt, expLCdt2, sol₁, sol₂, N₁, N₂, N₃, N₄, filter)
end


function stepforward!(s, ts::AbstractFilteredETDRK4TimeStepper, eq, v, p, g)
  # Substep 1
  eq.calcN!(ts.N₁, s.sol, s.t, s, v, p, g)
  @. ts.sol₁ = ts.expLCdt2*s.sol + ts.ζ*ts.N₁
  # Substep 2
  t2 = s.t + 0.5*ts.dt
  eq.calcN!(ts.N₂, ts.sol₁, t2, s, v, p, g)
  @. ts.sol₂ = ts.expLCdt2*s.sol + ts.ζ*ts.N₂
  # Substep 3
  eq.calcN!(ts.N₃, ts.sol₂, t2, s, v, p, g)
  @. ts.sol₂ = ts.expLCdt2*ts.sol₁ + ts.ζ*(2*ts.N₃ - ts.N₁)
  # Substep 4
  t3 = s.t + ts.dt
  eq.calcN!(ts.N₄, ts.sol₂, t3, s, v, p, g)

  # Update
  @. s.sol = ts.filter*(ts.expLCdt*s.sol +   ts.α * ts.N₁
                                         + 2*ts.β * (ts.N₂ + ts.N₃)
                                         +   ts.Γ * ts.N₄ )
  s.t += ts.dt
  s.step += 1

  nothing
end

struct DualFilteredETDRK4TimeStepper{Tc,Tr,dimc,dimr} <: AbstractDualFilteredETDRK4TimeStepper
  dt::Float64
  c::FilteredETDRK4TimeStepper{Tc,dimc}
  r::FilteredETDRK4TimeStepper{Tr,dimr}
end

function DualFilteredETDRK4TimeStepper(dt, LCc, LCr, g)
  c = FilteredETDRK4TimeStepper(dt, LCc, g)
  r = FilteredETDRK4TimeStepper(dt, LCr, g)
  DualFilteredETDRK4TimeStepper{cxeltype(LCc),cxeltype(LCr),ndims(LCc),ndims(LCr)}(dt, c, r)
end
FilteredETDRK4TimeStepper(dt, LCc, LCr, g) = DualFilteredETDRK4TimeStepper(dt, LCc, LCr, g)

function stepforward!(s::DualState, ts::AbstractDualFilteredETDRK4TimeStepper, eq, v, p, g)
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
  @. ts.c.sol₂ = (ts.c.expLCdt2*ts.c.sol₁ + ts.c.ζ*(2*ts.c.N₃ - ts.c.N₁))
  @. ts.r.sol₂ = (ts.r.expLCdt2*ts.r.sol₁ + ts.r.ζ*(2*ts.r.N₃ - ts.r.N₁))
  # Substep 4
  t3 = s.t + ts.c.dt
  eq.calcN!(ts.c.N₄, ts.r.N₄, ts.c.sol₂, ts.r.sol₂, t3, s, v, p, g)

  # Update
  @. s.solc = ts.c.filter*(ts.c.expLCdt.*s.solc +   ts.c.α * ts.c.N₁
                                    + 2*ts.c.β * (ts.c.N₂ + ts.c.N₃)
                                    +   ts.c.Γ * ts.c.N₄ )
  @. s.solr = ts.r.filter*(ts.r.expLCdt.*s.solr +   ts.r.α * ts.r.N₁
                                    + 2*ts.r.β * (ts.r.N₂ + ts.r.N₃)
                                    +   ts.r.Γ * ts.r.N₄ )
  s.t += ts.dt
  s.step += 1

  nothing
end



# ---
# RK4
# ---

"""
    RK4TimeStepper(dt, LC)

RK4 is the classical explicit 4th-order Runge-Kutta time-stepping
method. It uses a series of substeps/estimators to achieve 4th-order
accuracy over each individual time-step, at the cost of requiring
relatively more evaluations of the nonlinear right hand side.
It is described, among other places, in Bewley's Numerical
Renaissance.
"""
struct RK4TimeStepper{T,dim} <: AbstractRK4TimeStepper
  dt::Float64
  sol₁::Array{T,dim}
  RHS₁::Array{T,dim}
  RHS₂::Array{T,dim}
  RHS₃::Array{T,dim}
  RHS₄::Array{T,dim}
end

function RK4TimeStepper(dt, LC; cxsol=true)
  T = cxsol ? cxeltype(LC) : eltype(LC)
  @createarrays T size(LC) sol₁ RHS₁ RHS₂ RHS₃ RHS₄
  RK4TimeStepper{T,ndims(LC)}(dt, sol₁, RHS₁, RHS₂, RHS₃, RHS₄)
end

function stepforward!(s, ts::AbstractRK4TimeStepper, eq, v, p, g)
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
  @. s.sol += ts.dt*(1/6*ts.RHS₁ + 1/3*ts.RHS₂ + 1/3*ts.RHS₃ + 1/6*ts.RHS₄)
  s.t += ts.dt
  s.step += 1
  nothing
end

struct DualRK4TimeStepper{Tc,Tr,dimc,dimr} <: AbstractDualRK4TimeStepper
  dt::Float64
  c::RK4TimeStepper{Tc,dimc}
  r::RK4TimeStepper{Tr,dimr}
end

function DualRK4TimeStepper(dt, LCc, LCr)
  c = RK4TimeStepper(dt, LCc)
  r = RK4TimeStepper(dt, LCr)
  DualRK4TimeStepper{cxeltype(LCc),cxeltype(LCr),ndims(LCc),ndims(LCr)}(dt, c, r)
end

RK4TimeStepper(dt, LCc, LCr) = DualRK4TimeStepper(dt, LCc, LCr)

function stepforward!(s, ts::AbstractDualRK4TimeStepper, eq, v, p, g)
  eq.calcN!(ts.c.RHS₁, ts.r.RHS₁, s.solc, s.solr, s.t, s, v, p, g)
  @. ts.c.RHS₁ += eq.LCc*s.solc
  @. ts.r.RHS₁ += eq.LCr*s.solr
  # Substep 1
  t2 = s.t + 0.5*ts.c.dt
  @. ts.c.sol₁ = s.solc + 0.5*ts.c.dt*ts.c.RHS₁
  @. ts.r.sol₁ = s.solr + 0.5*ts.r.dt*ts.r.RHS₁
  eq.calcN!(ts.c.RHS₂, ts.r.RHS₂, ts.c.sol₁, ts.r.sol₁, t2, s, v, p, g)
  @. ts.c.RHS₂ += eq.LCc*ts.c.sol₁
  @. ts.r.RHS₂ += eq.LCr*ts.r.sol₁
  # Substep 2
  @. ts.c.sol₁ = s.solc + 0.5*ts.c.dt*ts.c.RHS₂
  @. ts.r.sol₁ = s.solr + 0.5*ts.r.dt*ts.r.RHS₂
  eq.calcN!(ts.c.RHS₃, ts.r.RHS₃, ts.c.sol₁, ts.r.sol₁, t2, s, v, p, g)
  @. ts.c.RHS₃ += eq.LCc*ts.c.sol₁
  @. ts.r.RHS₃ += eq.LCr*ts.r.sol₁
  # Substep 3
  t3 = s.t + ts.c.dt
  @. ts.c.sol₁ = s.solc + ts.c.dt*ts.c.RHS₃
  @. ts.r.sol₁ = s.solr + ts.r.dt*ts.r.RHS₃
  eq.calcN!(ts.c.RHS₄, ts.r.RHS₄, ts.c.sol₁, ts.r.sol₁, t3, s, v, p, g)
  @. ts.c.RHS₄ += eq.LCc*ts.c.sol₁
  @. ts.r.RHS₄ += eq.LCr*ts.r.sol₁

  # Substep 4 and final step
  @. s.solc += ts.c.dt*(1/6*ts.c.RHS₁ + 1/3*ts.c.RHS₂ + 1/3*ts.c.RHS₃ + 1/6*ts.c.RHS₄)
  @. s.solr += ts.r.dt*(1/6*ts.r.RHS₁ + 1/3*ts.r.RHS₂ + 1/3*ts.r.RHS₃ + 1/6*ts.r.RHS₄)
  s.t += ts.c.dt
  s.step += 1
  nothing
end





# ------------
# Filtered RK4
# ------------

struct FilteredRK4TimeStepper{T,dim} <: AbstractFilteredRK4TimeStepper
  dt::Float64
  sol₁::Array{T,dim}
  RHS₁::Array{T,dim}
  RHS₂::Array{T,dim}
  RHS₃::Array{T,dim}
  RHS₄::Array{T,dim}
  filter::Array{T,dim}    # Filter for solution
end

function FilteredRK4TimeStepper(dt, LC, g; cxsol=true, filterkwargs...)
  T = cxsol ? cxeltype(LC) : eltype(LC)
  @createarrays T size(LC) sol₁ RHS₁ RHS₂ RHS₃ RHS₄
  filter = makefilter(g, typeof(dt), size(LC); filterkwargs...)
  FilteredRK4TimeStepper{T,ndims(LC)}(dt, sol₁, RHS₁, RHS₂, RHS₃, RHS₄, filter)
end

function stepforward!(s, ts::AbstractFilteredRK4TimeStepper, eq, v, p, g)
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
  @. s.sol = ts.filter*(s.sol + ts.dt*(1/6*ts.RHS₁ + 1/3*ts.RHS₂ + 1/3*ts.RHS₃ + 1/6*ts.RHS₄))
  s.t += ts.dt
  s.step += 1

  nothing
end

struct DualFilteredRK4TimeStepper{Tc,Tr,dimc,dimr} <: AbstractDualFilteredRK4TimeStepper
  dt::Float64
  c::FilteredRK4TimeStepper{Tc,dimc}
  r::FilteredRK4TimeStepper{Tr,dimr}
end

function DualFilteredRK4TimeStepper(dt, LCc, LCr, g)
  c = FilteredRK4TimeStepper(dt, LCc, g)
  r = FilteredRK4TimeStepper(dt, LCr, g)
  DualFilteredRK4TimeStepper{cxeltype(LCc),cxeltype(LCr),ndims(LCc),ndims(LCr)}(dt, c, r)
end
FilteredRK4TimeStepper(dt, LCc, LCr, g) = DualFilteredRK4TimeStepper(dt, LCc, LCr, g)

function stepforward!(s, ts::AbstractDualFilteredRK4TimeStepper, eq, v, p, g)
  eq.calcN!(ts.c.RHS₁, ts.r.RHS₁, s.solc, s.solr, s.t, s, v, p, g)
  @. ts.c.RHS₁ += eq.LCc*s.solc
  @. ts.r.RHS₁ += eq.LCr*s.solr
  # Substep 1
  t2 = s.t + 0.5*ts.c.dt
  @. ts.c.sol₁ = s.solc + 0.5*ts.c.dt*ts.c.RHS₁
  @. ts.r.sol₁ = s.solr + 0.5*ts.r.dt*ts.r.RHS₁
  eq.calcN!(ts.c.RHS₂, ts.r.RHS₂, ts.c.sol₁, ts.r.sol₁, t2, s, v, p, g)
  @. ts.c.RHS₂ += eq.LCc*ts.c.sol₁
  @. ts.r.RHS₂ += eq.LCr*ts.r.sol₁
  # Substep 2
  @. ts.c.sol₁ = s.solc + 0.5*ts.c.dt*ts.c.RHS₂
  @. ts.r.sol₁ = s.solr + 0.5*ts.r.dt*ts.r.RHS₂
  eq.calcN!(ts.c.RHS₃, ts.r.RHS₃, ts.c.sol₁, ts.r.sol₁, t2, s, v, p, g)
  @. ts.c.RHS₃ += eq.LCc*ts.c.sol₁
  @. ts.r.RHS₃ += eq.LCr*ts.r.sol₁
  # Substep 3
  t3 = s.t + ts.c.dt
  @. ts.c.sol₁ = s.solc + ts.c.dt*ts.c.RHS₃
  @. ts.r.sol₁ = s.solr + ts.r.dt*ts.r.RHS₃
  eq.calcN!(ts.c.RHS₄, ts.r.RHS₄, ts.c.sol₁, ts.r.sol₁, t3, s, v, p, g)
  @. ts.c.RHS₄ += eq.LCc*ts.c.sol₁
  @. ts.r.RHS₄ += eq.LCr*ts.r.sol₁

  # Substep 4 and final step
  @. s.solc = ts.c.filter*(s.solc + ts.c.dt*(1/6*ts.c.RHS₁ + 1/3*ts.c.RHS₂ + 1/3*ts.c.RHS₃ + 1/6*ts.c.RHS₄))
  @. s.solr = ts.r.filter*(s.solr + ts.r.dt*(1/6*ts.r.RHS₁ + 1/3*ts.r.RHS₂ + 1/3*ts.r.RHS₃ + 1/6*ts.r.RHS₄))
  s.t += ts.c.dt
  s.step += 1
  nothing
end




# ---
# AB3
# ---

"""
    AB3TimeStepper(dt, LC)

3rd order Adams-Bashforth time stepping is an explicit scheme that uses
solutions from two previous time-steps to achieve 3rd order accuracy.
"""
struct AB3TimeStepper{T,dim} <: AbstractAB3TimeStepper
  dt::Float64
  RHS::Array{T,dim}
  RHS₋₁::Array{T,dim}
  RHS₋₂::Array{T,dim}
end

function AB3TimeStepper(dt, LC; cxsol=true)
  T = cxsol ? cxeltype(LC) : eltype(LC)
  @createarrays T size(LC) RHS RHS₋₁ RHS₋₂
  AB3TimeStepper{T,ndims(LC)}(dt, RHS, RHS₋₁, RHS₋₂)
end

function stepforward!(s, ts::AbstractAB3TimeStepper, eq, v, p, g)
  eq.calcN!(ts.RHS, s.sol, s.t, s, v, p, g)
  @. ts.RHS += eq.LC.*s.sol   # Add linear term to RHS

  if s.step < 3  # forward Euler steps to initialize AB3
    @. s.sol += ts.dt*ts.RHS    # Update
  else   # Otherwise, stepforward with 3rd order Adams Bashforth:
    @. s.sol += ts.dt*(23/12*ts.RHS - 16/12*ts.RHS₋₁ + 5/12*ts.RHS₋₂)
  end

  s.t += ts.dt
  s.step += 1
  ts.RHS₋₂ .= ts.RHS₋₁          # Store
  ts.RHS₋₁ .= ts.RHS            # ... previous values of RHS
  nothing
end


# ------------
# Filtered AB3
# ------------

"""
    FilteredAB3TimeStepper(dt, LC, g; filterkwargs...)

3rd order Adams-Bashforth time stepping is an explicit scheme that uses
solutions from two previous time-steps to achieve 3rd order accuracy.
"""
struct FilteredAB3TimeStepper{T,dim} <: AbstractFilteredAB3TimeStepper
  dt::Float64
  RHS::Array{T,dim}
  RHS₋₁::Array{T,dim}
  RHS₋₂::Array{T,dim}
  filter::Array{T,dim}    # Filter for solution
end

function FilteredAB3TimeStepper(dt, LC, g; cxsol=true, filterkwargs...)
  T = cxsol ? cxeltype(LC) : eltype(LC)
  @createarrays T size(LC) RHS RHS₋₁ RHS₋₂
  filter = makefilter(g, typeof(dt), size(LC); filterkwargs...)
  FilteredAB3TimeStepper{T,ndims(LC)}(dt, RHS, RHS₋₁, RHS₋₂, filter)
end

function stepforward!(s, ts::AbstractFilteredAB3TimeStepper, eq, v, p, g)
  eq.calcN!(ts.RHS, s.sol, s.t, s, v, p, g)
  @. ts.RHS += eq.LC.*s.sol   # Add linear term to RHS

  if s.step < 3  # forward Euler steps to initialize AB3
    @. s.sol = ts.filter*(s.sol + ts.dt*ts.RHS)    # Update
  else   # Otherwise, stepforward with 3rd order Adams Bashforth:
    @. s.sol = ts.filter*(s.sol + ts.dt*(23/12*ts.RHS - 16/12*ts.RHS₋₁ + 5/12*ts.RHS₋₂))
  end

  s.t += ts.dt
  s.step += 1
  ts.RHS₋₂ .= ts.RHS₋₁          # Store
  ts.RHS₋₁ .= ts.RHS            # ... previous values of RHS
  nothing
end
