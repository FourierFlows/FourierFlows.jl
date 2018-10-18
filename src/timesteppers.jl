"""
    stepforward!(prob)

Step forward `prob` one time step.
"""
stepforward!(prob) = stepforward!(prob.sol, prob.clock, prob.timestepper, prob.eqn, prob.vars, prob.params, prob.grid)

"""
    stepforward!(prob, nsteps)

Step forward `prob` for `nsteps`.
"""
function stepforward!(prob, nsteps) 
  for step = 1:nsteps
    stepforward!(prob)
  end
  nothing
end

"""
    stepforward!(prob, diags, nsteps)

Step forward `prob` for `nsteps`, incrementing diagnostics in the array `diags` along the way.
"""
function stepforward!(prob, diags, nsteps)
  for step = 1:nsteps

    stepforward!(prob)

    for diag in diags
      prob.step % diag.freq != 0 || increment!(diag)
    end

  end
  nothing
end

const fullyexplicitsteppers= [
  :ForwardEuler,
  :RK4,
  :AB3,
  :FilteredForwardEuler,
  :FilteredRK4,
  :FilteredAB3
]

isexplicit(stepper) = any(Symbol(stepper) .== fullyexplicitsteppers)

"""
    TimeStepper(stepper, eq, dt=nothing)
    
Generalized timestepper constructor. If `stepper` is explicit, `dt` is not used.
"""
function TimeStepper(stepper, eq, dt=nothing)
  fullsteppername = Symbol(stepper, :TimeStepper)
  args = isexplicit(stepper) ? (eq, dt) : eq
  eval(Expr(:call, fullsteppername, args...))
end


# The following time-steppers are implemented below
#
#   * Forward Euler
#   * Filtered Forward Euler
#   * RK4
#   * Filtered RK4
#   * ETDRK4
#   * Filtered ETDRK4
#   * AB3
#   * Filtered AB3

# --
# Forward Euler
# --

"""
    ForwardEulerTimeStepper(eq)

Initialize a forward Euler timestepper.
"""
struct ForwardEulerTimeStepper{T} <: AbstractTimeStepper{T}
  N::T # Explicit linear and nonlinear terms
  ForwardEulerTimeStepper(N::T) where T = new{T}(0*N)
end

function stepforward!(sol, cl, ts::ForwardEulerTimeStepper, eq, v, p, g)
  eq.calcN!(ts.N, sol, cl.t, cl, v, p, g)
  @. sol += cl.dt*(eq.L*sol + ts.N)
  cl.t += cl.dt
  cl.step += 1
  nothing
end

ForwardEulerTimeStepper(eq::Equation) = ForwardEulerTimeStepper(superzeros(eq.T, eq.dims))

"""
    FilteredForwardEulerTimeStepper(eq; filterkwargs...)

Construct a forward Euler timestepper with spectral filtering.
"""
struct FilteredForwardEulerTimeStepper{T,Tf} <: AbstractTimeStepper{T}
  N::T       
  filter::Tf
end

function FilteredForwardEulerTimeStepper(eq; filterkwargs...)
  filter = makefilter(eq; filterkwargs...)
  FilteredForwardEulerTimeStepper(superzeros(eq.T, eq.dims), filter)
end

function stepforward!(sol, cl, ts::FilteredForwardEulerTimeStepper, eq, v, p, g)
  eq.calcN!(ts.N, sol, cl.t, cl, v, p, g)
  @. sol = ts.filter*(sol + cl.dt*(ts.N + eq.L*sol))
  cl.t += cl.dt
  cl.step += 1
  nothing
end


# --
# RK4
# --

const sixth = 1/6
const third = 1/3

"""
    RK4TimeStepper(L, Tsol=cxeltype(L))

Construct a 4th-order Runge-Kutta time stepper.
"""
struct RK4TimeStepper{T} <: AbstractTimeStepper{T}
  sol₁::T 
  RHS₁::T
  RHS₂::T
  RHS₃::T
  RHS₄::T
end

struct FilteredRK4TimeStepper{T,Tf} <: AbstractTimeStepper{T}
  sol₁::T 
  RHS₁::T
  RHS₂::T
  RHS₃::T
  RHS₄::T
  filter::Tf
end

"""
    RK4TimeStepper(eq)

Construct a 4th-order Runge-Kutta time stepper with spectral filtering for the equation `eq`.
"""
function RK4TimeStepper(eq)
  @superzeros eq.T eq.dims N sol₁ RHS₁ RHS₂ RHS₃ RHS₄
  RK4TimeStepper(sol₁, RHS₁, RHS₂, RHS₃, RHS₄)
end

function FilteredRK4TimeStepper(eq; filterkwargs...)
  ts = RK4TimeStepper(eq)
  filter = makefilter(eq; filterkwargs...)
  FilteredRK4TimeStepper(getfield.(Ref(ts), fieldnames(typeof(ts)))..., filter)
end

function RK4substeps!(sol, cl, ts, eq, v, p, g)
  # Substep 1
  eq.calcN!(ts.RHS₁, sol, t, cl, v, p, g)
  @. ts.RHS₁ += eq.L*sol
  # Substep 2
  t1 = cl.t + 0.5*cl.dt
  @. ts.sol₁ = sol + 0.5*cl.dt*ts.RHS₁
  eq.calcN!(ts.RHS₂, ts.sol₁, t2, cl, v, p, g)
  @. ts.RHS₂ += eq.L*ts.sol₁
  # Substep 3
  @. ts.sol₁ = sol + 0.5*cl.dt*ts.RHS₂
  eq.calcN!(ts.RHS₃, ts.sol₁, t2, cl, v, p, g)
  @. ts.RHS₃ += eq.L*ts.sol₁
  # Substep 4
  t2 = cl.t + cl.dt
  @. ts.sol₁ = sol + cl.dt*ts.RHS₃
  eq.calcN!(ts.RHS₄, ts.sol₁, t3, cl, v, p, g)
  @. ts.RHS₄ += eq.L*ts.sol₁
  nothing
end

function RK4update!(sol, ts, dt)
  @. sol += dt*(sixth*ts.RHS₁ + third*ts.RHS₂ + third*ts.RHS₃ + sixth*ts.RHS₄)
  nothing
end

function RK4update!(sol, ts, dt, filter)
  @. sol = filter * (sol + dt*(sixth*ts.RHS₁ + third*ts.RHS₂ + third*ts.RHS₃ + sixth*ts.RHS₄))
  nothing
end

function stepforward!(sol, cl, ts::RK4TimeStepper, eq, v, p, g)
  RK4substeps!(sol, cl, ts, eq, v, p, g)
  RK4update!(sol, ts, dt)
  cl.t += cl.dt
  cl.step += 1
  nothing
end

function stepforward!(sol, cl, ts::FilteredRK4TimeStepper, eq, v, p, g)
  RK4substeps!(sol, cl, ts, eq, v, p, g)
  RK4update!(sol, ts, dt, ts.filter)
  cl.t += cl.dt
  cl.step += 1
  nothing
end



# ------
# ETDRK4
# ------

"""
    ETDRK4TimeStepper(eq, dt)

Construct a 4th-order exponential-time-differencing Runge-Kutta time stepper. The Rolls Royce of timestepping.
"""
struct ETDRK4TimeStepper{T,TL} <: AbstractTimeStepper{T}
  # ETDRK4 coefficents
  ζ::TL 
  α::TL
  β::TL
  Γ::TL
  expLdt::TL
  expLdt2::TL
  sol₁::T
  sol₂::T
  N₁::T 
  N₂::T
  N₃::T
  N₄::T
end

"""
    FilteredETDRK4TimeStepper(eq, dt; filterkwargs...)

Construct a 4th-order exponential-time-differencing Runge-Kutta time stepper with spectral filtering.
"""
struct FilteredETDRK4TimeStepper{T,TL,Tf} <: AbstractTimeStepper{T}
  # ETDRK4 coefficents:
  ζ::TL   
  α::TL
  β::TL
  Γ::TL
  expLdt::TL
  expLdt2::TL
  sol₁::T 
  sol₂::T
  N₁::T 
  N₂::T
  N₃::T
  N₄::T
  filter::Tf
end

function ETDRK4TimeStepper(eq, dt)
  dt = fltype(eq.T)(dt) # ensure dt is correct type.
  expLdt  = @. exp(dt*eq.L)
  expLdt2 = @. exp(0.5*dt*eq.L)
  ζ, α, β, Γ = getetdcoeffs(dt, eq.L)
  @superzeros eq.T eq.dims sol₁ sol₂ N₁ N₂ N₃ N₄
  ETDRK4TimeStepper(ζ, α, β, Γ, expLdt, expLdt2, sol₁, sol₂, N₁, N₂, N₃, N₄)
end

function FilteredETDRK4TimeStepper(eq, dt; filterkwargs...)
  ts = ETDRK4TimeStepper(eq, dt)
  filter = makefilter(eq; filterkwargs...)
  FilteredETDRK4TimeStepper(getfield.(Ref(ts), fieldnames(typeof(ts)))..., filter)
end

function ETDRK4update!(sol, ts)
  @. sol = (ts.expLdt*sol +   ts.α * ts.N₁
                          + 2*ts.β * (ts.N₂ + ts.N₃)
                          +   ts.Γ * ts.N₄ )
  nothing
end

function ETDRK4update!(sol, ts, filter)
  @. sol = filter*(ts.expLdt*sol +   ts.α * ts.N₁
                                 + 2*ts.β * (ts.N₂ + ts.N₃)
                                 +   ts.Γ * ts.N₄ )
  nothing
end

function ETDRK4substeps!(sol, cl, ts, eq, v, p, g)
  # Substep 1
  eq.calcN!(ts.N₁, sol, cl.t, cl, v, p, g)
  @. ts.sol₁ = ts.expLdt2*sol + ts.ζ*ts.N₁
  # Substep 2
  t2 = cl.t + 0.5*cl.dt
  eq.calcN!(ts.N₂, ts.sol₁, t2, cl, v, p, g)
  @. ts.sol₂ = ts.expLdt2*sol + ts.ζ*ts.N₂
  # Substep 3
  eq.calcN!(ts.N₃, ts.sol₂, t2, cl, v, p, g)
  @. ts.sol₂ = ts.expLdt2*ts.sol₁ + ts.ζ*(2*ts.N₃ - ts.N₁)
  # Substep 4
  t3 = cl.t + cl.dt
  eq.calcN!(ts.N₄, ts.sol₂, t3, cl, v, p, g)
  nothing
end

function stepforward!(sol, cl, ts::ETDRK4TimeStepper, eq, v, p, g)
  ETDRK4substeps(sol, cl, ts, eq, v, p, g)
  ETDRK4update(sol, ts) # update
  cl.t += cl.dt
  cl.step += 1
  nothing
end

function stepforward!(sol, cl, ts::FilteredETDRK4TimeStepper, eq, v, p, g)
  ETDRK4substeps!(sol, cl, ts, eq, v, p, g)
  ETDRK4update!(sol, ts, ts.filter) # update
  cl.t += cl.dt
  cl.step += 1
  nothing
end

# ---
# AB3
# ---

"""
    AB3TimeStepper(eq)

Construct a 3rd order Adams-Bashforth time stepper.
"""

const ab3h1 = 23/12
const ab3h2 = 16/12
const ab3h3 = 5/12

struct AB3TimeStepper{T} <: AbstractTimeStepper{T}
  RHS::T
  RHS₋₁::T
  RHS₋₂::T
end

function AB3TimeStepper(eq)
  @superzeros eq.T eq.dims RHS RHS₋₁ RHS₋₂
  AB3TimeStepper(RHS, RHS₋₁, RHS₋₂)
end

"""
    FilteredAB3TimeStepper(eq; filterkwargs...)_

Construct a 3rd order Adams-Bashforth time stepper with spectral filtering.
"""
struct FilteredAB3TimeStepper{T,Tf} <: AbstractTimeStepper{T}
  RHS::T
  RHS₋₁::T
  RHS₋₂::T
  filter::Tf
end

function FilteredAB3TimeStepper(eq; filterkwargs...)
  ts = AB3TimeStepper(eq)
  filter = makefilter(eq; filterkwargs...)
  FilteredAB3TimeStepper(getfield.(Ref(ts), fieldnames(typeof(ts)))..., filter)
end


function stepforward!(sol, cl, ts::AB3TimeStepper, eq, v, p, g)
  eq.calcN!(ts.RHS, sol, cl.t, cl, v, p, g)

  @. ts.RHS += eq.L*sol   # Add linear term to RHS

  if cl.step < 3  # forward Euler steps to initialize AB3
    @. sol += cl.dt*ts.RHS    # Update
  else   # Otherwise, stepforward with 3rd order Adams Bashforth:
    @. sol += cl.dt*(ab3h1*ts.RHS - ab3h2*ts.RHS₋₁ + ab3h3*ts.RHS₋₂)
  end

  cl.t += cl.dt
  cl.step += 1

  @. ts.RHS₋₂ = ts.RHS₋₁          # Store
  @. ts.RHS₋₁ = ts.RHS            # ... previous values of RHS

  nothing
end

function stepforward!(sol, cl, ts::FilteredAB3TimeStepper, eq, v, p, g)
  eq.calcN!(ts.RHS, sol, cl.t, cl, v, p, g)
  @. ts.RHS += eq.L*sol   # Add linear term to RHS

  if s.step < 3  # forward Euler steps to initialize AB3
    @. sol = ts.filter*(sol + cl.dt*ts.RHS)    # Update
  else   # Otherwise, stepforward with 3rd order Adams Bashforth:
    @. sol = ts.filter*(sol + cl.dt*(ab3h1*ts.RHS - ab3h2*ts.RHS₋₁ + ab3h3*ts.RHS₋₂))
  end

  cl.t += cl.dt
  cl.step += 1

  @. ts.RHS₋₂ = ts.RHS₋₁          # Store
  @. ts.RHS₋₁ = ts.RHS            # ... previous values of RHS
  nothing
end

# --
# Timestepper utils
# --

"""
    getetdcoeffs(dt, L; ncirc=32, rcirc=1)

Calculate ETDRK4 coefficients associated with the (diagonal) linear coefficient
L by integrating over a small circle in complex space.
"""
function getetdcoeffs(dt, L; ncirc=32, rcirc=1)

  shape = Tuple(cat(ncirc, ones(Int, ndims(L)), dims=1))

  circ = zeros(Complex{Float64}, shape) # use double precision for this calculation
  circ .= rcirc * exp.(2π*im/ncirc*(0.5:1:(ncirc-0.5)))
  circ = permutedims(circ, ndims(circ):-1:1)

  zc = @. dt*L + circ
  M = ndims(L)+1

  # Four coefficients: ζ, α, β, Γ
  ζc = @.             ( exp(zc/2)-1 ) / zc
  αc = @. ( -4 - zc + exp(zc)*(4 - 3zc + zc^2) ) / zc^3
  βc = @.    ( 2  + zc + exp(zc)*(-2 + zc) ) / zc^3
  Γc = @. ( -4 - 3zc - zc^2 + exp(zc)*(4 - zc) ) / zc^3

  ζ = dt*dropdims(mean(ζc, dims=M), dims=M)
  α = dt*dropdims(mean(αc, dims=M), dims=M)
  β = dt*dropdims(mean(βc, dims=M), dims=M)
  Γ = dt*dropdims(mean(Γc, dims=M), dims=M)

  if eltype(L) <: Real # this is conservative, but unclear if necessary
    ζ = real.(ζ)
    α = real.(α)
    β = real.(β)
    Γ = real.(Γ)
  end
 
  ζ, α, β, Γ
end
