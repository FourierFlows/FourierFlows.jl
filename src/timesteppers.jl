"""
    stepforward!(prob)

Step forward `prob` one time step.
"""
function stepforward!(prob::Problem)
  stepforward!(prob.sol, prob.clock, prob.timestepper, prob.eqn, prob.vars, prob.params, prob.grid)
end

"""
    stepforward!(prob, nsteps)

Step forward `prob` for `nsteps`.
"""
function stepforward!(prob::Problem, nsteps::Int)
  for step = 1:nsteps
    stepforward!(prob)
  end
  return nothing
end

"""
    stepforward!(prob, diags, nsteps)

Step forward `prob` for `nsteps`, incrementing `diags` along the way. `diags` may be a single `Diagnostic` or
a `Vector` of `Diagnostic`s.
"""
function stepforward!(prob::Problem, diags, nsteps::Int)
  for step = 1:nsteps
    stepforward!(prob)
    increment!(diags)
  end
  return nothing
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
    TimeStepper(stepper, eq, dt=nothing, dev=CPU(); kw...)

Instantiate the Time`stepper` for `eq`uation with timestep `dt` and
on the `dev`ice. The `kw` are passed to the timestepper constructor.
"""
function TimeStepper(stepper, eq, dt=nothing, dev::Device=CPU(); kw...)
  fullsteppername = Symbol(stepper, :TimeStepper)

  # Create expression that instantiates the time-stepper, depending on whether 
  # timestepper is explicit or not.
  expr = isexplicit(stepper) ? Expr(:call, fullsteppername, eq, dev) :
                               Expr(:call, fullsteppername, eq, dt, dev)

  # Add keyword arguments    
  length(kw) > 0 && push!(expr.args, Tuple(Expr(:kw, p.first, p.second) for p in kw)...)

  return eval(expr)
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
#
# Explicit time-steppers are constructed with the signature
#   ts = ExplicitTimeStepper(eq::Equation)
#
# Implicit time-steppers are constructed with the signature
#   ts = ImplicitTimeStepper(eq::Equation, dt)

# --
# Forward Euler
# --

"""
    ForwardEulerTimeStepper(eq, dev)

Initialize a forward Euler timestepper.
"""
struct ForwardEulerTimeStepper{T} <: AbstractTimeStepper{T}
  N::T # Explicit linear and nonlinear terms
  ForwardEulerTimeStepper(N::T) where T = new{T}(0N)
end

ForwardEulerTimeStepper(eq::Equation, dev::Device=CPU()) = ForwardEulerTimeStepper(devzeros(dev, eq.T, eq.dims))

function stepforward!(sol, cl, ts::ForwardEulerTimeStepper, eq, v, p, g)
  eq.calcN!(ts.N, sol, cl.t, cl, v, p, g)
  @. sol += cl.dt*(eq.L*sol + ts.N)
  cl.t += cl.dt
  cl.step += 1
  return nothing
end

"""
    FilteredForwardEulerTimeStepper(eq, dev; filterkwargs...)

Construct a forward Euler timestepper with spectral filtering.
"""
struct FilteredForwardEulerTimeStepper{T,Tf} <: AbstractTimeStepper{T}
  N::T
  filter::Tf
end

function FilteredForwardEulerTimeStepper(eq::Equation, dev::Device=CPU(); filterkwargs...)
  filter = makefilter(eq; filterkwargs...)
  FilteredForwardEulerTimeStepper(devzeros(dev, eq.T, eq.dims), filter)
end


function stepforward!(sol, cl, ts::FilteredForwardEulerTimeStepper, eq, v, p, g)
  eq.calcN!(ts.N, sol, cl.t, cl, v, p, g)
  @. sol = ts.filter*(sol + cl.dt*(ts.N + eq.L*sol))
  cl.t += cl.dt
  cl.step += 1
  return nothing
end


# --
# RK4
# --

"""
    RK4TimeStepper(eq, dev)

Construct a 4th-order Runge-Kutta time stepper.
"""
struct RK4TimeStepper{T} <: AbstractTimeStepper{T}
  sol₁::T
  RHS₁::T
  RHS₂::T
  RHS₃::T
  RHS₄::T
end

"""
    FilteredRK4TimeStepper(eq, dev; filterkwargs...)

Construct a 4th-order Runge-Kutta time stepper with spectral filtering for the equation `eq`.
"""
struct FilteredRK4TimeStepper{T,Tf} <: AbstractTimeStepper{T}
  sol₁::T
  RHS₁::T
  RHS₂::T
  RHS₃::T
  RHS₄::T
  filter::Tf
end

function RK4TimeStepper(eq::Equation, dev::Device=CPU())
  @devzeros typeof(dev) eq.T eq.dims N sol₁ RHS₁ RHS₂ RHS₃ RHS₄
  RK4TimeStepper(sol₁, RHS₁, RHS₂, RHS₃, RHS₄)
end

function FilteredRK4TimeStepper(eq::Equation, dev::Device=CPU(); filterkwargs...)
  ts = RK4TimeStepper(eq, dev)
  filter = makefilter(eq; filterkwargs...)
  FilteredRK4TimeStepper(getfield.(Ref(ts), fieldnames(typeof(ts)))..., filter)
end

function addlinearterm!(RHS, L, sol)
  @. RHS += L*sol
  return nothing
end

function substepsol!(newsol, sol, RHS, dt)
  @. newsol = sol + dt*RHS
  return nothing
end

function RK4substeps!(sol, cl, ts, eq, v, p, g, t, dt)
  # Substep 1
  eq.calcN!(ts.RHS₁, sol, t, cl, v, p, g)
  addlinearterm!(ts.RHS₁, eq.L, sol)
  # Substep 2
  substepsol!(ts.sol₁, sol, ts.RHS₁, dt/2)
  eq.calcN!(ts.RHS₂, ts.sol₁, t+dt/2, cl, v, p, g)
  addlinearterm!(ts.RHS₂, eq.L, ts.sol₁)
  # Substep 3
  substepsol!(ts.sol₁, sol, ts.RHS₂, dt/2)
  eq.calcN!(ts.RHS₃, ts.sol₁, t+dt/2, cl, v, p, g)
  addlinearterm!(ts.RHS₃, eq.L, ts.sol₁)
  # Substep 4
  substepsol!(ts.sol₁, sol, ts.RHS₃, dt)
  eq.calcN!(ts.RHS₄, ts.sol₁, t+dt, cl, v, p, g)
  addlinearterm!(ts.RHS₄, eq.L, ts.sol₁)
  return nothing
end

function RK4update!(sol, RHS₁, RHS₂, RHS₃, RHS₄, dt)
  @. sol += dt*(RHS₁ / 6 + RHS₂ / 3  + RHS₃ / 3 + RHS₄ / 6)
  return nothing
end

function RK4update!(sol, RHS₁, RHS₂, RHS₃, RHS₄, filter, dt)
  @. sol = filter * (sol + dt*(RHS₁ / 6 + RHS₂ / 3  + RHS₃ / 3 + RHS₄ / 6))
  return nothing
end

function stepforward!(sol, cl, ts::RK4TimeStepper, eq, v, p, g)
  RK4substeps!(sol, cl, ts, eq, v, p, g, cl.t, cl.dt)
  RK4update!(sol, ts.RHS₁, ts.RHS₂, ts.RHS₃, ts.RHS₄, cl.dt)
  cl.t += cl.dt
  cl.step += 1
  return nothing
end

function stepforward!(sol, cl, ts::FilteredRK4TimeStepper, eq, v, p, g)
  RK4substeps!(sol, cl, ts, eq, v, p, g, cl.t, cl.dt)
  RK4update!(sol, ts.RHS₁, ts.RHS₂, ts.RHS₃, ts.RHS₄, ts.filter, cl.dt)
  cl.t += cl.dt
  cl.step += 1
  return nothing
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

function ETDRK4TimeStepper(eq::Equation, dt, dev::Device=CPU())
  dt = fltype(eq.T)(dt) # ensure dt is correct type.
  expLdt, expLdt2 = getexpLs(dt, eq)
  ζ, α, β, Γ = getetdcoeffs(dt, eq.L)
  @devzeros typeof(dev) eq.T eq.dims sol₁ sol₂ N₁ N₂ N₃ N₄
  ETDRK4TimeStepper(ζ, α, β, Γ, expLdt, expLdt2, sol₁, sol₂, N₁, N₂, N₃, N₄)
end

function FilteredETDRK4TimeStepper(eq::Equation, dt, dev::Device=CPU(); filterkwargs...)
  ts = ETDRK4TimeStepper(eq, dt, dev)
  filter = makefilter(eq; filterkwargs...)
  FilteredETDRK4TimeStepper(getfield.(Ref(ts), fieldnames(typeof(ts)))..., filter)
end

function ETDRK4update!(sol, expLdt, α, β, Γ, N₁, N₂, N₃, N₄)
  @. sol = (expLdt*sol +  α * N₁
                       + 2β * (N₂ + N₃)
                       +  Γ * N₄ )
  return nothing
end

function ETDRK4update!(sol, ts, filter)
  @. sol = filter*(ts.expLdt*sol +   ts.α * ts.N₁
                                 + 2*ts.β * (ts.N₂ + ts.N₃)
                                 +   ts.Γ * ts.N₄ )
  return nothing
end

function ETDRK4substep12!(sol₁, expLdt2, sol, ζ, N)
  @. sol₁ = expLdt2*sol + ζ*N
  return nothing
end

function ETDRK4substep3!(sol₂, expLdt2, sol₁, ζ, N₁, N₃)
  @. sol₂ = expLdt2*sol₁ + ζ*(2N₃ - N₁)
  return nothing
end

function ETDRK4substeps!(sol, cl, ts, eq, v, p, g)
  # Substep 1
  eq.calcN!(ts.N₁, sol, cl.t, cl, v, p, g)
  ETDRK4substep12!(ts.sol₁, ts.expLdt2, sol, ts.ζ, ts.N₁)
  @. ts.sol₁ = ts.expLdt2*sol + ts.ζ*ts.N₁
  # Substep 2
  t2 = cl.t + cl.dt/2
  eq.calcN!(ts.N₂, ts.sol₁, t2, cl, v, p, g)
  ETDRK4substep12!(ts.sol₂, ts.expLdt2, sol, ts.ζ, ts.N₂)
  # Substep 3
  eq.calcN!(ts.N₃, ts.sol₂, t2, cl, v, p, g)
  ETDRK4substep3!(ts.sol₂, ts.expLdt2, ts.sol₁, ts.ζ, ts.N₁, ts.N₃)
  # Substep 4
  t3 = cl.t + cl.dt
  eq.calcN!(ts.N₄, ts.sol₂, t3, cl, v, p, g)
  return nothing
end

function stepforward!(sol, cl, ts::ETDRK4TimeStepper, eq, v, p, g)
  ETDRK4substeps!(sol, cl, ts, eq, v, p, g)
  ETDRK4update!(sol, ts.expLdt, ts.α, ts.β, ts.Γ, ts.N₁, ts.N₂, ts.N₃, ts.N₄)
  cl.t += cl.dt
  cl.step += 1
  return nothing
end

function stepforward!(sol, cl, ts::FilteredETDRK4TimeStepper, eq, v, p, g)
  ETDRK4substeps!(sol, cl, ts, eq, v, p, g)
  ETDRK4update!(sol, ts, ts.filter) # update
  cl.t += cl.dt
  cl.step += 1
  return nothing
end


# ---
# AB3
# ---

"""
    AB3TimeStepper(eq, dev)

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

function AB3TimeStepper(eq::Equation, dev::Device=CPU())
  @devzeros typeof(dev) eq.T eq.dims RHS RHS₋₁ RHS₋₂
  AB3TimeStepper(RHS, RHS₋₁, RHS₋₂)
end


"""
    FilteredAB3TimeStepper(eq, dev; filterkwargs...)_

Construct a 3rd order Adams-Bashforth time stepper with spectral filtering.
"""
struct FilteredAB3TimeStepper{T,Tf} <: AbstractTimeStepper{T}
  RHS::T
  RHS₋₁::T
  RHS₋₂::T
  filter::Tf
end

function FilteredAB3TimeStepper(eq::Equation, dev::Device=CPU(); filterkwargs...)
  ts = AB3TimeStepper(eq, dev)
  filter = makefilter(eq; filterkwargs...)
  FilteredAB3TimeStepper(getfield.(Ref(ts), fieldnames(typeof(ts)))..., filter)
end


function AB3update!(sol, ts, cl)
  if cl.step < 3  # forward Euler steps to initialize AB3
    @. sol += cl.dt*ts.RHS    # Update
  else   # Otherwise, stepforward with 3rd order Adams Bashforth:
    @. sol += cl.dt*(ab3h1*ts.RHS - ab3h2*ts.RHS₋₁ + ab3h3*ts.RHS₋₂)
  end
  return nothing
end

function AB3update!(sol, ts::FilteredAB3TimeStepper, cl)
  if cl.step < 3  # forward Euler steps to initialize AB3
    @. sol = ts.filter*(sol + cl.dt*ts.RHS)    # Update
  else   # Otherwise, stepforward with 3rd order Adams Bashforth:
    @. sol = ts.filter*(sol + cl.dt*(ab3h1*ts.RHS - ab3h2*ts.RHS₋₁ + ab3h3*ts.RHS₋₂))
  end
  return nothing
end

function stepforward!(sol, cl, ts::AB3TimeStepper, eq, v, p, g)
  eq.calcN!(ts.RHS, sol, cl.t, cl, v, p, g)
  addlinearterm!(ts.RHS, eq.L, sol)
  AB3update!(sol, ts, cl)
  cl.t += cl.dt
  cl.step += 1
  @. ts.RHS₋₂ = ts.RHS₋₁          # Store
  @. ts.RHS₋₁ = ts.RHS            # ... previous values of RHS
  return nothing
end

function stepforward!(sol, cl, ts::FilteredAB3TimeStepper, eq, v, p, g)
  eq.calcN!(ts.RHS, sol, cl.t, cl, v, p, g)
  addlinearterm!(ts.RHS, eq.L, sol)
  AB3update!(sol, ts, cl)
  cl.t += cl.dt
  cl.step += 1
  @. ts.RHS₋₂ = ts.RHS₋₁          # Store
  @. ts.RHS₋₁ = ts.RHS            # ... previous values of RHS
  return nothing
end

# --
# Timestepper utils
# --

function getexpLs(dt, eq)
  expLdt  = @. exp(dt*eq.L)
  expLdt2 = @. exp(dt*eq.L/2)
  expLdt, expLdt2
end

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

  zc = dt*L .+ circ
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

"""
    step_until!(prob, stop_time)

Step forward `prob` up to time `stop_time`.
"""
function step_until!(prob::Problem, stop_time::Union{Float64, Float32, Float16})
  step_until!(prob::Problem, prob.timestepper::AbstractTimeStepper, stop_time)
  return nothing
end

function step_until!(prob::Problem, ::Union{ETDRK4TimeStepper, FilteredETDRK4TimeStepper}, stop_time)
  @error "step_until! requires fully explicit time stepper; does not work with ETDRK4"
  return nothing
end

function step_until!(prob::Problem, timestepper, stop_time)
  Δt = stop_time - prob.clock.t
  dt = prob.clock.dt
  nsteps = Int(floor(Δt/dt))
  stepforward!(prob, nsteps)
  t_remaining = Δt - nsteps * dt
  prob.clock.dt = t_remaining
  stepforward!(prob)
  prob.clock.dt = dt
  return nothing
end
