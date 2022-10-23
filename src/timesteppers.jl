"""
    stepforward!(prob)

Step forward `prob` one time step.
"""
stepforward!(prob::Problem) =
  stepforward!(prob.sol, prob.clock, prob.timestepper, prob.eqn, prob.vars, prob.params, prob.grid)

"""
    stepforward!(prob, nsteps::Int)

Step forward `prob` for `nsteps`.
"""
function stepforward!(prob::Problem, nsteps::Int)
  for _ in 1:nsteps
    stepforward!(prob)
  end
  
  return nothing
end

"""
    stepforward!(prob::Problem, diags, nsteps::Int)

Step forward `prob` for `nsteps`, incrementing `diags` along the way. `diags` may be a 
single `Diagnostic` or a `Vector` of `Diagnostic`s.
"""
function stepforward!(prob::Problem, diags, nsteps::Int)
  for _ in 1:nsteps
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
  :FilteredAB3,
  :LSRK54
]

isexplicit(stepper) = any(Symbol(stepper) .== fullyexplicitsteppers)

"""
    TimeStepper(stepper, equation, dt=nothing, dev=CPU(); kw...)

Instantiate the Time`stepper` for `equation` with timestep `dt` and
on the `dev`ice. The `kw` are passed to the timestepper constructor.
"""
function TimeStepper(stepper, equation, dt=nothing, dev::Device=CPU(); kw...)
  fullsteppername = Symbol(stepper, :TimeStepper)

  # Create expression that instantiates the time-stepper, depending on whether 
  # timestepper is explicit or not.
  expr = isexplicit(stepper) ? Expr(:call, fullsteppername, equation, dev) :
                               Expr(:call, fullsteppername, equation, dt, dev)

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
#   * LSRK54
#   * Filtered ETDRK4
#   * AB3
#   * Filtered AB3
#
# Explicit time-steppers are constructed with the signature
#   ts = ExplicitTimeStepper(equation::Equation)
#
# Implicit time-steppers are constructed with the signature
#   ts = ImplicitTimeStepper(equation::Equation, dt)

# --
# Forward Euler
# --

"""
    struct ForwardEulerTimeStepper{T} <: AbstractTimeStepper{T}

A Forward Euler timestepper for time-stepping `∂u/∂t = RHS(u, t)` via:
```
uⁿ⁺¹ = uⁿ + dt * RHS(uⁿ, tⁿ)
```
"""
struct ForwardEulerTimeStepper{T} <: AbstractTimeStepper{T}
  N::T # Explicit linear and nonlinear terms
  ForwardEulerTimeStepper(N::T) where T = new{T}(0N)
end

"""
    ForwardEulerTimeStepper(equation::Equation, dev::Device=CPU())

Construct a Forward Euler timestepper for `equation` on device `dev`.
"""
ForwardEulerTimeStepper(equation::Equation, dev::Device=CPU()) = 
  ForwardEulerTimeStepper(zeros(dev, equation.T, equation.dims))

function stepforward!(sol, clock, ts::ForwardEulerTimeStepper, equation, vars, params, grid)
  equation.calcN!(ts.N, sol, clock.t, clock, vars, params, grid)
  @. sol += clock.dt * (equation.L * sol + ts.N)
  clock.t += clock.dt
  clock.step += 1
  
  return nothing
end

"""
    struct FilteredForwardEulerTimeStepper{T,Tf} <: AbstractTimeStepper{T}

A Forward Euler timestepper with spectral filtering. See [`ForwardEulerTimeStepper`](@ref).
"""
struct FilteredForwardEulerTimeStepper{T,Tf} <: AbstractTimeStepper{T}
       N :: T
  filter :: Tf
end

"""
    FilteredForwardEulerTimeStepper(equation, dev; filterkwargs...)

Construct a Forward Euler timestepper with spectral filtering for `equation` on device `dev`.
"""
function FilteredForwardEulerTimeStepper(equation::Equation, dev::Device=CPU(); filterkwargs...)
  filter = makefilter(equation; filterkwargs...)
  
  return FilteredForwardEulerTimeStepper(zeros(dev, equation.T, equation.dims), filter)
end

function stepforward!(sol, clock, ts::FilteredForwardEulerTimeStepper, equation, vars, params, grid)
  equation.calcN!(ts.N, sol, clock.t, clock, vars, params, grid)
  @. sol = ts.filter * (sol + clock.dt * (ts.N + equation.L * sol))
  clock.t += clock.dt
  clock.step += 1
  
  return nothing
end


# --
# RK4
# --

"""
    struct RK4TimeStepper{T} <: AbstractTimeStepper{T}

A 4th-order Runge-Kutta timestepper for time-stepping `∂u/∂t = RHS(u, t)` via:
```
uⁿ⁺¹ = uⁿ + dt/6 * (k₁ + 2 * k₂ + 2 * k₃ + k₄)
```
where
```
k₁ = RHS(uⁿ, tⁿ)
k₂ = RHS(uⁿ + k₁ * dt/2, tⁿ + dt/2)
k₃ = RHS(uⁿ + k₂ * dt/2, tⁿ + dt/2)
k₄ = RHS(uⁿ + k₃ * dt, tⁿ + dt)
```
"""
struct RK4TimeStepper{T} <: AbstractTimeStepper{T}
  sol₁ :: T
  RHS₁ :: T
  RHS₂ :: T
  RHS₃ :: T
  RHS₄ :: T
end

"""
    RK4TimeStepper(equation::Equation, dev::Device=CPU())

Construct a 4th-order Runge-Kutta timestepper for `equation` on device `dev`.
"""
function RK4TimeStepper(equation::Equation, dev::Device=CPU())
  @devzeros typeof(dev) equation.T equation.dims N sol₁ RHS₁ RHS₂ RHS₃ RHS₄
  
  return RK4TimeStepper(sol₁, RHS₁, RHS₂, RHS₃, RHS₄)
end

"""
    struct FilteredRK4TimeStepper{T,Tf} <: AbstractTimeStepper{T}

A 4th-order Runge-Kutta timestepper with spectral filtering. See [`RK4TimeStepper`](@ref).
"""
struct FilteredRK4TimeStepper{T,Tf} <: AbstractTimeStepper{T}
    sol₁ :: T
    RHS₁ :: T
    RHS₂ :: T
    RHS₃ :: T
    RHS₄ :: T
  filter :: Tf
end

"""
    FilteredRK4TimeStepper(equation::Equation, dev::Device=CPU(); filterkwargs...)

Construct a 4th-order Runge-Kutta timestepper with spectral filtering for `equation` on device `dev`.
"""
function FilteredRK4TimeStepper(equation::Equation, dev::Device=CPU(); filterkwargs...)
  ts = RK4TimeStepper(equation, dev)
  filter = makefilter(equation; filterkwargs...)
  
  return FilteredRK4TimeStepper(getfield.(Ref(ts), fieldnames(typeof(ts)))..., filter)
end

function addlinearterm!(RHS, L, sol)
  @. RHS += L*sol
  
  return nothing
end

function substepsol!(newsol, sol, RHS, dt)
  @. newsol = sol + dt*RHS
  
  return nothing
end

function RK4substeps!(sol, clock, ts, equation, vars, params, grid, t, dt)
  # Substep 1
  equation.calcN!(ts.RHS₁, sol, t, clock, vars, params, grid)
  addlinearterm!(ts.RHS₁, equation.L, sol)
  
  # Substep 2
  substepsol!(ts.sol₁, sol, ts.RHS₁, dt/2)
  equation.calcN!(ts.RHS₂, ts.sol₁, t+dt/2, clock, vars, params, grid)
  addlinearterm!(ts.RHS₂, equation.L, ts.sol₁)
  
  # Substep 3
  substepsol!(ts.sol₁, sol, ts.RHS₂, dt/2)
  equation.calcN!(ts.RHS₃, ts.sol₁, t+dt/2, clock, vars, params, grid)
  addlinearterm!(ts.RHS₃, equation.L, ts.sol₁)
  
  # Substep 4
  substepsol!(ts.sol₁, sol, ts.RHS₃, dt)
  equation.calcN!(ts.RHS₄, ts.sol₁, t+dt, clock, vars, params, grid)
  addlinearterm!(ts.RHS₄, equation.L, ts.sol₁)
  
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

function stepforward!(sol, clock, ts::RK4TimeStepper, equation, vars, params, grid)
  RK4substeps!(sol, clock, ts, equation, vars, params, grid, clock.t, clock.dt)
  RK4update!(sol, ts.RHS₁, ts.RHS₂, ts.RHS₃, ts.RHS₄, clock.dt)
  clock.t += clock.dt
  clock.step += 1
  
  return nothing
end

function stepforward!(sol, clock, ts::FilteredRK4TimeStepper, equation, vars, params, grid)
  RK4substeps!(sol, clock, ts, equation, vars, params, grid, clock.t, clock.dt)
  RK4update!(sol, ts.RHS₁, ts.RHS₂, ts.RHS₃, ts.RHS₄, ts.filter, clock.dt)
  clock.t += clock.dt
  clock.step += 1
  
  return nothing
end

# ------
# LSRK(5)4
# ------
"""
    struct LSRK54TimeStepper{T} <: AbstractTimeStepper{T}

A 4th-order 5-stages 2-storage Runge-Kutta timestepper for time-stepping
`∂u/∂t = RHS(u, t)` via:
```
S² = 0
for i = 1:5
  S²   = Aᵢ * S² + dt * RHS(uⁿ, t₀ + Cᵢ * dt)
  uⁿ⁺¹ = uⁿ + Bᵢ * S²
end
```

where `Aᵢ`, `Bᵢ`, and `Cᵢ` are the A, B, C coefficients from LSRK tableau
table at the ``i``-th stage. For details, please refer to

> Carpenter, M. H. and Kennedy, C. A. (1994). Fourth-order 2N-storage Runge–Kutta schemes,
  Technical Report NASA TM-109112, NASA Langley Research Center, VA.
"""
struct LSRK54TimeStepper{T,V} <: AbstractTimeStepper{T}
   S² :: T
  RHS :: T
    A :: V
    B :: V
    C :: V
end

"""
    LSRK54TimeStepper(equation::Equation, dev::Device=CPU())

Construct a 4th-order 5-stages low storage Runge-Kutta timestepper for `equation` on device `dev`.
"""
function LSRK54TimeStepper(equation::Equation, dev::Device=CPU())
  @devzeros typeof(dev) equation.T equation.dims S² RHS
  
  T = equation.T
  A = T[0.0,               -0.417890474499852, -1.19215169464268 , -1.69778469247153 , -1.51418344425716 ]
  B = T[0.149659021999229,  0.379210312999627,  0.822955029386982,  0.699450455949122,  0.153057247968152]
  C = T[0.0,                0.149659021999229,  0.3704009573642045, 0.6222557631344415, 0.95828213067469 ]

  return LSRK54TimeStepper(S², Fⁱ, Tuple(A), Tuple(B), Tuple(C))
end

function LSRK54!(sol, clock, ts, equation, vars, params, grid, t, dt)
    T = equation.T
    A, B, C = ts.A, ts.B, ts.C

    # initialize the S² term
    @. ts.S² = 0

    for i = 1:5
      equation.calcN!(ts.RHS, sol, t + C[i] * dt , clock, vars, params, grid)
      addlinearterm!(ts.RHS, equation.L, sol)

      @. ts.S² = A[i] * ts.S² + dt * ts.RHS
      @.  sol += B[i] * ts.S²
    end

    return nothing
end

function stepforward!(sol, clock, ts::LSRK54TimeStepper, equation, vars, params, grid)
  LSRK54!(sol, clock, ts, equation, vars, params, grid, clock.t, clock.dt)
  clock.t    += clock.dt
  clock.step += 1
  return nothing
end

# ------
# ETDRK4
# ------

"""
    struct ETDRK4TimeStepper{T,TL} <: AbstractTimeStepper{T}

A 4th-order exponential-time-differencing Runge-Kutta timestepper for time-stepping
`∂u/∂t = L * u + N(u)`. The scheme treats the linear term `L` exact while for the
nonlinear terms `N(u)` it uses a 4th-order Runge-Kutta scheme. That is,
```
uⁿ⁺¹ = exp(L * dt) * uⁿ + RK4(N(uⁿ))
```
For more info refer to 

> Kassam, A. K., & Trefethen, L. N. (2005). Fourth-order time-stepping for stiff PDEs. _SIAM Journal on Scientific Computing_, **26(4)**, 1214-1233.
"""
struct ETDRK4TimeStepper{T,TL} <: AbstractTimeStepper{T}
  # ETDRK4 coefficents
        ζ :: TL
        α :: TL
        β :: TL
        Γ :: TL
  expLdt  :: TL
  exp½Ldt :: TL
     sol₁ :: T
     sol₂ :: T
       N₁ :: T
       N₂ :: T
       N₃ :: T
       N₄ :: T
end

"""
    ETDRK4TimeStepper(equation::Equation, dt, dev::Device=CPU())

Construct a 4th-order exponential-time-differencing Runge-Kutta timestepper with timestep `dt`
for `equation` on device `dev`.
"""
function ETDRK4TimeStepper(equation::Equation, dt, dev::Device=CPU())
  dt = fltype(equation.T)(dt) # ensure dt is correct type.
  expLdt, exp½Ldt = getexpLs(dt, equation)
  ζ, α, β, Γ = getetdcoeffs(dt, equation.L)
  @devzeros typeof(dev) equation.T equation.dims sol₁ sol₂ N₁ N₂ N₃ N₄
  
  return ETDRK4TimeStepper(ζ, α, β, Γ, expLdt, exp½Ldt, sol₁, sol₂, N₁, N₂, N₃, N₄)
end

"""
    FilteredETDRK4TimeStepper{T,TL,Tf} <: AbstractTimeStepper{T}

A 4th-order exponential-time-differencing Runge-Kutta timestepper with spectral filtering.
See [`ETDRK4TimeStepper`](@ref).
"""
struct FilteredETDRK4TimeStepper{T,TL,Tf} <: AbstractTimeStepper{T}
  # ETDRK4 coefficents:
        ζ :: TL
        α :: TL
        β :: TL
        Γ :: TL
  expLdt  :: TL
  exp½Ldt :: TL
     sol₁ :: T
     sol₂ :: T
       N₁ :: T
       N₂ :: T
       N₃ :: T
       N₄ :: T
   filter :: Tf
end

"""
    FilteredETDRK4TimeStepper(equation, dt; filterkwargs...)

Construct a 4th-order exponential-time-differencing Runge-Kutta timestepper with timestep `dt` and with
spectral filtering for `equation` on device `dev`.
"""
function FilteredETDRK4TimeStepper(equation::Equation, dt, dev::Device=CPU(); filterkwargs...)
  timestepper = ETDRK4TimeStepper(equation, dt, dev)
  filter = makefilter(equation; filterkwargs...)
  
  return FilteredETDRK4TimeStepper(getfield.(Ref(timestepper), fieldnames(typeof(timestepper)))..., filter)
end

function ETDRK4update!(sol, expLdt, α, β, Γ, N₁, N₂, N₃, N₄)
  @. sol = (expLdt * sol +  α * N₁
                         + 2β * (N₂ + N₃)
                         +  Γ * N₄ )

  return nothing
end

function ETDRK4update!(sol, ts, filter)
  @. sol = filter * (ts.expLdt * sol +     ts.α * ts.N₁
                                     + 2 * ts.β * (ts.N₂ + ts.N₃)
                                     +     ts.Γ * ts.N₄ )

  return nothing
end

function ETDRK4substep12!(sol₁, exp½Ldt, sol, ζ, N)
  @. sol₁ = exp½Ldt * sol + ζ * N

  return nothing
end

function ETDRK4substep3!(sol₂, exp½Ldt, sol₁, ζ, N₁, N₃)
  @. sol₂ = exp½Ldt * sol₁ + ζ * (2N₃ - N₁)

  return nothing
end

function ETDRK4substeps!(sol, clock, ts, equation, vars, params, grid)
  # Substep 1
  equation.calcN!(ts.N₁, sol, clock.t, clock, vars, params, grid)
  ETDRK4substep12!(ts.sol₁, ts.exp½Ldt, sol, ts.ζ, ts.N₁)

  # Substep 2
  t2 = clock.t + clock.dt/2
  equation.calcN!(ts.N₂, ts.sol₁, t2, clock, vars, params, grid)
  ETDRK4substep12!(ts.sol₂, ts.exp½Ldt, sol, ts.ζ, ts.N₂)

  # Substep 3
  equation.calcN!(ts.N₃, ts.sol₂, t2, clock, vars, params, grid)
  ETDRK4substep3!(ts.sol₂, ts.exp½Ldt, ts.sol₁, ts.ζ, ts.N₁, ts.N₃)

  # Substep 4
  t3 = clock.t + clock.dt
  equation.calcN!(ts.N₄, ts.sol₂, t3, clock, vars, params, grid)

  return nothing
end

function stepforward!(sol, clock, ts::ETDRK4TimeStepper, equation, vars, params, grid)
  ETDRK4substeps!(sol, clock, ts, equation, vars, params, grid)
  ETDRK4update!(sol, ts.expLdt, ts.α, ts.β, ts.Γ, ts.N₁, ts.N₂, ts.N₃, ts.N₄)
  clock.t += clock.dt
  clock.step += 1

  return nothing
end

function stepforward!(sol, clock, ts::FilteredETDRK4TimeStepper, equation, vars, params, grid)
  ETDRK4substeps!(sol, clock, ts, equation, vars, params, grid)
  ETDRK4update!(sol, ts, ts.filter) # update
  clock.t += clock.dt
  clock.step += 1

  return nothing
end


# ---
# AB3
# ---

const ab3h1 = 23/12
const ab3h2 = 16/12
const ab3h3 = 5/12

"""
    struct AB3TimeStepper{T} <: AbstractTimeStepper{T}

A 3rd-order Adams-Bashforth timestepper for time-stepping `∂u/∂t = RHS(u, t)` via:
```
uⁿ⁺¹ = uⁿ + dt/12 * (23 * RHS(uⁿ, tⁿ) - 16 * RHS(uⁿ⁻¹, tⁿ⁻¹) + 5 * RHS(uⁿ⁻², tⁿ⁻²))
```

Adams-Bashforth is a multistep method, i.e., it not only requires information from the `n`-th time-step
(`uⁿ`) but also from the previous two timesteps (`uⁿ⁻¹` and `uⁿ⁻²`). For the first two timesteps, it
falls back to a forward Euler timestepping scheme:
```
uⁿ⁺¹ = uⁿ + dt * RHS(uⁿ, tⁿ)
```
"""
struct AB3TimeStepper{T} <: AbstractTimeStepper{T}
  RHS::T
  RHS₋₁::T
  RHS₋₂::T
end

"""
    AB3TimeStepper(equation::Equation, dev::Device=CPU())

Construct a 3rd order Adams-Bashforth timestepper for `equation` on device `dev`.
"""
function AB3TimeStepper(equation::Equation, dev::Device=CPU())
  @devzeros typeof(dev) equation.T equation.dims RHS RHS₋₁ RHS₋₂

  return AB3TimeStepper(RHS, RHS₋₁, RHS₋₂)
end


"""
    struct FilteredAB3TimeStepper{T} <: AbstractTimeStepper{T}

A 3rd order Adams-Bashforth timestepper with spectral filtering. See [`AB3TimeStepper`](@ref).
"""
struct FilteredAB3TimeStepper{T, Tf} <: AbstractTimeStepper{T}
   RHS   :: T
   RHS₋₁ :: T
   RHS₋₂ :: T
  filter :: Tf
end

"""
    FilteredAB3TimeStepper(equation::Equation, dev::Device=CPU(); filterkwargs...)

Construct a 3rd order Adams-Bashforth timestepper with spectral filtering for `equation` on device `dev`.
"""
function FilteredAB3TimeStepper(equation::Equation, dev::Device=CPU(); filterkwargs...)
  timestepper = AB3TimeStepper(equation, dev)
  filter = makefilter(equation; filterkwargs...)
  
  return FilteredAB3TimeStepper(getfield.(Ref(timestepper), fieldnames(typeof(timestepper)))..., filter)
end

function AB3update!(sol, ts, clock)
  if clock.step < 3  # forward Euler steps to initialize AB3
    @. sol += clock.dt * ts.RHS    # Update
  else   # Otherwise, stepforward with 3rd order Adams Bashforth:
    @. sol += clock.dt * (ab3h1 * ts.RHS - ab3h2 * ts.RHS₋₁ + ab3h3 * ts.RHS₋₂)
  end
  
  return nothing
end

function AB3update!(sol, ts::FilteredAB3TimeStepper, clock)
  if clock.step < 3  # forward Euler steps to initialize AB3
    @. sol = ts.filter * (sol + clock.dt * ts.RHS)    # Update
  else   # Otherwise, stepforward with 3rd order Adams Bashforth:
    @. sol = ts.filter * (sol + clock.dt * (ab3h1 * ts.RHS - ab3h2 * ts.RHS₋₁ + ab3h3 * ts.RHS₋₂))
  end
  
  return nothing
end

function stepforward!(sol, clock, ts::AB3TimeStepper, equation, vars, params, grid)
  equation.calcN!(ts.RHS, sol, clock.t, clock, vars, params, grid)
  addlinearterm!(ts.RHS, equation.L, sol)
  
  AB3update!(sol, ts, clock)
  clock.t += clock.dt
  clock.step += 1
  
  @. ts.RHS₋₂ = ts.RHS₋₁          # Store
  @. ts.RHS₋₁ = ts.RHS            # ... previous values of RHS
  
  return nothing
end

function stepforward!(sol, clock, ts::FilteredAB3TimeStepper, equation, vars, params, grid)
  equation.calcN!(ts.RHS, sol, clock.t, clock, vars, params, grid)
  addlinearterm!(ts.RHS, equation.L, sol)
  
  AB3update!(sol, ts, clock)
  clock.t += clock.dt
  clock.step += 1
  
  @. ts.RHS₋₂ = ts.RHS₋₁          # Store
  @. ts.RHS₋₁ = ts.RHS            # ... previous values of RHS
  
  return nothing
end

# --
# Timestepper utils
# --

function getexpLs(dt, equation)
  expLdt  = @. exp(dt * equation.L)
  exp½Ldt = @. exp(dt * equation.L/2)
  
  return expLdt, exp½Ldt
end

"""
    getetdcoeffs(dt, L; ncirc=32, rcirc=1)

Calculate ETDRK4 coefficients associated with the (diagonal) linear coefficient
`L` by integrating over a unit circle in the complex space.
"""
function getetdcoeffs(dt, L; ncirc=32, rcirc=1)
  shape = Tuple(cat(ncirc, ones(Int, ndims(L)), dims=1))

  circ = zeros(Complex{Float64}, shape) # use double precision for this calculation
  circ .= rcirc * exp.(2π * im/ncirc * (0.5:1:(ncirc-0.5)))
  circ = permutedims(circ, ndims(circ):-1:1)

  zc = dt * L .+ circ
  M = ndims(L) + 1

  # Four coefficients: ζ, α, β, Γ
  ζc = @. ( exp(zc/2)-1 ) / zc
  αc = @. ( -4 - zc + exp(zc) * (4 - 3zc + zc^2) ) / zc^3
  βc = @. ( 2  + zc + exp(zc) * (-2 + zc) ) / zc^3
  Γc = @. ( -4 - 3zc - zc^2 + exp(zc) * (4 - zc) ) / zc^3

  ζ = dt * dropdims(mean(ζc, dims=M), dims=M)
  α = dt * dropdims(mean(αc, dims=M), dims=M)
  β = dt * dropdims(mean(βc, dims=M), dims=M)
  Γ = dt * dropdims(mean(Γc, dims=M), dims=M)

  if eltype(L) <: Real # this is conservative, but unclear if necessary
    ζ = real.(ζ)
    α = real.(α)
    β = real.(β)
    Γ = real.(Γ)
  end

  return ζ, α, β, Γ
end

getetdcoeffs(dt, L::CuArray; kwargs...) = 
    (CuArray(ζ) for ζ in getetdcoeffs(dt, Array(L); kwargs...))

"""
    step_until!(prob, stop_time)

Step forward `prob` until `stop_time`. Cannot be used with ETDRK4 time steppers.

See also: [`stepforward!`](@ref)
"""
step_until!(prob, stop_time) = step_until!(prob, prob.timestepper, stop_time)

step_until!(prob, ::Union{ETDRK4TimeStepper, FilteredETDRK4TimeStepper}, stop_time) =
  error("step_until! requires fully explicit time stepper; does not work with ETDRK4")

function step_until!(prob, timestepper, stop_time)
  # Throw an error if stop_time is not greater than the current problem time
  stop_time > prob.clock.t || error("stop_time must be greater than prob.clock.t")

  # Extract current time step
  dt = prob.clock.dt

  # Step forward until just before stop_time
  time_interval = stop_time - prob.clock.t
  nsteps = floor(Int, time_interval / dt)
  stepforward!(prob, nsteps)

  # Take one final small step so that prob.clock.t = stop_time
  t_remaining = time_interval - prob.clock.t
  prob.clock.dt = t_remaining
  stepforward!(prob)

  # Restore previous time-step
  prob.clock.dt = dt

  return nothing
end
