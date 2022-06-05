# [Diagnostics](@id diagnostics_docs)

We can add diagnostics to a FourierFlows's problem using [`Diagnostic`](@ref)
functionality.

```@meta
DocTestSetup = quote
    using FourierFlows
    using LinearAlgebra: mul!, ldiv!
    nx, Lx = 32, 2.0
    grid = OneDGrid(nx, Lx)
    struct Params <: AbstractParams
    α :: Float64
    end
    params = Params(0.1)
    struct Vars <: AbstractVars
        u :: Array{Float64,1}
    uh :: Array{Complex{Float64}, 1}
    end
    vars = Vars(zeros(Float64, (grid.nx,)), zeros(Complex{Float64}, (grid.nkr,)))
    L = - params.α * ones(grid.nkr)
    function calcN!(N, sol, t, clock, vars, params, grid)
    @. N = 0
    return nothing
    end
    equation = FourierFlows.Equation(L, calcN!, grid)
    stepper, dt = "ForwardEuler", 0.02
    prob = FourierFlows.Problem(equation, stepper, dt, grid, vars, params)
    u0 = @. cos(π * grid.x)
    mul!(prob.sol, grid.rfftplan, u0)
    function energy(prob)
        ldiv!(prob.vars.u, grid.rfftplan, prob.sol)
        return sum(prob.vars.u.^2) * prob.grid.dx
    end
end
```

```@setup 3
using FourierFlows
using CairoMakie
set_theme!(Theme(linewidth = 3, fontsize = 20))

using LinearAlgebra: mul!

nx, Lx = 32, 2.0
grid = OneDGrid(nx, Lx)

struct Params <: AbstractParams
  α :: Float64
end

params = Params(0.1)

struct Vars <: AbstractVars
    u :: Array{Float64,1}
   uh :: Array{Complex{Float64}, 1}
end

vars = Vars(zeros(Float64, (grid.nx,)), zeros(Complex{Float64}, (grid.nkr,)))

L = - params.α * ones(grid.nkr)

function calcN!(N, sol, t, clock, vars, params, grid)
  @. N = 0
  
  return nothing
end

equation = FourierFlows.Equation(L, calcN!, grid)

stepper, dt = "ForwardEuler", 0.02

prob = FourierFlows.Problem(equation, stepper, dt, grid, vars, params)

u0 = @. cos(π * grid.x)

mul!(prob.sol, grid.rfftplan, u0)
```

To demonstrate how we add diagnostics to a PDE problem, let's try to add
one to the simple PDE problem we constructed in the [Problem](@ref problem_docs)
section. For example, say we'd like to add a diagnostic we refer to as the "energy"
and which we define to be:

```math
E = \int u^2 \, \mathrm{d} x .
```

After we have constructed the problem (`prob`) (see [Problem](@ref problem_docs) section),
we then create a function that takes `prob` as its argument returns the diagnostic:

```@example 3
using LinearAlgebra: ldiv!

function energy(prob)
    ldiv!(prob.vars.u, grid.rfftplan, prob.sol)

    return sum(prob.vars.u.^2) * prob.grid.dx
end

nothing #hide
```

and then we create a [`Diagnostic`](@ref) using the
[`Diagnostic`](@ref Diagnostic(calc, prob; freq=1, nsteps=100, ndata=ceil(Int, (nsteps+1)/freq)))
constructor. Say we want to save energy every 2 time-steps, then:

```jldoctest; output = false
E = Diagnostic(energy, prob, freq=2, nsteps=200)

# output

Diagnostic
  ├─── calc: energy
  ├─── prob: FourierFlows.Problem{DataType, Vector{ComplexF64}, Float64, Vector{Float64}}
  ├─── data: 101-element Vector{Float64}
  ├────── t: 101-element Vector{Float64}
  ├── steps: 101-element Vector{Int64}
  ├─── freq: 2
  └────── i: 1
```
```@example 3
E = Diagnostic(energy, prob, freq=2, nsteps=200) #hide
```

Now, when we step forward the problem we provide the diagnostic as the second positional 
argument in [`stepforward!`](@ref stepforward!(prob, diags, nsteps)):

```@example 3
stepforward!(prob, E, 200)
```

Doing so, the diagnostic is computed and saved at the appropriate frequency (prescribed
by `E.freq`).

!!! info "Multiple diagnostics"
    If we want to include multiple diagnostics we can gather all of them in an array, e.g.,
    ```julia
    diag1 = Diagnostic(foo, prob)
    diag2 = Diagnostic(bar, prob)

    stepforward!(prob, [diag1, diag2], 1)
    ```

The times that the diagnostic was saved are gathered in `E.t`. Thus, we can easily
plot the energy time-series, e.g., 

```@example 3
using CairoMakie

lines(E.t, E.data, axis = (xlabel = "time", ylabel = "energy"))

current_figure() # hide
```
