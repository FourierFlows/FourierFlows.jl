# [Problem](@id problem_docs)

```@setup 2
using FourierFlows
using LinearAlgebra: mul!, ldiv!
using CairoMakie
CairoMakie.activate!(type = "svg")
set_theme!(Theme(linewidth = 3, fontsize = 20))
```

Everything needed to solve a PDE in `FourierFlows.jl` is gathered in a composite type
named [`Problem`](@ref FourierFlows.Problem). [`Problem`](@ref FourierFlows.Problem) contains
various other composite types (see [`Problem`](@ref FourierFlows.Problem) for details).

Here, we demonstrate how we can construct a [`Problem`](@ref FourierFlows.Problem) to solve
the simple 1D equation:

```math
\partial_t u(x, t) = - \alpha \, u(x, t) ,
```

on domain ``x \in [-1, 1]``.

First, we construct our grid

```@example 2
using FourierFlows

nx, Lx = 32, 2.0

grid = OneDGrid(; nx, Lx)
```

Our problem has a parameter ``\alpha``. Thus, we create a `Params` as:

```@example 2
struct Params <: AbstractParams
  α :: Float64
end
```

and then we use the `Params`'s constructor to populate our `params` with the parameter value, 
e.g., ``\alpha = 0.1``:

```@example 2
α = 0.1

params = Params(α)
```

The particular equation is so simple that it makes no difference performance-wise whether 
we time-step it in physical or in wavenumber space. For PDEs with nonlinear terms, 
time-stepping in wavenumber space is much more efficient. Thus, for demonstration purposes, 
we will time-step the equation in wavenumber space, i.e.,

```math
\partial_t \hat{u}(k, t) = - \alpha \, \hat{u}(k, t) .
```

The variables involved are ``u`` and its Fourier transform ``\hat{u}``. Thus, we 
construct the `vars` as:

```@example 2
struct Vars <: AbstractVars
    u :: Array{Float64, 1}
   uh :: Array{Complex{Float64}, 1}
end
```

and, like before, we use the `Vars`'s constructor to populate the `vars` with 
zero arrays,

```@example 2
vars = Vars(zeros(Float64, (grid.nx,)), zeros(Complex{Float64}, (grid.nkr,)))
```

Note that the Fourier transform of a real-valued array `u` is complex-valued. Also
because we use the real Fourier transform, the array `uh` is smaller.

In this simple example our state variable is simply `uh`, i.e., `sol = uh`.

Next we need to construct the equation. Equation contains the linear coefficients 
for the linear part of the PDE, stored in an array `L`, and the function `calcN!()`
that  calculates the nonlinear terms from the state variable `sol`. In our case,
our equation is linear and, therefore,

```@example 2
L = - params.α * ones(grid.nkr)

nothing # hide
```

and

```@example 2
function calcN!(N, sol, t, clock, vars, params, grid)
  @. N = 0
  
  return nothing
end

nothing # hide
```

Note that `calcN!()` needs to have the above argument structure. With `L` and `calcN!`
in hand we can construct our problem's equation:

```@example 2
equation = FourierFlows.Equation(L, calcN!, grid)
```

Last, we have to pick a time-stepper and a time-step `dt` and gather everything 
a FourierFlows's [`Problem`](@ref FourierFlows.Problem). Time-steppers are prescribed
via a string that corresponds to the name of the implemented time-steppers _without_
the `TimeStepper` ending. (See [Time-stepping section](@ref timestepping) for a list
of implemented time-stepping schemes.) For example, a problem that uses a
[`ForwardEulerTimeStepper`](@ref) with time step of 0.02 is constructed via:

```@example 2
stepper, dt = "ForwardEuler", 0.02

prob = FourierFlows.Problem(equation, stepper, dt, grid, vars, params)
```

By default, the `Problem` constructor takes `sol` a complex valued array filled with zeros
with same size as `L`.

The [`prob.clock`](@ref FourierFlows.Clock) contains the time-step `dt` and the current `step`
and time `t` of the simulation:

```@example 2
prob.clock
```

Let's initiate our problem with, e.g., ``u(x, 0) = \cos(\pi x)``, integrate up 
to ``t = 4`` and compare our numerical solution with the analytic solution 
``u(x, t) = e^{-\alpha t} \cos(\pi x)``.

```@example 2
u0 = @. cos(π * grid.x)

using LinearAlgebra: mul!

mul!(prob.sol, grid.rfftplan, u0)

nothing # hide
```

Since our time-step is chosen `dt = 0.02`, we need to step forward `prob` for ``200`` 
time-steps to reach ``t = 4``.

```@example 2
stepforward!(prob, 200)
```

Now let's transform our state vector `sol` back in physical space

```@example 2
using LinearAlgebra: ldiv!

ldiv!(prob.vars.u, grid.rfftplan, prob.sol)

nothing # hide
```

and finally, let's plot our solution and compare with the analytic solution:

```@example 2
using CairoMakie, Printf

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "x", title = @sprintf("u(x, t=%1.2f)", prob.clock.t))

x = range(-Lx/2, Lx/2, length=200)
lines!(ax, x, @. cos(π * x) * exp(-prob.params.α * 4); label = "analytical")
lines!(ax, x, @. cos(π * x); linestyle = :dash, color = :gray, label = "initial condition")

scatter!(ax, grid.x, prob.vars.u; markersize = 14, color = :salmon, label = "numerical")

axislegend()

current_figure() # hide
```

A good practice is to encompass all functions and type definitions related with a PDE under 
a single module, e.g.,

```julia
module MyPDE

  ...

end # end module
```

For a more elaborate example we urge you to have a look at the [`Diffusion`](@ref FourierFlows.Diffusion) 
module located at [`src/diffusion.jl`](https://github.com/FourierFlows/FourierFlows.jl/blob/main/src/diffusion.jl)
and also the modules included in the child package
[GeophysicalFlows.jl](https://github.com/FourierFlows/GeophysicalFlows.jl).
