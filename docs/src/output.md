# Output

To save output we use [`Output`](@ref). Let's see how we can use the example developed
in [Problem](@ref problem_docs) and [Diagnostics](@ref diagnostics_docs) sections to
demonstrate how we can save some output to disk and then load it.

```@meta
DocTestSetup = quote
    using FourierFlows
    using LinearAlgebra: mul!, ldiv!
    nx, Lx = 32, 2.0
    grid = OneDGrid(; nx, Lx)
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
    E = Diagnostic(energy, prob, freq=2, nsteps=200)
    filepath = "."
    filename = joinpath(filepath, "simplestpde.jld2")
end
```

```@setup 4
using FourierFlows
using CairoMakie
CairoMakie.activate!(type = "svg")
set_theme!(Theme(linewidth = 3, fontsize = 20))
using LinearAlgebra: mul!, ldiv!

nx, Lx = 32, 2.0
grid = OneDGrid(; nx, Lx)

struct Params <: AbstractParams
  α :: Float64
end

params = Params(0.1)

struct Vars <: AbstractVars
     u :: Array{Float64, 1}
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

E = Diagnostic(energy, prob, freq=2, nsteps=200)
```

After we have created the problem (`prob`) and the energy diagnostic (`E`), we
can construct an [`Output`](@ref). We first choose where we'd like to output
the `.jld2` file(s):

```@example 4
filepath = "."
filename = joinpath(filepath, "simplestpde.jld2")
```

and the construct the [`Output`](@ref). To do so, we provide tuples of fields
that we want to output and a function what takes `prob` as its argument and
returns the corresponding value of that field. For this example, let's save
the energy `E` and the state vector `sol`.

```jldoctest; output = false, filter = r"path:.*"
get_uh(prob) = prob.sol

out = Output(prob, filename, (:uh, get_uh), (:E, energy))

# output

Output
  ├──── prob: FourierFlows.Problem{DataType, Vector{ComplexF64}, Float64, Vector{Float64}}
  ├──── path: ./simplestpde.jld2
  └── fields: Dict{Symbol, Function}(:uh => get_uh, :E => energy)
```

```@example 4
get_uh(prob) = prob.sol #hide

out = Output(prob, filename, (:uh, get_uh), (:E, energy)) #hide
```

Note that we haven't saved anything to disk yet!

By calling [`saveproblem`](@ref)

```@example 4
saveproblem(out)
```

we write all certain aspects of the problem on the `.jld2` file. For example, doing
so saves the grid parameters (`Lx`, `nx`, ...), everything in `prob.params`, the
the linear operator and its attributes from `prob.eqn`, and the time-step (`prob.clock.dt`).
All these are very useful in case we'd like to re-create the problem.

Let's open the file and have a quick look what's been written there.

```@example 4
using JLD2
file = jldopen(out.path)
```

Accessing the saved values is done using usual HDF5 way, e.g.,

```@example 4
file["grid/nx"]
```

Let's now close the `file` and move on with our demonstration.

```@example 4
close(file)
```

Now we can start saving our output fields. We do so via [`saveoutput`](@ref), which
will go through the list of fields we provided to [`Output`](@ref), call the functions
that compute the fields, and write their results on disk.

```@example 4
saveoutput(out)
```

The output fields are saved under `"snapshots"` group, e.g.,

```@example 4
file = jldopen(out.path)
file["snapshots"]
```

Let's close the `file` again.

```@example 4
close(file)
```

We can now time-step the problem forward and save the output files on disk every so often.

```@example 4
for _ in 1:40
    stepforward!(prob, E, 5)
    saveoutput(out)
end
```

All right! Now we have simulated 200 time-steps saving output every 5 time-steps.

Let's see now how we can load the output we saved in the `.jld2` file and visualize it.

We load the `.jld2` file and extract the saved iterations.

```@example 4
file = jldopen(out.path)

iterations = parse.(Int, keys(file["snapshots/t"]))
```

Then we can construct a vector with all saved times and all saved values for energy.

```@example 4
times = [file["snapshots/t/$iteration"] for iteration in iterations]

energies = [file["snapshots/E/$iteration"] for iteration in iterations]

nothing #hide
```

and plot:

```@example 4
using CairoMakie

lines(times, energies; axis = (xlabel = "t", ylabel = "energy"))

current_figure() # hide
```

Lastly, let's load the saved `uh` fields, process them (get `u` by convert to physical space),
and animate them.

We use Makie's `Observable` to animate the data. To dive into how `Observable`s work we refer to
[Makie.jl's Documentation](https://makie.juliaplots.org/stable/documentation/nodes/index.html).

```@example 4
using CairoMakie, Printf

nx = file["grid/nx"]
 x = file["grid/x"]

n = Observable(1)

u = @lift irfft(file[string("snapshots/uh/", iterations[$n])], nx)
title = @lift @sprintf("u(x, t=%1.2f)", file[string("snapshots/t/", iterations[$n])])

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "x", title)
scatterlines!(ax, x, u; markersize = 14)

frames = 1:length(iterations)

record(fig, "animation.mp4", frames, framerate=16) do i
    n[] = i
end
```

![](animation.mp4)
