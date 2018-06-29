# Code Basics



## Basic Notation

The code solves partial differential equations of the general form:

$\partial_t u = \mathcal{L}u + \mathcal{N}(u)\ .$

Boundary conditions in all spatial dimensions are periodic.

Of course, ODEs are special cases of the above. Thus the code also solves ODEs.

We decompose the right hand side of the above in a linear part ($\mathcal{L}u$)
and a nonlinear part ($\mathcal{N}(u)$). The time steppers treat the linear and
nonlinear parts differently.

The coefficients for the linear operator $\mathcal{L}$ are stored in an array called `LC`. The term $\mathcal{N}(u)$ is computed for by calling the function `calcN!`.




## Abstract SuperTypes

The code is divided along conceptual lines into problem-agnostic and
problem-specific components. Files that contain problem-agnostic parts
of the code are stored in `/src`. Files in `/src` define the domain,
'AbstractTypes' that supertype problem-specific types, and
time-stepper types and routines. Problem-specific modules are stores
in `/src/physics`.

Below is a list of all Abstract Supertypes used by the code:

- `AbstractGrid`: Includes all variables that have to do with the grid, both in physical space as well as in wavenumber space. Currently implemented are: `ZeroGrid` for ODEs, `OneGrid` for PDEs with one spatial dimension, and `TwoGrid` for PDEs with two spatial dimensions. Grids are generic and work for any problem of that particular dimension.

- `AbstractParams`: Includes all parameters or functions related with the problem do not vary throughout the integration.

- `AbstractVars`: Includes all variables of the problem that change along the integration.

- `AbstractEquation`: Includes the array with the coefficients of the linear part of the equation, `LC` as well as function `calcN!` that computes the nonlinear part of the equation.

- `AbstractState`: Includes the solution `sol` at current time-step as well as the time-step `dt`, the time `t`, and `step` which counts the number of time-steps taken.

- `AbstractTimeStepper`: Includes all details for the time-stepper (e.g., `dt`, various coefficients, `sol` at intermediate time-step values). Time-steppers are generic and work for any problem.

- `AbstractProblem`: A super-supertype that includes all of the above. That is `problem` includes `grid`, `vars`, `params`, `eqn`, `ts`, `state`, and also `t` and `step`.

Grids and time-steppers are generic and work for any problem of that particular dimension. State and Problem just gathers things together. Thus, to code up a new equation you only need to prescribe `params`, `vars`, `LC`, and `calcN!`.




## Source code organization

The code is divided along conceptual lines into problem-agnostic and
problem-specific components. Files that contain problem-agnostic parts
of the code are stored in `/src`. Files in `/src` define the domain,
'AbstractTypes' that supertype problem-specific types, and
time-stepper types and routines. Problem-specific modules are stores
in `/src/physics`.

Here's an overview of the code structure:

- `/src/`
    - `FourierFlows.jl`
        - Defines supertyping AbstractParams, AbstractGrid, etc.
        - Defines a `Problem` type to organize the grid, vars, params,
            equation, and timestepper into a single structure.
        - Includes all sources files and physics files.
   - `timesteppers.jl`: defines modules and `stepforward!` routines for
        various time-steppers. Current implemented time-steppers are:
        - Forward Euler
        - 3rd-order Adams-Bashforth (AB3)
        - 4th-order Runge-Kutta (RK4)
        - 4th-order Runge-Kutta Exponential Time Differencing (ETDRK4)
        - 4th-order Dual Runge-Kutta (DualRK4)
        - 4th-order Dual Runge-Kutta Exponential Time Differencing (DualETDRK4)

        For each time-stepper exists also a "filtered" version that filters
        out high-wavenumber spectral components of the solution. The `Dual`
        time-steppers evolve a state variable that comprises both of real valued
        and complex valued fields.

    - `physics/`
        - `twodturb.jl`: Defines a `TwoDTurb` module that provides a
                solver for the two-dimensional vorticity equation.
        - `barotropicqg.jl`: Defines a `BarotropicQG` module that provides
                several solvers for the barotropic QG model that permit beta,
                topography, beta + topography, and forcing.
        - `kuramotosivashinsky.jl`: Defines a `KuramotoSivashinsky` module that
                solves the Kuramoto-Sivashinsky.
        - `verticallyfourierboussinesq.jl`: Defines a `VerticallyFourierBoussinesq` module that
                solves the two-mode truncation of the Fourier series thin-layer approximation to the hydrostatic Boussinesq equations.
        - `verticallycosinerboussinesq.jl`: Defines a `VerticallyCosineBoussinesq` module that
                solves the two-mode truncation of the Sin/Cos series thin-layer approximation to the hydrostatic Boussinesq equations.
        - `traceradvdiff.jl`: Defines a `TracerAdvDiff` module that
                provides a solver for a two-dimensional and periodic tracer
                field in a given 2D flow (u, w), which can be an arbitrary
                function of x, z, and t.

## Examples

- `examples/twodturb/McWilliams.jl`: A script that simulates decaying two-dimensional turbulence reproducing the results of the paper by

  > McWilliams, J. C. (1984). The emergence of isolated coherent vortices in turbulent flow. *J. Fluid Mech.*, **146**, 21-43

- `examples/barotropicqg/decayingbetaturb.jl`: An example script that simulates decaying quasi-geostrophic flow on a beta-plane demonstrating zonation.

- `examples/barotropicqg/ACConelayer.jl`: A script that simulates barotropic quasi-geostrophic flow above topography reproducing the results of the paper by

  > Constantinou, N. C. (2018). A barotropic model of eddy saturation. *J. Phys. Oceanogr.*, **48 (2)**, 397-411



## Tutorials

```@contents
Pages = [
    "modules/twodturb.md",
    "modules/barotropicqg.md"
    "modules/kuramotosivashinsky.md"
    "modules/traceradvdiff.md"
        ]
Depth = 1
```



## DocStrings

```@contents
Pages = [
    "general/general.md",
    "modules/twodturb.md",
    "modules/barotropicqg.md",
    "modules/boussinesq.md",
    "modules/kuramotosivashinsky.md",
    "modules/traceradvdiff.md",
    "man/docstrings.md",
    ]
Depth = 2
```

## Index

```@index
Pages = [
    "modules/twodturb.md",
    "modules/barotropicqg.md",
    "modules/boussinesq.md",
    "modules/kuramotosivashinsky.md",
    "modules/traceradvdiff.md",
    "man/docstrings.md",
    ]
```
