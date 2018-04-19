# FourierFlows.jl Documentation


## Overview

FourierFlows provides solvers for partial differential equations on
doubly-periodic domains using Fourier-based pseudospectral methods.
A central intent of the software's design is also to provide a framework
for writing new, fast solvers for new physical problems.
The code is written in [Julia](https://julialang.org).


## Installation

FourierFlows is a registered package and so can be installed via `Pkg.add`.

```julia
Pkg.add("FourierFlows")
```

For now, this package supports Julia `0.6`. Support for version `0.7` is on its
way.


## Basic Notation

The code solves partial differential equations of the general form:

$\partial_t u = \mathcal{L}u + \mathcal{N}(u)\ .$

We decompose the right hand side of the above in a linear part ($\mathcal{L}u$)
and a nonlinear part ($\mathcal{N}(u)$). The time steppers treat the linear and
nonlinear parts differently.

The coefficients for the linear operator $\mathcal{L}$ are stored in array `LC`.
The term $\mathcal{N}(u)$ is computed for by calling the function `calcN!`.


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

## Writing fast solvers

The performance-intensive part of the code involves just two functions: the
time-stepping scheme `stepforward!`, and the function `calcN!` that
calculates the nonlinear part of the given equation's right-hand side.
Optimization of these two functions for a given problem will produce the
fastest possible code.


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


## Future work

The code is in the chaotic stage of development. A main goal for the future
is to permit the use of shared memory parallelism in the compute-intensive
routines (shared-memory parallelism provided already by FFTW/MKLFFT, but
is not yet native to Julia for things like element-wise matrix multiplication,
addition, and assignment). This feature may possibly be enabled by
Intel Lab's [ParallelAccelerator](https://github.com/IntelLabs/ParallelAccelerator.jl)
package.


## Developers

FourierFlows is currently being developed by [Gregory L. Wagner](https://glwagner.github.io) and
[Navid C. Constantinou](http://www.navidconstantinou.com).


## DocStrings

```@contents
Pages = [
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
