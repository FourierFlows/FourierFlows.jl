[![Build Status](https://travis-ci.org/FourierFlows/FourierFlows.jl.svg?branch=master)](https://travis-ci.org/FourierFlows/FourierFlows.jl) [![Build status](https://ci.appveyor.com/api/projects/status/3hm86k8d4qdch730?svg=true)](https://ci.appveyor.com/project/navidcy/fourierflows-jl) [![codecov](https://codecov.io/gh/FourierFlows/FourierFlows.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/FourierFlows/FourierFlows.jl) [![Coverage Status](https://coveralls.io/repos/github/FourierFlows/FourierFlows.jl/badge.svg?branch=master)](https://coveralls.io/github/FourierFlows/FourierFlows.jl?branch=master)

[![FourierFlows](http://pkg.julialang.org/badges/FourierFlows_0.6.svg)](http://pkg.julialang.org/detail/FourierFlows)
[![FourierFlows](http://pkg.julialang.org/badges/FourierFlows_0.7.svg)](http://pkg.julialang.org/detail/FourierFlows)


[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://FourierFlows.github.io/FourierFlows.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://FourierFlows.github.io/FourierFlows.jl/latest)

[![DOI](https://zenodo.org/badge/95982414.svg)](https://zenodo.org/badge/latestdoi/95982414)

# FourierFlows

## Overview

This software provides solvers for partial differential equations on
doubly-periodic domains using Fourier-based pseudospectral methods.
A central intent of the software's design is also to provide a framework
for writing new, fast solvers for new physical problems.
The code is written in [Julia][].

For more details refer to the [documentation](https://fourierflows.github.io/FourierFlows.jl/latest/).


## Source code organization

The code is divided along conceptual lines into problem-agnostic and
problem-specific components. Files that contain problem-agnostic parts
of the code are stored in `/src`. Files in `/src` define the domain,
'AbstractTypes' that supertype problem-specific types, and
time-stepper types and routines. Problem-specific modules are stores in
`/src/physics`.

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
        - 4th-order Runge-Kutta Exponential Time Differencing

        For each time-stepper exists also a "filtered" version that filters
        out high-wavenumber spectral components of the solution.
    - `physics/`
        - `kuramotosivashinsky.jl`: Defines a `KuramotoSivashinsky` module that
                provides the one-dimensional Kuramoto-Sivashinsky equation.
        - `twodturb.jl`: Defines a `TwoDTurb` module that provides a
                solver for the two-dimensional vorticity equation.
        - `barotropicqg.jl`: Defines a `BarotropicQG` module that provides
                several solvers for the barotropic QG model that permit beta,
                topography, beta + topography, and forcing.
        - `traceradvdiff.jl`: Defines a `TracerAdvDiff` module that
                provides a solver for a two-dimensional and periodic tracer
                field in a given 2D flow (u, w), which can be an arbitrary
                function of x, z, and t.


## Basic Notation

The code solves partial differential equations of the general
form: `∂u/∂t = L*u + N(u)`, where `u` denotes the solution, which is
typically an array of complex coefficients for each Fourier mode. In general,
`L` is an array of coeffients that describe the linear part of the equation
governing `u`, while `N(u)` is an arbitrary function that may contain terms
both linear and nonlinear in `u`, as well forcing terms.
The time-steppers currently implemented only accepted
diagonal `L` arrays, which means that `L` and `u` have the same size.
Currently, the ETDRK4 time-stepper is the only time-stepper implemented that
makes special use of `L`. Both `L` and the function that calculates `N(u)` are
stored as fields of the type `AbstractEquation`.


## Writing fast solvers

The performance-intensive part of the code involves just two functions: the
timestepping scheme `stepforward!`, and the function `calcN!` that
calculates the nonlinear part of the given equation's right-hand side.
Optimization of these two functions for a given problem will produce the
fastest possible code.


## Future work

The code is in the chaotic stage of development. A main goal for the future
is to permit the use of shared memory parallelism in the compute-intensive
routines (shared-memory parallelism provided already by FFTW/MKLFFT, but
is not yet native to Julia for things like element-wise matrix multiplication,
addition, and assignment). This feature may possibly be enabled by
Intel Lab's [ParallelAccelerator][] package.


# Developers

FourierFlows is currently being developed by [Gregory L. Wagner][] (@glwagner)
and [Navid C. Constantinou][] (@navidcy).


[Julia]: https://julialang.org/
[ParallelAccelerator]: https://github.com/IntelLabs/ParallelAccelerator.jl
[Navid C. Constantinou]: http://www.navidconstantinou.com
[Gregory L. Wagner]: https://glwagner.github.io
