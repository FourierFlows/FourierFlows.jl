# FourierFlows

## Overview

This software provides solvers for partial differential equations on
doubly-periodic domains using Fourier-based pseudospectral methods.
A central intent of the software's design is also to provide a framework 
for writing new, fast solvers for new physical problems. Achieving this goal
requires a simultaneous capacity for accepting 'slow' but high-level code, as well 
as lower-level, optimized routines. The code is written in the [Julia][].

## Tests

Coming soon...

## Examples

Coming soon...

## Source code organization

The code is divided along conceptual lines into problem-agnostic and
problem-specific components. Files that contain problem-agnostic parts
of the code are stored in ``/src``. Files in ``/src`` define the domain,
'AbstractTypes' that supertype problem-specific types, and
time-stepper types and routines. Problem-specific modules are stores in
``/src/physics``.

Here's an overview of the code structure:

- ``/src/``
    - ``fourierflows.jl``
        - Module ``FourierFlowTypes``: *defines supertyping AbstractTypes.*
        - Module ``TwoDGrid``: *defines a composite type for 2D physical and Fourier
            spectral grids.*
        - Modules contained in ``timesteppers.jl`` via inclusion at the end of the file.
   - ``timesteppers.jl``: *defines modules and ``stepforward!`` routines for various
        time-steppers. Current implemented time-steppers are:*
        - Forward Euler
        - 3rd-order Adams-Bashforth (AB3)
        - 4th-order Runge-Kutta (RK4)
        - 4th-order Runge-Kutta Exponential Time Differencing (ETDRK4)
    - ``physics/``
        - ``twodturb.jl`` *Defines a ``TwoDTurb`` module that provides a solver for the
                two-dimensional vorticity equation.*
        - ``barotropicqg.jl`` *Defines a ``BarotropicQG`` module that provides several
                solvers for the barotropic QG model that permit beta, topography,
                beta + topography, and forcing.*


## Writing fast solvers

The performance-intensive part of the code involves just two functions: the timestepping
scheme ``stepforward!``, and the function ``calcNL!`` that calculates the nonlinear part
of the given equation's right-hand side. Optimization of these two functions for a given
problem will produce the fastest possible code.

## Future work

The code is an a very early stage of development. A main goal for the future
is to permit the use of shared memory parallelism in the compute-intensive
routines (shared-memory parallelism provided already by FFTW/MKLFFT, but
is not yet native to Julia for things like element-wise matrix multiplication,
addition, and assignment). This feature may possibly be enabled by
Intel Lab's [ParallelAccelerator][] package.

# Authors

Fourier flows is currently being developed by [Gregory L. Wagner][] (@glwagner)
and [Navid C. Constantinou][] (@navidcy)


[Julia]: https://julialang.org/
[ParallelAccelerator]: https://github.com/IntelLabs/ParallelAccelerator.jl
[Navid C. Constantinou]: http://www.navidconstantinou.com
[Gregory L. Wagner]: https://glwagner.github.io
