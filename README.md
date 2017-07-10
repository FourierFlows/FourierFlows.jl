# FourierFlows

## Overview

This software provides solvers for partial differential equations on
doubly-periodic domains using Fourier-based pseudospectral methods.
A central intent of the software's design is also to provide a framework
for writing new, fast solvers for new phyiscal problems. Achieving this goal
requires a simultaneous capacity for accepting 'slow' but high-level code, as well 
as lower-level, optimized routines. The code is written in the [Julia][].

## Code organization

The code is divided along conceptual lines into problem-agnostic and 
problem-specific components. The problem-agnostic parts of the code are 
contained in ``/src``. Files in ``/src`` define the domain, master 'AbstractTypes'
for the model, time-stepper types and time-stepping routines.

The problem-specific modules, which provide the problem-specific ``Vars`` and ``Params``
types, as well as routines for calculating linear and nonlinear parts of the 
specified equation, are in the directory ``/src/physics``.

## Future work 

The code is an a very early stage of development. A main goal for the future
is to permit the use of shared memory parallelism in the compute-intensive 
routines (shared-memory parallelism provided already by FFTW/MKLFFT, but 
is not yet native to Julia for things like element-wise matrix multiplication, 
addition, and assignment). This feature may possibly be enabled by 
Intel Lab's [ParallelAccelerator][] package.

# Authors

Fourier flows is currently being developed by [Gregory L. Wagner] (@glwagner) 
and [Navid Constantinou][] (@navidcy)


[Julia]: https://julialang.org/
[ParallelAccelerator]: https://github.com/IntelLabs/ParallelAccelerator.jl
[Navid Constantinou]: http://www.navidconstantinou.com/
[Gregory L. Wagner]: https://glwagner.github.io
