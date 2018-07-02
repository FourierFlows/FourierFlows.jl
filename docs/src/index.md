# FourierFlows.jl Documentation


## Overview

FourierFlows provides solvers for partial differential equations on doubly-periodic domains using Fourier-based pseudospectral methods. A central intent of the software's design is also to provide a framework for writing new, fast solvers for new physical problems. The code is written in [Julia](https://julialang.org).


## Usage

The code solves partial differential equations of the general form:

$\partial_t u = \mathcal{L}u + \mathcal{N}(u)\ .$

We decompose the right hand side of the above in a linear part ($\mathcal{L}u$) and a nonlinear part ($\mathcal{N}(u)$). The nonlinear part may include external forcing, e.g., $\mathcal{N}(u) = -u\partial_x u + f$.


## Installation

FourierFlows is a registered package and so can be installed via `Pkg.add`.

```julia
Pkg.add("FourierFlows")
```

For now, this package supports Julia `0.6`. Support for version `0.7` is on its way.


## Writing fast solvers

The performance-intensive part of the code involves just two functions: the time-stepping scheme `stepforward!`, and the function `calcN!` that calculates the nonlinear part of the given equation's right-hand side. Optimization of these two functions for a given problem will produce the fastest possible code.


## Future work

The code is in the chaotic stage of development. A main goal for the future is to permit the use of shared memory parallelism in the compute-intensive routines (shared-memory parallelism provided already by FFTW/MKLFFT, but is not yet native to Julia for things like element-wise matrix multiplication, addition, and assignment). This feature may possibly be enabled by Intel Lab's [ParallelAccelerator](https://github.com/IntelLabs/ParallelAccelerator.jl) package.


## Developers

FourierFlows is currently being developed by [Gregory L. Wagner](https://glwagner.github.io) and [Navid C. Constantinou](http://www.navidconstantinou.com).


## Cite

The code is citable via [zenodo](https://zenodo.org). Please cite as:

> Gregory L. Wagner & Navid C. Constantinou. (2018). FourierFlows/FourierFlows.jl: FourierFlows v0.1.1 (Version v0.1.1). Zenodo.  [http://doi.org/10.5281/zenodo.1302136](http://doi.org/10.5281/zenodo.1302136)
