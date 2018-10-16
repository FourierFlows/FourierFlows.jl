# FourierFlows.jl

<p align="left">
    <a href="https://travis-ci.org/FourierFlows/FourierFlows.jl">
        <img src="https://travis-ci.org/FourierFlows/FourierFlows.jl.svg?branch=master" title="Build Status">
    </a>
    <a href="https://ci.appveyor.com/project/navidcy/fourierflows-jl">
        <img src="https://ci.appveyor.com/api/projects/status/3hm86k8d4qdch730?svg=true" title="Build Status">
    </a>
    <a href="https://FourierFlows.github.io/FourierFlows.jl/stable">
        <img src="https://img.shields.io/badge/docs-stable-blue.svg">
    </a>
    <a href="https://FourierFlows.github.io/FourierFlows.jl/latest">
        <img src="https://img.shields.io/badge/docs-latest-blue.svg">
    </a>
    <a href="http://pkg.julialang.org/detail/FourierFlows">
        <img src="http://pkg.julialang.org/badges/FourierFlows_0.7.svg" title="FourierFlows">
    </a>
    <a href="https://codecov.io/gh/FourierFlows/FourierFlows.jl">
        <img src="https://codecov.io/gh/FourierFlows/FourierFlows.jl/branch/master/graph/badge.svg" title="codecov">
    </a>
    <a href="https://doi.org/10.5281/zenodo.1161724">
        <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.1161724.svg" alt="DOI">
    </a>
</p>

## Overview

This software provides tools for partial differential equations on
doubly-periodic domains using Fourier-based pseudospectral methods.
A central intent of the software's design is also to provide a framework
for writing new, fast solvers for new physical problems.
The code is written in [Julia][].

For more details refer to the [documentation](https://fourierflows.github.io/FourierFlows.jl/latest/).

## Installation

The is a registered package and can be installed by:

```julia
julia> using Pkg
julia> Pkg.add("FourierFlows")
```

## Example

For a simple example involving the advection and diffusion of a passive tracer,
see `examples/tracers_cellularflow.jl`.

# Developers

FourierFlows is currently being developed by [Gregory L. Wagner][] (@glwagner)
and [Navid C. Constantinou][] (@navidcy).


## Cite

The code is citable via [zenodo](https://zenodo.org). Please cite as:

> Gregory L. Wagner & Navid C. Constantinou. (2018). FourierFlows/FourierFlows.jl: FourierFlows v0.2.0 (Version v0.2.0). Zenodo.  [http://doi.org/10.5281/zenodo.1161724](http://doi.org/10.5281/zenodo.1161724)



[Julia]: https://julialang.org/
[Navid C. Constantinou]: http://www.navidconstantinou.com
[Gregory L. Wagner]: https://glwagner.github.io
