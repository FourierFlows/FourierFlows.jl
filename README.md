# FourierFlows.jl

<p align="left">
    <a href="https://travis-ci.com/FourierFlows/FourierFlows.jl">
        <img alt="Build Status for CPU" src="https://img.shields.io/travis/com/FourierFlows/FourierFlows.jl/master?label=CPU&logo=travis&logoColor=white">
    </a>
    <a href="https://gitlab.com/JuliaGPU/FourierFlows.jl/commits/master">
      <img alt="Build Status for GPU" src="https://img.shields.io/gitlab/pipeline/JuliaGPU/FourierFlows.jl/master?label=GPU&logo=gitlab&logoColor=white">
    </a>
    <a href="https://ci.appveyor.com/project/navidcy/fourierflows-jl">
        <img alt="Build Status for Windows" src="https://img.shields.io/appveyor/ci/navidcy/fourierflows-jl/master?label=Window&logo=appveyor&logoColor=white">
    </a>
    <a href="https://FourierFlows.github.io/FourierFlows.jl/stable">
        <img alt="stable docs" src="https://img.shields.io/badge/docs-stable-blue.svg">
    </a>
    <a href="https://FourierFlows.github.io/FourierFlows.jl/dev">
        <img alt="latest docs" src="https://img.shields.io/badge/docs-dev-blue.svg">
    </a>
    <a href="https://coveralls.io/github/FourierFlows/FourierFlows.jl?branch=master">
        <img src="https://coveralls.io/repos/github/FourierFlows/FourierFlows.jl/badge.svg?branch=master" title="coveralls">
    </a>
    <a href="https://codecov.io/gh/FourierFlows/FourierFlows.jl">
        <img src="https://codecov.io/gh/FourierFlows/FourierFlows.jl/branch/master/graph/badge.svg" title="codecov">
    </a>
    <a href="https://doi.org/10.5281/zenodo.1161724">
        <img alt="doi" src="https://zenodo.org/badge/DOI/10.5281/zenodo.1161724.svg" alt="DOI">
    </a>
</p>

## Overview

This software provides tools for partial differential equations on
periodic domains using Fourier-based pseudospectral methods.
A central intent of the software's design is also to provide a framework
for writing new, fast solvers for new physical problems.
The code is written in [Julia][].

For more details refer to the [documentation](https://fourierflows.github.io/FourierFlows.jl/latest/).

## Installation

But it is simple:

```julia
] add FourierFlows
```

and no more.

## Example(s)

See the child package [GeophysicalFlows.jl][] for example usage of `FourierFlows.jl` for problems in 
Geophysical Fluid Dynamics.

## Developers

FourierFlows is currently being developed by [Gregory L. Wagner][] (@glwagner)
and [Navid C. Constantinou][] (@navidcy).


## Cite

The code is citable via [zenodo](https://zenodo.org). Please cite as:

> Navid C. Constantinou & Gregory L. Wagner. (2019). FourierFlows/FourierFlows.jl: FourierFlows v0.3.2 (Version v0.3.2). Zenodo.  [http://doi.org/10.5281/zenodo.1161724](http://doi.org/10.5281/zenodo.1161724)



[Julia]: https://julialang.org/
[Navid C. Constantinou]: http://www.navidconstantinou.com
[Gregory L. Wagner]: https://glwagner.github.io
[GeophysicalFlows.jl]: https://github.com/FourierFlows/GeophysicalFlows.jl
