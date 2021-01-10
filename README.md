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
    <a href="https://FourierFlows.github.io/FourierFlowsDocumentation/stable">
        <img alt="stable docs" src="https://img.shields.io/badge/documentation-stable%20release-blue">
    </a>
    <a href="https://FourierFlows.github.io/FourierFlowsDocumentation/dev">
        <img alt="latest docs" src="https://img.shields.io/badge/documentation-in%20development-orange">
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

For more details refer to the [documentation](https://fourierflows.github.io/FourierFlowsDocumentation/dev/).

## Installation

It is simple:

```julia
] add FourierFlows
```

and no more.

## Usage

See the documentation for tutorials on (i) how [construct grids](https://fourierflows.github.io/FourierFlowsDocumentation/stable/grids/) and use Fourier transform to compute derivatives and (ii) how to [set up a PDE](https://fourierflows.github.io/FourierFlowsDocumentation/stable/problem/), time-step it forward. and visualize the output.

## Example(s)

An example for coding up and solving the linear shallow water equations is [documented](https://fourierflows.github.io/FourierFlowsDocumentation/stable/generated/OneDShallowWaterGeostrophicAdjustment/).

See also the child packages [GeophysicalFlows.jl][] for example usage of `FourierFlows.jl` for problems in Geophysical Fluid Dynamics.

## Developers

FourierFlows is currently being developed by [Gregory L. Wagner][] (@glwagner)
and [Navid C. Constantinou][] (@navidcy).


## Cite

The code is citable via [zenodo](https://zenodo.org). Please cite as:

> Navid C. Constantinou & Gregory L. Wagner. (2021). FourierFlows/FourierFlows.jl: FourierFlows v0.6.9 (Version v0.6.9). Zenodo.  [http://doi.org/10.5281/zenodo.1161724](http://doi.org/10.5281/zenodo.1161724)


[Julia]: https://julialang.org/
[Navid C. Constantinou]: http://www.navidconstantinou.com
[Gregory L. Wagner]: https://glwagner.github.io
[GeophysicalFlows.jl]: https://github.com/FourierFlows/GeophysicalFlows.jl
