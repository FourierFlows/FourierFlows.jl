# FourierFlows.jl

<p align="left">
    <a href="https://buildkite.com/julialang/fourierflows-dot-jl">
        <img alt="Buildkite CPU+GPU build status" src="https://img.shields.io/buildkite/4d921fc17b95341ea5477fb62df0e6d9364b61b154e050a123/master?logo=buildkite&label=Buildkite%20CPU%2BGPU">
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
        <img src="https://codecov.io/gh/FourierFlows/FourierFlows.jl/branch/main/graph/badge.svg?token=BrgeSmKJHD"/>
    </a>
    <a href="https://github.com/SciML/ColPrac">
      <img alt="ColPrac: Contributor's Guide on Collaborative Practices for Community Packages" src="https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet">
    </a>
    <a href="https://doi.org/10.5281/zenodo.1161724">
        <img alt="doi" src="https://zenodo.org/badge/DOI/10.5281/zenodo.1161724.svg" alt="DOI">
    </a>
</p>


## Overview

This software provides tools for partial differential equations on periodic domains using 
Fourier-based pseudospectral methods. A central intent of the software's design is also to 
provide a framework for writing new, fast solvers for new physical problems. The code is 
written in [Julia][].

For more details refer to the [documentation](https://fourierflows.github.io/FourierFlowsDocumentation/stable/).


## Installation

To install, use Julia's built-in package manager to add the package and also to instantiate/build all the required dependencies

```julia
julia> using Pkg
julia> Pkg.add("FourierFlows")
julia> Pkg.instantiate()
```

The most recent version of FourierFlows.jl requires Julia v1.6 (the current long-term-release) or later.

Last version that works with Julia v1.5 is FourierFlows.jl v0.7.2.


## Usage

See the documentation for tutorials on (i) how [construct grids](https://fourierflows.github.io/FourierFlowsDocumentation/stable/grids/) and use Fourier transform to compute derivatives and (ii) how to [set up a PDE](https://fourierflows.github.io/FourierFlowsDocumentation/stable/problem/), time-step it forward, and visualize the output, (iii) how to [add diagnostics](https://fourierflows.github.io/FourierFlowsDocumentation/stable/diagnostics/), and (iv) how to [write and read output](https://fourierflows.github.io/FourierFlowsDocumentation/stable/output/) from disk.


## Scalability

For now, FourierFlows.jl is restricted to run on either a single CPU or single GPU. Multi-threading
can enhance performance for the Fourier transforms. By default, FourierFlows.jl will use the 
maximum number of threads available on your machine. You can set the number of threads used by
FourierFlows.jl by setting the environment variable, e.g.,

```
$ export JULIA_NUM_THREADS=4
```

For more information on multi-threading users are directed to the [Julia Documentation](https://docs.julialang.org/en/v1/manual/multi-threading/).

If your machine has more than one GPU available, then functionality within CUDA.jl package 
enables the user to choose the GPU device that FourierFlows.jl should use. The user is referred
to the [CUDA.jl Documentation](https://juliagpu.github.io/CUDA.jl/stable/lib/driver/#Device-Management);
in particular, [`CUDA.devices`](https://juliagpu.github.io/CUDA.jl/stable/lib/driver/#CUDA.devices) 
and [`CUDA.CuDevice`](https://juliagpu.github.io/CUDA.jl/stable/lib/driver/#CUDA.CuDevice).


## Example(s)

An example for coding up and solving the linear shallow water equations is [documented](https://fourierflows.github.io/FourierFlowsDocumentation/stable/literated/OneDShallowWaterGeostrophicAdjustment/).

See also the child packages [GeophysicalFlows.jl][] for example usage of FourierFlows.jl for 
problems in Geophysical Fluid Dynamics.


## Getting help

Interested in using FourierFlows.jl or trying to figure out how to use it? Please feel free 
to ask us questions and get in touch! The [documentation](https://fourierflows.github.io/FourierFlowsDocumentation/stable) 
is always the best place to start. Check out the [examples](https://github.com/FourierFlows/FourierFlows.jl/tree/master/examples) and [open an issue](https://github.com/FourierFlows/FourierFlows.jl/issues/new) 
or [start a discussion](https://github.com/FourierFlows/FourierFlows.jl/discussions/new) 
if you have any questions, comments, suggestions, etc.


## Developers

FourierFlows is currently being developed by [Gregory L. Wagner][] (@glwagner)
and [Navid C. Constantinou][] (@navidcy).


## Citing

The code is citable via [zenodo](https://zenodo.org). Please cite as:

> Gregory L. Wagner, Navid C. Constantinou, and contributors. (2022). FourierFlows/FourierFlows.jl: FourierFlows v0.9.4 (Version v0.9.4). Zenodo. [http://doi.org/10.5281/zenodo.1161724](http://doi.org/10.5281/zenodo.1161724)


## Contributing

We are excited to get more people involved in contributing to the development of FourierFlows.jl! 
We welcome any contribution, no matter how big or small! It's always great to have new people 
look at the code with fresh eyes: you will see errors that other developers have missed.

Let us know by [open an issue](https://github.com/FourierFlows/FourierFlows.jl/issues/new) 
or [start a discussion](https://github.com/FourierFlows/FourierFlows.jl/discussions/new) 
if you'd like to work on a new feature or implement a new module, if you're new to open-source 
and want to find a cool little project or issue to work on that fits your interests! We're more 
than happy to help along the way.

For more information, check out our [contributor's guide](https://github.com/FourierFlows/FourierFlows.jl/blob/master/CONTRIBUTING.md).


[Julia]: https://julialang.org/
[Navid C. Constantinou]: http://www.navidconstantinou.com
[Gregory L. Wagner]: https://glwagner.github.io
[GeophysicalFlows.jl]: https://github.com/FourierFlows/GeophysicalFlows.jl
