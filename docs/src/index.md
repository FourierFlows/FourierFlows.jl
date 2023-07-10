# FourierFlows.jl Documentation

## Overview

FourierFlows provides a framework to write solvers for partial differential equations on
periodic domains with Fourier-based pseudospectral methods that run seamlessly on CPUs and GPUs.
We support 1-, 2-, and 3-dimensional domains.

The next few sections step through the basics of building and a module to solve a partial
differential equation.

## Examples

A demonstration for how to code up and solve the linear shallow water equations is found
in the [Examples](literated/OneDShallowWaterGeostrophicAdjustment/) section of the documentation.

For more examples of `FourierFlows.jl` in action, see the child packages
[`GeophysicalFlows.jl`](https://github.com/FourierFlows/GeophysicalFlows.jl)
or [`PassiveTracerFlows.jl`](https://github.com/FourierFlows/PassiveTracerFlows.jl).

!!! info "Unicode"
    Often unicode symbols appear in modules for variables or parameters. For example,
    `κ` appears as the diffusivity in the `Diffusion` module. Unicode symbols can be entered 
    in the Julia REPL by typing, e.g., `\kappa` followed by `tab` key. Read more about Unicode 
    symbols in the [Julia Documentation](https://docs.julialang.org/en/v1/manual/unicode-input/).

## Developers

FourierFlows.jl started during the on Atmospheric and Oceanic Fluid Dynamics Meeting 2017 by
[Greg Wagner](https://glwagner.github.io) and [Navid Constantinou](https://www.navidconstantinou.com).
Since then [various people have contributed](https://github.com/FourierFlows/FourierFlows.jl/graphs/contributors).

## Cite

The software is citable via [zenodo](https://doi.org/10.5281/zenodo.1161724).
