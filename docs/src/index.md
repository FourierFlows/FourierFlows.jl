# FourierFlows.jl Documentation

## Overview

FourierFlows provides a framework to write solvers for partial differential equations on periodic domains with
Fourier-based pseudospectral methods that run seamlessly on CPUs and GPUs. We support 1-, 2-, and 3-dimensional domains.

## Examples

An example demonstrating how to code up and solve the linear shallow water equations is found
in the [Examples](https://fourierflows.github.io/FourierFlowsDocumentation/stable/generated/OneDShallowWaterGeostrophicAdjustment/) section of the documentation.

!!! warning "Animations and Mac OS X Safari"
    Often animations included the documentation, e.g., those in the [linear shallow water example](https://fourierflows.github.io/FourierFlowsDocumentation/stable/generated/OneDShallowWaterGeostrophicAdjustment/)
    don't show up in Safari. Until this issue is fixed we suggest using alternative browsers 
    to view the documentation.

For more examples of `FourierFlows.jl` in action, see the child packages
[`GeophysicalFlows.jl`](https://github.com/FourierFlows/GeophysicalFlows.jl) or [`PassiveTracerFlows.jl`](https://github.com/FourierFlows/PassiveTracerFlows.jl).

## Developers

FourierFlows is currently being developed by [Gregory L. Wagner](https://glwagner.github.io) and 
[Navid C. Constantinou](http://www.navidconstantinou.com).

## Cite

The code is citable via [zenodo](https://doi.org/10.5281/zenodo.1161724).
