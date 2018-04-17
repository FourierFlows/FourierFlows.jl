# TracerAdvDiff Module

This module solves the advection diffusion equation for a passive tracer
concentration $c(x,y,t)$ in two-dimensions:

$$\partial_t c + \boldsymbol{u} \boldsymbol{\cdot} \boldsymbol{\nabla} c = \kappa \nabla^2 c\ ,$$

where $\boldsymbol{u} = (u,\upsilon)$ is the two-dimensional advecting velocity and $\kappa$ is the diffusivity.
