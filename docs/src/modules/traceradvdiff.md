# TracerAdvDiff Module

### Basic Equations

This module solves the advection diffusion equation for a passive tracer
concentration $c(x, y, t)$ in two-dimensions:

$$\partial_t c + \boldsymbol{u} \boldsymbol{\cdot} \boldsymbol{\nabla} c = \underbrace{\eta \partial_x^2 c + \kappa \partial_y^2 c}_{\textrm{diffusivity}} + \underbrace{\kappa_h (-1)^{n_{h}} \nabla^{2n_{h}}c}_{\textrm{hyper-diffusivity}}\ ,$$

where $\boldsymbol{u} = (u,v)$ is the two-dimensional advecting flow, $\eta$ the $x$-diffusivity and $\kappa$ is the $y$-diffusivity. If $\eta$ is not defined then the code uses isotropic diffusivity, i.e., $\eta \partial_x^2 c + \kappa \partial_y^2 c\mapsto\kappa\nabla^2$. The advecting flow could be either compressible or incompressible. 


### Implementation

The equation is time-stepped forward in Fourier space:

$$\partial_t \widehat{c} = - \widehat{\boldsymbol{u} \boldsymbol{\cdot} \boldsymbol{\nabla} c} - \left[ (\eta k_x^2 + \kappa k_y^2) +\kappa_h k^{2\nu_h} \right]\widehat{c}\ .$$

Thus:

```math
\begin{align*}
\mathcal{L} &= -\eta k_x^2 - \kappa k_y^2 - \kappa_h k^{2\nu_h}\ , \\
\mathcal{N}(\widehat{c}) &= - \mathrm{FFT}(u \partial_x c + \upsilon \partial_y c)\ .
\end{align*}
```
