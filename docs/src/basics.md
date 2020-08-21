# Code Basics

## Notation

The code solves differential equations of the form

```math
 \partial_t u = \mathcal{L}u + \mathcal{N}(u)\ ,
```
using Fourier transforms on periodic domains. Above, $u(\bm{x}, t)$ is the state 
variable. On the right-hand-side, term $\mathcal{L}u$ is a 'linear' part of the 
equation. The term $\mathcal{N}(u)$ is, in general, a 'nonlinear' part.

In FourierFlows.jl, $\mathcal{L}u$ is specified by physical modules as an array 
with the same dimension as $u$. The nonlinear term $\mathcal{N}(u)$ is specified 
via a function.

Boundary conditions in all spatial dimensions are periodic. That allows us to 
expand all variables using a Fourier decomposition. For example, if $u$ depends 
only in one spatial dimension and it's defined over the domain $x\in[0, L_x]$, 
then it can be written as:

```math
u(x, t) = \sum_{k} \hat{u}(k, t)\,\mathrm{e}^{\mathrm{i} k x}\ ,
```

where wavenumbers $k$ take the values $\tfrac{2\pi}{L_x}[0,\pm 1,\pm 2,\dots]$.

Equations are oftentimes time-stepped forward in Fourier space. That way $u$ 
becomes the array with all Fourier coefficients of the solution. Although this is
by any means restictive, it usually enhances performance.