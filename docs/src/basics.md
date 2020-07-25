# Code Basics

## Notation

The code solves differential equations of the form

```math
 \partial_t u = \mathcal{L}u + \mathcal{N}(u)\ ,
```
using Fourier transforms on periodic domains. The right side term $\mathcal{L}u$ is a 'linear' part of the equation.
The term $\mathcal{N}(u)$ is, in general, a 'nonlinear' part. In FourierFlows, $\mathcal{L}u$ is specified by 
physical modules in Fourier space as an array with the same size as $\hat u$, the Fourier transform of $u$.
The nonlinear term $\mathcal{N}u$ is specified by a function.

Boundary conditions in all spatial dimensions are periodic. That allows us to expand all variables using a Fourier 
decomposition. For example, a variable $\phi(x, t)$ that depends in one spatial dimension is expanded as:


```math
\phi(x, t) = \sum_{k} \widehat{\phi}(k, t)\,e^{\mathrm{i} k x}\ ,
```
where wavenumbers $k$ take the values $\tfrac{2\pi}{L_x}[0,\pm 1,\pm 2,\dots]$. The equation is time-stepped 
forward in Fourier space. That way $u$ becomes the array with all Fourier coefficients of the solution.

The coefficients for the linear operator $\mathcal{L}$ are stored in an array called `L`. The term $\mathcal{N}(u)$ 
is computed through the function `calcN!`.
