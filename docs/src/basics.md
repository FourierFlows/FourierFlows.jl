# Code Basics

## Notation

The code solves partial differential equations of the form

```math
 \partial_t u = \mathcal{L}u + \mathcal{N}(u)\ ,
```
using Fourier transforms on periodic domains. Above, ``u(\bm{x}, t)`` is the state variable. 
On the right-hand-side, term ``\mathcal{L} u`` is a 'linear' part of the equation. The term 
``\mathcal{N}(u)`` is, in general, a 'nonlinear' part.

In FourierFlows.jl, ``\mathcal{L} u`` is specified by the various modules by prescribing the
linear operator ``\mathcal{L}`` as an array of the same dimension as ``u``. The nonlinear term 
``\mathcal{N}(u)`` is specified via a function that takes ``u`` as its argument.

Boundary conditions in all spatial dimensions are periodic. That allows us to expand all 
variables using a Fourier decomposition. For example, if ``u`` depends only in one spatial 
dimension, ``x``, defined in the domain ``x \in [0, L_x]``, then:

```math
u(x, t) = \sum_{k} \hat{u}(k, t)\,\mathrm{e}^{\mathrm{i} k x}\ ,
```

where wavenumbers ``k`` take the values ``\tfrac{2\pi}{L_x}\{0,\pm 1,\pm 2,\dots\}``. In this 
way, the differential equation transforms to:

```math
 \partial_t \hat{u} = \hat{\mathcal{L}} \hat{u} + \widehat{ \mathcal{N}(u) }\ ,
```

where ``\hat{\mathcal{L}}`` above denotes the linear operator in Fourier space.

Equations are oftentimes time-stepped forward in Fourier space. Doing so, ``\hat{u}(k, t)`` 
becomes our state variable, i.e., the array with all Fourier coefficients of the solution 
``u``. Although time-stepping in Fourier space is by no means a restriction of the code, it 
usually enhances performance.