# Code Basics

## Notation

The code solves partial differential equations of the form

```math
 \partial_t u = Lu + N(u)\ ,
```
using Fourier transforms on periodic domains. Above, ``u(\bm{x}, t)`` is the state variable. 
On the right-hand-side, term ``L u`` is the 'linear' part of the equation. The term 
``N(u)`` is, in general, a 'nonlinear' part.

In FourierFlows.jl, ``L u`` is specified by the various modules by prescribing the
linear operator ``L`` as an array of the same dimension as ``u``. The nonlinear term 
``N(u)`` is specified via a function that takes ``u`` as its argument.

Boundary conditions in all spatial dimensions are periodic. That allows us to expand all 
variables using a Fourier decomposition. For example, if ``u`` depends only in one spatial 
dimension, ``x``, defined over the domain ``x \in [0, L_x]``, then:

```math
u(x, t) = \sum_{k_x} \hat{u}(k_x, t) \, e^{i k_x x} \ ,
```

where wavenumbers ``k_x`` take the values ``\tfrac{2\pi}{L_x}\{0,\pm 1,\pm 2,\dots\}``. When we 
further consider that ``x`` takes discrete values over ``[0, L_x]``, e.g., ``x_j``,
``j = 0, 1, \dots, n_x``, then only ``n_x`` wavenumbers are indepented. If we denote ``u_j(t) \equiv u(x_j, t)`` and ``\hat{u}_k(t) \equiv \hat{u}(\tfrac{2\pi}{L_x} k, t)`` with ``k`` integer, then the 
Fourier sum above truncates to:

```math
u_j(t) = \sum_{k=-n_x/2}^{n_x/2-1} \hat{u}_k(t)\,e^{2\pi i k x_j / L_x}\ .
```

(Above we assumed that ``n_x`` is even.)

Applying the Fourier transform as above, the partial differential equation transforms to:

```math
 \partial_t \hat{u} = \hat{L} \hat{u} + \widehat{ N(u) }\ ,
```

where ``\hat{L}`` above denotes the linear operator in Fourier space.

Equations are oftentimes time-stepped forward in Fourier space. Doing so, ``\hat{u}_k(t)`` 
becomes our state variable, i.e., the array with all Fourier coefficients of the solution 
``u``. Although time-stepping in Fourier space is by no means a restriction of the code, it 
usually enhances performance.