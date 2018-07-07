# TwoDTurb Module

```math
\newcommand{\J}{\mathsf{J}}
```

### Basic Equations

This module solves two-dimensional incompressible turbulence. The flow is given
through a streamfunction $\psi$ as $(u,\upsilon) = (-\partial_y\psi, \partial_x\psi)$.
The dynamical variable used here is the component of the vorticity of the flow
normal to the plane of motion, $q=\partial_x \upsilon- \partial_y u = \nabla^2\psi$.
The equation solved by the module is:

$$\partial_t q + \J(\psi, q) = \underbrace{-\left[\mu(-1)^{n_\mu} \nabla^{2n_\mu}
+\nu(-1)^{n_\nu} \nabla^{2n_\nu}\right] q}_{\textrm{dissipation}} + f\ .$$

where $\J(a, b) = (\partial_x a)(\partial_y b)-(\partial_y a)(\partial_x b)$. On
the right hand side, $f(x,y,t)$ is forcing, $\mu$ is hypoviscosity, and $\nu$ is
hyperviscosity. Plain old linear drag corresponds to $n_{\mu}=0$, while normal
viscosity corresponds to $n_{\nu}=1$.

### Implementation

The equation is time-stepped forward in Fourier space:

$$\partial_t \widehat{q} = - \widehat{J(\psi, q)} -\left(\mu k^{2n_\mu}
+\nu k^{2n_\nu}\right) \widehat{q}  + \widehat{f}\ .$$

In doing so the Jacobian is computed in the conservative form: $\J(a,b) =
\partial_y [ (\partial_x a) b] -\partial_x[ (\partial_y a) b]$.

Thus:

$$\mathcal{L} = -\mu k^{-2n_\mu} - \nu k^{2n_\nu}\ ,$$
$$\mathcal{N}(\widehat{q}) = - \mathrm{i}k_x \mathrm{FFT}(u q)-
	\mathrm{i}k_y \mathrm{FFT}(\upsilon q) + \widehat{f}\ .$$


### AbstractTypes and Functions

**Params**

For the unforced case ($f=0$) parameters AbstractType is build with `Params` and it includes:
- `nu`:   Float; viscosity or hyperviscosity coefficient.
- `nnu`: Integer$>0$; the order of viscosity $n_\nu$. Case $n_\nu=1$ give normal viscosity.
- `mu`: Float; bottom drag or hypoviscosity coefficient.
- `nmu`: Integer$\ge 0$; the order of hypodrag $n_\mu$. Case $n_\mu=0$ give plain linear drag $\mu$.

For the forced case ($f\ne 0$) parameters AbstractType is build with `ForcedParams`. It includes all parameters in `Params` and additionally:
- `calcF!`: Function that calculates the forcing $\widehat{f}$


**Vars**

For the unforced case ($f=0$) variables AbstractType is build with `Vars` and it includes:
- `q`: Array of Floats; relative vorticity.
- `U`: Array of Floats; $x$-velocity, $u$.
- `V`: Array of Floats; $y$-velocity, $v$.
- `sol`: Array of Complex; the solution, $\widehat{q}$.
- `qh`: Array of Complex; the Fourier transform $\widehat{q}$.
- `Uh`: Array of Complex; the Fourier transform $\widehat{u}$.
- `Vh`: Array of Complex; the Fourier transform $\widehat{v}$.

For the forced case ($f\ne 0$) variables AbstractType is build with `ForcedVars`. It includes all variables in `Vars` and additionally:
- `Fh`: Array of Complex; the Fourier transform $\widehat{f}$.
- `prevsol`: Array of Complex; the values of the solution `sol` at the previous time-step (useful for calculating the work done by the forcing).



**`calcN!` function**

The nonlinear term $\mathcal{N}(\widehat{q})$ is computed via functions:

- `calcN_advection!`: computes $- \widehat{J(\psi, q)}$ and stores it in array `N`.

```julia
function calcN_advection!(N, sol, t, s, v, p, g)
  @. v.Uh =  im * g.l  * g.invKKrsq * sol
  @. v.Vh = -im * g.kr * g.invKKrsq * sol
  @. v.qh = sol

  A_mul_B!(v.U, g.irfftplan, v.Uh)
  A_mul_B!s(v.V, g.irfftplan, v.Vh)
  A_mul_B!(v.q, g.irfftplan, v.qh)

  @. v.U *= v.q # U*q
  @. v.V *= v.q # V*q

  A_mul_B!(v.Uh, g.rfftplan, v.U) # \hat{U*q}
  A_mul_B!(v.Vh, g.rfftplan, v.V) # \hat{U*q}

  @. N = -im*g.kr*v.Uh - im*g.l*v.Vh
  nothing
end
```

- `calcN_forced!`: computes $- \widehat{J(\psi, q)}$ via `calcN_advection!` and then adds to it the forcing $\widehat{f}$ computed via `calcF!` function. Also saves the solution $\widehat{q}$ of the previous time-step in array `prevsol`.

```julia
function calcN_forced!(N, sol, t, s, v, p, g)
  calcN_advection!(N, sol, t, s, v, p, g)
  if t == s.t # not a substep
    v.prevsol .= s.sol # used to compute budgets when forcing is stochastic
    p.calcF!(v.Fh, sol, t, s, v, p, g)
  end
  @. N += v.Fh
  nothing
end
```
- `updatevars!`: uses `sol` to compute $q$, $u$, $v$, $\widehat{u}$, and $\widehat{v}$ and stores them into corresponding arrays of `Vars`/`ForcedVars`.


## Examples

- `examples/twodturb/McWilliams.jl`: A script that simulates decaying two-dimensional turbulence reproducing the results of the paper by

  > McWilliams, J. C. (1984). The emergence of isolated coherent vortices in turbulent flow. *J. Fluid Mech.*, **146**, 21-43.

- `examples/twodturb/IsotropicRingForcing.jl`: A script that simulates stochastically forced two-dimensional turbulence. The forcing is temporally delta-corraleted and its spatial structure is isotropic with power in a narrow annulus of total radius $k_f$ in wavenumber space.
